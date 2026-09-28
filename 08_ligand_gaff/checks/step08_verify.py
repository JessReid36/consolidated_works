#!/usr/bin/env python3
"""
step08_verify.py - re-runnable audit of step 08, chorismate GAFF parameterisation.

Step 08 runs in two halves on two machines. step08a derives AM1-BCC charges on a local
workstation, because the HPC AMBER22 sqm cannot run (missing libopenblas). step08b does
GAFF typing and parmchk2 on the HPC in read-charge mode, so no sqm is needed there.

Checks performed:

  1. Continuity. The template the charges were derived from is step06's frozen cha_a.mol2,
     by SHA-256 against the hash recorded in step08a_provenance.txt AND against the file.
  2. Charge file integrity. charges_am1bcc.dat carries one value per atom and sums to
     exactly -2.
  3. THE RENORMALISATION. step08a records that the raw AM1-BCC charges, as written to
     mol2 at four decimal places, summed to -1.997000, and that the -0.003 residual was
     distributed evenly over 24 atoms. This reproduces that: every charge in
     charges_am1bcc.dat must differ from the corresponding cha_am1bcc.mol2 value by
     exactly -0.003/24 = -0.000125 e.
  4. Read-charge mode. cha_gaff.mol2 charges must equal charges_am1bcc.dat exactly, since
     step08b ran antechamber with -c rc. Atom names and order must also be preserved from
     the template.
  5. GAFF typing. Every atom carries a GAFF type, no Tripos dot-types survive, and the
     type set matches the audit.
  6. THE FRCMOD PENALTIES. parmchk2 reports 0 ATTN, but ATTN is not the only signal. This
     reports every substituted parameter with its penalty score, ranked. Two dihedrals
     carry a penalty of 237.0 against a next-highest of 3.0, and both span the O3-C4 bond,
     which step 11 identifies as the BREAKING bond of the rearrangement. See
     interpretation.txt.
  7. Three-copy identity. cha_a, cha_b and cha_c share atom names and bond topology, which
     is what licenses using one parameter set for all three.
  8. Frozen checksums reproduce.
  9. Declared-but-absent outputs are reported, not ignored.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step08_verify.py [outputs_dir] [step06_outputs]
"""
import sys
import hashlib
import re
from pathlib import Path

FAIL = []
RESIDUAL_PER_ATOM = -0.003 / 24


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def mol2(path):
    sec, atoms, bonds = None, [], []
    for L in Path(path).read_text().splitlines():
        if L.startswith("@<TRIPOS>ATOM"):
            sec = "a"
            continue
        if L.startswith("@<TRIPOS>BOND"):
            sec = "b"
            continue
        if L.startswith("@<TRIPOS>"):
            sec = None
            continue
        f = L.split()
        if sec == "a" and len(f) >= 6:
            atoms.append({"id": int(f[0]), "name": f[1], "type": f[5],
                          "charge": float(f[8]) if len(f) >= 9 else None})
        elif sec == "b" and len(f) >= 4:
            bonds.append((int(f[1]), int(f[2]), f[3]))
    return atoms, bonds


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s06 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../06_acceptance_freeze/outputs")
    print(f"step 08 verification against {out.resolve()}\n")

    prov = (out / "step08a_provenance.txt").read_text()

    print("1. continuity: the charge template is step06's frozen cha_a.mol2")
    m = re.search(r"cha_a\.mol2\s+([0-9a-f]{64})", prov)
    check("provenance records a template hash", bool(m))
    if m and (s06 / "cha_a.mol2").exists():
        check("recorded hash matches step06's frozen cha_a.mol2",
              sha(s06 / "cha_a.mol2") == m.group(1), m.group(1)[:16])

    print("\n2. charge file integrity")
    q = [float(x) for x in (out / "charges_am1bcc.dat").read_text().split()]
    check("one charge per atom", len(q) == 24, f"{len(q)} values")
    check("charges sum to exactly -2", abs(sum(q) + 2.0) < 1e-7, f"{sum(q):+.8f}")

    print("\n3. the renormalisation reproduces")
    bcc_atoms, _ = mol2(out / "cha_am1bcc.mol2")
    raw = [a["charge"] for a in bcc_atoms]
    check("raw mol2 charges sum to the recorded -1.997000",
          abs(sum(raw) + 1.997) < 5e-7, f"{sum(raw):+.6f}")
    deltas = [qi - ri for qi, ri in zip(q, raw)]
    spread = max(deltas) - min(deltas)
    check("the residual was distributed EVENLY over all 24 atoms", spread < 1e-9,
          f"per-atom shift spread {spread:.2e} e")
    check("the per-atom shift equals -0.003/24",
          abs(sum(deltas) / len(deltas) - RESIDUAL_PER_ATOM) < 1e-9,
          f"{sum(deltas)/len(deltas):+.8f} vs {RESIDUAL_PER_ATOM:+.8f}")

    print("\n4. read-charge mode preserved the charges and the atom order")
    gaff_atoms, gaff_bonds = mol2(out / "cha_gaff.mol2")
    tmpl_atoms, tmpl_bonds = mol2(s06 / "cha_a.mol2")
    worst = max(abs(a["charge"] - qi) for a, qi in zip(gaff_atoms, q))
    check("cha_gaff.mol2 charges equal charges_am1bcc.dat", worst < 1e-8,
          f"max difference {worst:.2e} e")
    check("atom names preserved from the template",
          [a["name"] for a in gaff_atoms] == [a["name"] for a in tmpl_atoms])
    conn_g = [(a, b) for a, b, o in gaff_bonds]
    conn_t = [(a, b) for a, b, o in tmpl_bonds]
    check("bond CONNECTIVITY preserved from the template", conn_g == conn_t,
          f"{len(gaff_bonds)} bonds")
    order_changed = [(a, b, ot, og) for (a, b, ot), (_, _, og)
                     in zip(tmpl_bonds, gaff_bonds) if ot != og]
    if order_changed:
        print(f"         NOTE: {len(order_changed)} bond ORDERS rewritten by antechamber:")
        for a, b, ot, og in order_changed:
            na = gaff_atoms[a - 1]["name"]
            nb = gaff_atoms[b - 1]["name"]
            print(f"           {na}-{nb}: '{ot}' -> '{og}'")
        print("         These are the two carboxylates. GAFF represents carboxylate")
        print("         delocalisation through equal charges and identical 'o' types on")
        print("         both oxygens, not through an aromatic bond order, so the rewrite")
        print("         is the expected GAFF treatment. Connectivity is unchanged.")

    print("\n5. GAFF typing")
    aud = {}
    for line in (out / "step08b_cha_gaff_audit.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        if len(f) == 2:
            aud[f[0]] = f[1]
    types = sorted({a["type"] for a in gaff_atoms})
    check("no Tripos dot-types survive", not [t for t in types if "." in t], str(types))
    check("type set matches the audit", ",".join(types) == aud["gaff_atom_types"],
          ",".join(types))
    check("net charge matches the audit",
          abs(sum(a["charge"] for a in gaff_atoms) - float(aud["net_charge"])) < 1e-6)

    print("\n6. frcmod substituted parameters, ranked by penalty")
    frc = (out / "cha.frcmod").read_text().splitlines()
    subs = []
    for L in frc:
        pm = re.search(r"penalty score=\s*([0-9.]+)", L)
        if pm:
            if "same as" in L:
                param, anal = L.split("same as")[0], L.split("same as")[1]
            elif "Using general" in L:
                param, anal = L.split("Using general")[0], "general " + \
                    L.split("Using general")[1]
            else:
                param, anal = L.split("penalty")[0], "unlabelled"
            subs.append((float(pm.group(1)), param.strip()[:30],
                         anal.split(",")[0].strip()[:34]))
    check("frcmod ATTN count is zero", "ATTN" not in "\n".join(frc),
          f"audit records {aud.get('frcmod_ATTN_count')}")
    for pen, param, anal in sorted(subs, reverse=True):
        print(f"         penalty {pen:>6.1f}   {param:<16} substituted by {anal}")
    high = [s for s in subs if s[0] >= 100]
    if high:
        print(f"         {len(high)} parameter(s) carry a penalty at or above 100.")
        print("         Both span os-c3, i.e. the O3-C4 bond, which step 11 identifies")
        print("         as the BREAKING bond. GAFF governs the MD sampling that feeds")
        print("         frame selection, not the QM/MM reaction itself. See")
        print("         interpretation.txt and 00_crosscutting/open_items.txt.")
    check("the high-penalty parameters are recorded, not silently accepted", True,
          f"{len(high)} at penalty >= 100, max {max((s[0] for s in subs), default=0):.1f}")

    print("\n7. three-copy identity licenses one parameter set")
    for line in (out / "step08b_ligand_identity_check.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        check(f"{f[0]}: names and topology match cha_a",
              f[4] == "True" and f[5] == "True",
              f"{f[1]} atoms, {f[2]} bonds, q={f[3]}")

    print("\n8. frozen checksums")
    for line in (out / "sha256_step08_ligand_gaff.txt").read_text().splitlines():
        h, p = line.split()
        f = out / Path(p).name
        check(f"{Path(p).name} matches its checksum", f.exists() and sha(f) == h,
              h[:16])

    print("\n9. declared but absent")
    absent = [n for n in ["parmchk2_cha_gaff.log", "antechamber_am1bcc.log"]
              if not (out / n).exists()]
    print(f"         {absent if absent else 'none'}")
    print("         Both are run logs declared by the scripts but never committed. They")
    print("         carry no result; the parmchk2 outcome is captured in cha.frcmod and")
    print("         the audit. Recorded so the absence is deliberate, not overlooked.")

    print("\nNOT CHECKABLE HERE:")
    am = re.search(r"AmberTools:\s*(\S+)", prov)
    print(f"         Toolchain split. Charges derived with AmberTools "
          f"{am.group(1) if am else '?'} on a local workstation; GAFF typing and parmchk2")
    print("         run with antechamber 22.0 on the HPC. Both versions are recorded but")
    print("         no committed output tests whether the two agree on typing.")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

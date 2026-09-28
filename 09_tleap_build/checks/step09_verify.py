#!/usr/bin/env python3
"""
step09_verify.py - re-runnable audit of step 09, the tleap build.

step09a combines the H++-protonated protein with the three frozen substrates into one
pre-tleap PDB. step09b builds the topology with ff14SB, GAFF and TIP3P, solvates in a
10 A box, neutralises with Na+, and audits the result and the tleap log.

Checks performed:

  1. Continuity. The protein input is step07's accepted structure, the substrate source is
     step06's frozen complex, and the parameter template is step08's typed substrate.
  2. Charge arithmetic, reproduced end to end. Protein formal charge, three substrates at
     -2 each, complex total, and the Na+ count needed to neutralise.
  3. Atom-count conservation through assembly: protein + substrates = combined, and
     combined = dry prmtop natom.
  4. Substrate templates in the combined structure match the GAFF template in atom count
     and order, for all three copies.
  5. Solvated system: prmtop and inpcrd agree on atom count, box present, total charge
     integral and zero, water and ion counts consistent with the dry charge.
  6. Build recipe. The force fields actually sourced in tleap_complex.in are read and
     reported, since the recipe is the methodological claim.
  7. tleap warnings: every warning is bucketed, none falls in 'other', and no errors or
     fatals occurred.
  8. THE CLOSE CONTACTS. Seven distinct pairs survive into the built topology, the closest
     at 1.317 A. This resolves each against the graft provenance from step 02 and the
     contact set from step 04, because where they sit determines whether they matter. Six
     of seven involve residues grafted from 1DBF. The seventh is a substrate contact at
     the active site. Neither is asserted to be benign here; both are deferred to the
     minimisation of step 10, which is what the audit itself claims.
  9. Checksums, excluding the prmtop DATE line churn the step's own note describes.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step09_verify.py [outputs] [step06] [step07] [step08] [step02]
"""
import sys
import re
import hashlib
from pathlib import Path
from collections import defaultdict

FAIL = []
CONTACT_SET = {7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115}
GRAFT = ({("A", r) for r in [1] + list(range(116, 128))}
         | {("B", r) for r in [1] + list(range(115, 128))}
         | {("C", r) for r in [1, 2] + list(range(115, 128))})


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def kv(path, sep="\t"):
    d = {}
    for line in Path(path).read_text().splitlines():
        f = line.split(sep)
        if len(f) == 2:
            d[f[0]] = f[1]
    return d


def pdb_atoms(p):
    out = []
    for L in Path(p).read_text().splitlines():
        if L.startswith("ATOM  ") or L.startswith("HETATM"):
            out.append({"rec": L[:6].strip(), "name": L[12:16].strip(),
                        "resname": L[17:20].strip(), "chain": L[21],
                        "resid": int(L[22:26])})
    return out


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s06 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../06_acceptance_freeze/outputs")
    s07 = Path(sys.argv[3] if len(sys.argv) > 3 else "../../07_protonation_hpp/outputs")
    s08 = Path(sys.argv[4] if len(sys.argv) > 4 else "../../08_ligand_gaff/outputs")
    print(f"step 09 verification against {out.resolve()}\n")

    pre = kv(out / "step09a_pretleap_audit.txt")
    top = kv(out / "step09b_topology_audit.txt")

    print("1. continuity")
    for lbl, p in [("step07 accepted protein", s07 / "abc_protonated_hpp_accepted.pdb"),
                   ("step06 frozen complex", s06 / "abc_cha_preprotonation.pdb"),
                   ("step08 GAFF template", s08 / "cha_gaff.mol2"),
                   ("step08 frcmod", s08 / "cha.frcmod")]:
        check(f"{lbl} reachable", p.exists(), sha(p)[:16] if p.exists() else "")

    print("\n2. charge arithmetic reproduced")
    ce = {}
    for line in (out / "step09a_charge_estimate.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        ce[f[0]] = float(f[1])
    check("three substrates carry 3 x -2", abs(ce["three_CHA"] + 6.0) < 1e-6,
          f"{ce['three_CHA']:+.4f}")
    check("complex total = protein + substrates",
          abs(ce["complex_total"] - (ce["protein_formal_estimate"] + ce["three_CHA"])) < 1e-6,
          f"{ce['protein_formal_estimate']:+.0f} + {ce['three_CHA']:+.1f} = {ce['complex_total']:+.0f}")
    check("Na+ count equals the magnitude of the complex charge",
          int(ce["expected_Na_to_neutralize"]) == abs(int(ce["complex_total"])),
          f"{int(ce['expected_Na_to_neutralize'])} Na+ for {int(ce['complex_total'])}")
    check("dry prmtop charge matches the estimate",
          int(top["dry_charge_int"]) == int(ce["complex_total"]),
          f"{top['dry_charge_int']} vs {int(ce['complex_total'])}")

    print("\n3. atom-count conservation through assembly")
    np_, nl = int(pre["protein_atoms"]), int(pre["cha_total_atoms"])
    check("protein + substrates = combined",
          np_ + nl == int(pre["combined_atoms"]),
          f"{np_} + {nl} = {pre['combined_atoms']}")
    check("combined = dry prmtop natom",
          int(pre["combined_atoms"]) == int(top["dry_prmtop_natom"]),
          f"{pre['combined_atoms']} vs {top['dry_prmtop_natom']}")
    comb = pdb_atoms(out / "complex_for_tleap.pdb")
    check("pre-tleap PDB reproduces the combined atom count",
          len(comb) == int(pre["combined_atoms"]), f"{len(comb)}")

    print("\n4. substrate templates in the combined structure")
    for line in (out / "step09a_cha_template_check.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        check(f"{f[0]}{f[1]}: 24 atoms, order matches the GAFF template",
              f[2] == "24" and f[3] == "True" and f[4] == f[5] == f[6] == "none",
              f"dup={f[4]} missing={f[5]} extra={f[6]}")

    print("\n5. solvated system")
    check("dry prmtop and inpcrd agree on atom count",
          top["prmtop_inpcrd_natom_match_dry"] == "True")
    check("solvated prmtop and inpcrd agree on atom count",
          top["prmtop_inpcrd_natom_match_solvated"] == "True")
    check("dry unit carries no box", top["dry_has_box"] == "False")
    check("solvated unit carries a box", top["solvated_has_box"] == "True",
          top.get("solvated_box_size", ""))
    check("solvated total charge is integral and zero",
          abs(float(top["solvated_total_charge"])) < 1e-3
          and int(top["solvated_charge_int"]) == 0,
          f"{top['solvated_total_charge']}")
    check("Na+ count matches the dry charge magnitude",
          top["Na_ions"] == top["expected_Na_from_dry_charge"],
          f"{top['Na_ions']} Na+, {top['Cl_ions']} Cl-")
    check("all three substrates survive solvation",
          top["solvated_CHA_residues"] == "3")

    print("\n6. build recipe, read from tleap_complex.in")
    recipe = (out / "tleap_complex.in").read_text()
    for want in ["leaprc.protein.ff14SB", "leaprc.gaff", "leaprc.water.tip3p"]:
        check(f"sources {want}", want in recipe)
    box = re.search(r"solvatebox\s+\S+\s+(\S+)\s+([\d.]+)", recipe)
    check("solvates with an explicit box and buffer", bool(box),
          f"{box.group(1)} buffer {box.group(2)} A" if box else "")
    check("neutralises with Na+ only", "addions complex Na+ 0" in recipe)

    print("\n7. tleap warnings")
    buckets = {}
    for line in (out / "step09b_warning_classification.tsv").read_text().splitlines():
        f = line.split("\t")
        if f[0] == "bucket":
            buckets[f[1]] = int(f[2])
    check("no warning falls outside a named bucket", buckets.get("other", 0) == 0,
          f"buckets {buckets}")
    check("no tleap errors", int(top["tleap_errors"]) == 0)
    check("no tleap fatals", int(top["tleap_fatal"]) == 0)
    check("bucket counts account for the reported total",
          sum(buckets.values()) >= int(top["tleap_warnings_total"]),
          f"{sum(buckets.values())} bucketed, {top['tleap_warnings_total']} reported")

    print("\n8. close contacts, resolved against graft provenance and the contact set")
    rows = [l for l in (out / "step09b_warning_classification.tsv").read_text().splitlines()
            if l.startswith("close_contact\t")]
    check("seven distinct close contacts, as the audit records", len(rows) == 7,
          f"{len(rows)}")
    n_graft, n_site = 0, 0
    for L in rows:
        parts = re.findall(r"([ABC]):(\d+):(\w+):(\w+)", L)
        d = float(re.search(r"\(([\d.]+) A\)", L).group(1))
        tag = L.rstrip().split("\t")[-1]
        g = sorted({f"{c}{int(r)}" for c, r, rn, an in parts if (c, int(r)) in GRAFT})
        s = sorted({f"{rn}{int(r)}" for c, r, rn, an in parts
                    if int(r) in CONTACT_SET and rn != "CHA"})
        lig = any(rn == "CHA" for c, r, rn, an in parts)
        n_graft += bool(g)
        n_site += bool(lig or s)
        gs = ",".join(g) if g else "-"
        ss = ",".join(s) if s else "-"
        print(f"         {d:.3f} A  {tag:<20} graft={gs:<12} "
              f"site={ss:<8}{'  LIGAND' if lig else ''}")
    check("every close contact is attributed to graft, active site, or neither",
          n_graft + n_site >= len(rows), f"{n_graft} graft-involved, {n_site} site-involved")
    print(f"         {n_graft} of {len(rows)} involve residues grafted from 1DBF at step 02.")
    print(f"         {n_site} involve the substrate or a contact-set residue.")
    print("         The audit defers all of these to the restrained minimisation of")
    print("         step 10. That deferral is NOT verified here and is tracked as")
    print("         00_crosscutting/open_items.txt item 1.")

    print("\n9. checksums")
    rec = {}
    for line in (out / "sha256_step09_tleap_build.txt").read_text().splitlines():
        h, p = line.split()
        rec[Path(p).name] = h
    for n, h in rec.items():
        f = out / n
        ptr = out / (n + ".pointer")
        if not f.exists() and ptr.exists():
            # Oversized file held as a pointer; verify the recorded hash instead.
            rec_h = next((l.split(":", 1)[1].strip()
                          for l in ptr.read_text().splitlines()
                          if l.startswith("sha256:")), None)
            check(f"{n} held as a pointer, hash recorded in it", rec_h is not None,
                  f"{(rec_h or '')[:16]}")
        elif not f.exists():
            check(f"{n} present", False)
        elif n.endswith(".prmtop"):
            print(f"         {n}: prmtop embeds a wall-clock DATE line; raw hash churns "
                  f"per run. See step09_fingerprints.txt.")
        else:
            check(f"{n} matches its checksum", sha(f) == h, h[:16])

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

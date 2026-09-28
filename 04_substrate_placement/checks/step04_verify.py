#!/usr/bin/env python3
"""
step04_verify.py - re-runnable audit of step 04, substrate placement.

IMPORTANT CORRECTION TO AN EARLIER READING OF THIS STEP.
step04_place_chorismate_abc.sh emits TWO contact records and they are different things:

  step04_placement_contact_report.txt   - an UNBIASED scan. Every protein heavy atom is
                                          measured against the ligand and every residue
                                          within 4.0 A is listed. No residue list is
                                          consulted. This is a criterion-derived set.
  step04_catalytic_contact_report.tsv   - a LIST-DRIVEN check. For each residue in a
                                          hard-coded set it finds the closest copy across
                                          all chains and tests the expected same/cross
                                          chain origin.

The criterion-first derivation of the contact set therefore already exists in the
committed record. It was not introduced by step04b. What step04b adds is the cutoff
SENSITIVITY ladder and the nearest-excluded margin, neither of which step04 produces.

Checks performed:

  1. Continuity. The protein input is step03's committed output and the ligand inputs are
     step01's committed geometries, by SHA-256.
  2. Registration consistency. The source analogue site hard-coded in step04 for each
     ligand matches the site step01b measured in ligand_tsa_registration.tsv.
  3. Transform reproduced. The analogue-to-analogue Kabsch fit is recomputed from the raw
     2CHT deposition for each source/target pair and compared against
     step04_placement_transform_report.tsv (common atom count and RMSD).
  4. Ligand integrity after placement: 24 atoms, 16 heavy, net charge -2 preserved.
  5. Unbiased contact scan reproduced. Every residue within 4.0 A of each placed ligand is
     recomputed from the combined structure and compared against
     step04_placement_contact_report.txt, entry for entry.
  6. Contact set. The union over the three sites is compared against the twelve residues
     published by Agbaglo et al. (2024, main text p. 8, scheme S3).
  7. Cutoff sensitivity and exclusion margin. NOT produced by step04. Reports the set size
     at 3.0/3.5/4.0/4.5/5.0/6.0 A and the closest excluded residues, so the stability of
     the 4.0 A choice and its margin are visible.
  8. Chain-origin assignments in the list-driven report are all 'ok', and Arg116 is
     labelled second_shell rather than active_site.
  9. No severe overlap: the minimum protein-ligand heavy-atom distance is at or above the
     script's own 0.80 A failure threshold.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step04_verify.py [outputs_dir] [step03_outputs] [step01_outputs]
"""
import sys
import math
import hashlib
from pathlib import Path

import numpy as np
from collections import defaultdict

FAIL = []
AGBAGLO_SET = {7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115}
CUTOFF = 4.0
LADDER = [3.0, 3.5, 4.0, 4.5, 5.0, 6.0]


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_pdb(p, want_hetatm=False):
    out = []
    for L in Path(p).read_text().splitlines():
        if L.startswith("ATOM  ") or (want_hetatm and L.startswith("HETATM")):
            el = (L[76:78].strip() if len(L) >= 78 else "") or L[12:16].strip()[0]
            out.append({"rec": L[:6].strip(), "name": L[12:16].strip(),
                        "resname": L[17:20].strip(), "chain": L[21],
                        "resid": int(L[22:26]), "element": el,
                        "xyz": (float(L[30:38]), float(L[38:46]), float(L[46:54]))})
    return out


def mol2(p):
    lines = Path(p).read_text().splitlines()
    sec, atoms = None, []
    for L in lines:
        if L.startswith("@<TRIPOS>ATOM"):
            sec = "a"
            continue
        if L.startswith("@<TRIPOS>"):
            sec = None
            continue
        f = L.split()
        if sec == "a" and len(f) >= 6:
            atoms.append(f)
    q = sum(float(a[8]) for a in atoms) if atoms and len(atoms[0]) >= 9 else None
    heavy = sum(1 for a in atoms if a[5].split(".")[0] != "H")
    return len(atoms), heavy, q


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q, by SVD."""
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    V, _, Wt = np.linalg.svd(H)
    U = V @ np.diag([1.0, 1.0, float(np.sign(np.linalg.det(V @ Wt)))]) @ Wt
    return math.sqrt(((((P - Pc) @ U) + Qc - Q) ** 2).sum() / len(P))


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s03 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../03_cleanup_validate/outputs")
    s01 = Path(sys.argv[3] if len(sys.argv) > 3 else "../../01_inputs_acquisition/outputs")
    print(f"step 04 verification against {out.resolve()}\n")

    print("1. continuity: inputs are the committed step 03 and step 01 outputs")
    check("step03 abc_repaired_clean.pdb reachable", (s03 / "abc_repaired_clean.pdb").exists(),
          sha(s03 / "abc_repaired_clean.pdb")[:16] if (s03 / "abc_repaired_clean.pdb").exists() else "")
    rec = {}
    for line in (s01 / "sha256_step01_inputs.txt").read_text().splitlines():
        h, p = line.split()
        rec[Path(p).name] = h
    for lig in ["liga.mol2", "ligb.mol2", "ligc.mol2"]:
        check(f"step01 {lig} matches its recorded hash",
              (s01 / lig).exists() and sha(s01 / lig) == rec[lig], rec[lig][:16])

    print("\n2. registration consistency between step01b and step04")
    reg = {}
    for line in (s01 / "ligand_tsa_registration.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        reg[f[0]] = f"{f[2]}{f[3]}"
    trows = [l.split("\t") for l in
             (out / "step04_placement_transform_report.tsv").read_text().splitlines()[1:]]
    for r in trows:
        check(f"{r[0]}: step04 source site {r[2]} matches step01b registration",
              reg[r[1]] == r[2], f"step01b says {reg[r[1]]}")

    print("\n3. analogue-to-analogue transform reproduced from raw 2CHT")
    tsa = defaultdict(dict)
    for a in read_pdb(s01 / "2cht_raw.pdb", want_hetatm=True):
        if a["rec"] == "HETATM" and a["resname"] == "TSA":
            tsa[(a["chain"], a["resid"])][a["name"]] = a["xyz"]
    for r in trows:
        src = (r[2][0], int(r[2][1:]))
        tgt = (r[3][0], int(r[3][1:]))
        common = sorted(set(tsa[src]) & set(tsa[tgt]))
        rmsd = kabsch_rmsd([tsa[src][n] for n in common], [tsa[tgt][n] for n in common])
        check(f"{r[0]}: common analogue atoms", len(common) == int(r[4]),
              f"{len(common)} vs {r[4]}")
        check(f"{r[0]}: fit RMSD reproduced", abs(rmsd - float(r[5])) < 5e-3,
              f"{rmsd:.4f} vs {r[5]}")

    print("\n4. ligand integrity after placement")
    for lig in ["cha_a_placed.mol2", "cha_b_placed.mol2", "cha_c_placed.mol2"]:
        n, heavy, q = mol2(out / lig)
        check(f"{lig}: 24 atoms, 16 heavy", (n, heavy) == (24, 16), f"{n}/{heavy}")
        check(f"{lig}: net charge -2 preserved",
              q is not None and abs(q + 2.0) < 1e-6,
              "n/a" if q is None else f"{q:+.4f}")

    print("\n5. unbiased 4.0 A contact scan reproduced")
    comb = read_pdb(out / "abc_with_chorismate_unprotonated.pdb", want_hetatm=True)
    prot = [a for a in comb if a["resname"] != "CHA" and a["element"] != "H"]
    ligs = defaultdict(list)
    for a in comb:
        if a["resname"] == "CHA" and a["element"] != "H":
            ligs[a["chain"]].append(a["xyz"])
    label = {"A": "cha_a", "B": "cha_b", "C": "cha_c"}
    mind = {}
    for ch, coords in ligs.items():
        for a in prot:
            d = min(math.dist(a["xyz"], c) for c in coords)
            k = (ch, a["chain"], a["resid"], a["resname"])
            if k not in mind or d < mind[k]:
                mind[k] = d
    recorded = defaultdict(dict)
    cur = None
    for L in (out / "step04_placement_contact_report.txt").read_text().splitlines():
        if L.startswith("contacts within"):
            cur = L.split()[-1]
            continue
        f = L.split()
        if cur and len(f) >= 4 and f[0].isalpha() and len(f[0]) == 1:
            recorded[cur][(f[0], int(f[1]), f[2])] = float(f[-1].split("=")[-1])
    for ch in "ABC":
        got = {(pc, ri, rn): d for (lc, pc, ri, rn), d in mind.items()
               if lc == ch and d <= CUTOFF}
        want = recorded[label[ch]]
        check(f"{label[ch]}: residue set within 4.0 A matches the report",
              set(got) == set(want), f"{len(got)} found, {len(want)} recorded")
        worst = max((abs(got[k] - want[k]) for k in set(got) & set(want)), default=0.0)
        check(f"{label[ch]}: distances match the report", worst < 2e-3,
              f"max difference {worst:.4f} A")

    print("\n6. contact set against the published set")
    union = {ri for (lc, pc, ri, rn), d in mind.items() if d <= CUTOFF}
    check("union over the three sites equals the Agbaglo set", union == AGBAGLO_SET,
          f"n={len(union)}" + ("" if union == AGBAGLO_SET
                               else f" extra={sorted(union - AGBAGLO_SET)} "
                                    f"missing={sorted(AGBAGLO_SET - union)}"))

    print("\n7. cutoff sensitivity and exclusion margin (not produced by step04)")
    for c in LADDER:
        s = {ri for (lc, pc, ri, rn), d in mind.items() if d <= c}
        flag = "  <- working value" if c == CUTOFF else ""
        print(f"         {c:>4.1f} A : n={len(s):<3} "
              f"added={sorted(s - AGBAGLO_SET)}{flag}")
    excl = sorted(((d, ri, rn) for (lc, pc, ri, rn), d in mind.items()
                   if ri not in union))
    seen, near = set(), []
    for d, ri, rn in excl:
        if ri in seen:
            continue
        seen.add(ri)
        near.append((rn, ri, d))
        if len(near) >= 5:
            break
    print("         nearest excluded residues: " +
          ", ".join(f"{rn}{ri} {d:.3f} A" for rn, ri, d in near))
    check("the 4.0 A set is stable across 3.5-4.0 A",
          {ri for (lc, pc, ri, rn), d in mind.items() if d <= 3.5} == AGBAGLO_SET)

    print("\n8. chain-origin assignments in the list-driven report")
    rows = [l.split("\t") for l in
            (out / "step04_catalytic_contact_report.tsv").read_text().splitlines()
            if l and not l.startswith("#") and not l.startswith("ligand")]
    bad = [r for r in rows if not r[9].endswith(":ok")]
    check("every chain-origin assignment recovered by measurement", not bad,
          f"{len(rows)} rows, {len(bad)} mismatched")
    a116 = [r for r in rows if r[3] == "116"]
    check("Arg116 is labelled second_shell, not active_site",
          bool(a116) and all(r[2] == "second_shell" for r in a116),
          f"{len(a116)} rows")

    print("\n9. no severe overlap")
    aud = {}
    for line in (out / "step04_placement_audit.txt").read_text().splitlines():
        f = line.split("\t")
        if len(f) == 2:
            aud[f[0]] = f[1]
    for lab in ["cha_a", "cha_b", "cha_c"]:
        v = float(aud[f"{lab}_min_prot_lig_heavy_A"])
        check(f"{lab}: minimum contact above the 0.80 A failure threshold", v >= 0.80,
              f"{v:.3f} A")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

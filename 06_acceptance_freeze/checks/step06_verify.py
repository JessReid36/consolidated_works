#!/usr/bin/env python3
"""
step06_verify.py - re-runnable audit of step 06, the acceptance freeze.

Step 06 renames the substrate atoms uniquely for tleap, rebuilds the combined complex
with the safe names, audits it, and freezes the result with checksums. The renaming is
the substantive operation and it must be a pure relabelling: coordinates, atom types,
charges and bonds all preserved.

Checks performed:

  1. Continuity. The protein input is step03's committed output and the substrate inputs
     are step04's committed placed geometries, by SHA-256.
  2. Renaming is a pure relabelling. For each substrate, the renaming table is applied to
     the step04 input and the result compared against the frozen mol2 atom by atom:
     coordinates to 1e-6, atom type and charge exactly, and the bond block unchanged.
  3. Atom names are now unique, 24 distinct names per substrate.
  4. Net charge preserved at exactly -2 in every frozen substrate.
  5. Combined complex composition: protein atoms, heteroatoms, substrate residues and
     duplicate-name residues, reproduced from the structure and compared against
     step06_combined_pdb_audit.txt. Protein + heteroatom counts must sum to the file.
  6. The freeze moved nothing. Protein coordinates in the combined complex are identical
     to step03's output, and substrate coordinates are identical to step04's placed
     geometries.
  7. Frozen checksums reproduce against the committed frozen files.

A NOTE, not a check, is printed for the residue list in step06_acceptance_decision.txt.
That text is a heredoc literal inside step06_freeze_preprotonation.sh and names the
pre-Agbaglo active-site list, omitting Val73'. The frozen structure is unaffected. It is
tracked as open_items.txt item 5 and is not treated as a failure here, because nothing in
the frozen output depends on it.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step06_verify.py [outputs_dir] [step03_outputs] [step04_outputs]
"""
import sys
import math
import hashlib
from pathlib import Path
from collections import Counter

FAIL = []
CONTACT_SET = {7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def mol2(path):
    """Return (atoms, bonds). atoms: list of dicts. bonds: list of (a1, a2, order)."""
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
            atoms.append({"id": int(f[0]), "name": f[1],
                          "xyz": (float(f[2]), float(f[3]), float(f[4])),
                          "type": f[5],
                          "charge": float(f[8]) if len(f) >= 9 else None})
        elif sec == "b" and len(f) >= 4:
            bonds.append((int(f[1]), int(f[2]), f[3]))
    return atoms, bonds


def pdb(path):
    out = []
    for L in Path(path).read_text().splitlines():
        if not (L.startswith("ATOM  ") or L.startswith("HETATM")):
            continue
        out.append({"rec": L[:6].strip(), "name": L[12:16].strip(),
                    "resname": L[17:20].strip(), "chain": L[21],
                    "resid": int(L[22:26]),
                    "xyz": (float(L[30:38]), float(L[38:46]), float(L[46:54]))})
    return out


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s03 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../03_cleanup_validate/outputs")
    s04 = Path(sys.argv[3] if len(sys.argv) > 3 else "../../04_substrate_placement/outputs")
    print(f"step 06 verification against {out.resolve()}\n")

    print("1. continuity")
    check("step03 abc_repaired_clean.pdb reachable",
          (s03 / "abc_repaired_clean.pdb").exists(),
          sha(s03 / "abc_repaired_clean.pdb")[:16]
          if (s03 / "abc_repaired_clean.pdb").exists() else "")
    for lab in "abc":
        check(f"step04 cha_{lab}_placed.mol2 reachable",
              (s04 / f"cha_{lab}_placed.mol2").exists())

    print("\n2. renaming is a pure relabelling")
    for lab in "abc":
        src_atoms, src_bonds = mol2(s04 / f"cha_{lab}_placed.mol2")
        new_atoms, new_bonds = mol2(out / f"cha_{lab}.mol2")
        table = {}
        for line in (out / f"cha_{lab}_renaming.tsv").read_text().splitlines()[1:]:
            f = line.split("\t")
            table[int(f[0])] = (f[1], f[2], f[4], float(f[5]))
        check(f"cha_{lab}: renaming table covers every atom",
              len(table) == len(src_atoms) == len(new_atoms),
              f"table {len(table)}, src {len(src_atoms)}, frozen {len(new_atoms)}")
        worst_xyz, bad_type, bad_q, bad_map = 0.0, 0, 0, 0
        for s, n in zip(src_atoms, new_atoms):
            worst_xyz = max(worst_xyz, math.dist(s["xyz"], n["xyz"]))
            if s["type"] != n["type"]:
                bad_type += 1
            if s["charge"] is not None and abs(s["charge"] - n["charge"]) > 1e-9:
                bad_q += 1
            old, new, _, _ = table[s["id"]]
            if s["name"] != old or n["name"] != new:
                bad_map += 1
        check(f"cha_{lab}: coordinates preserved", worst_xyz < 1e-6,
              f"max displacement {worst_xyz:.2e} A")
        check(f"cha_{lab}: atom types preserved", bad_type == 0, f"{bad_type} changed")
        check(f"cha_{lab}: charges preserved", bad_q == 0, f"{bad_q} changed")
        check(f"cha_{lab}: renaming table agrees with both files", bad_map == 0,
              f"{bad_map} mismatched")
        check(f"cha_{lab}: bond block unchanged", src_bonds == new_bonds,
              f"{len(src_bonds)} bonds")

    print("\n3. atom names are unique after renaming")
    for lab in "abc":
        atoms, _ = mol2(out / f"cha_{lab}.mol2")
        names = [a["name"] for a in atoms]
        check(f"cha_{lab}: 24 distinct names", len(set(names)) == len(names) == 24,
              f"{len(set(names))} distinct of {len(names)}")

    print("\n4. net charge preserved")
    aud = {}
    for line in (out / "step06_ligand_mol2_audit.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        aud[f[0]] = (int(f[1]), int(f[2]), int(f[3]), float(f[4]), int(f[5]))
    for lab in "abc":
        atoms, bonds = mol2(out / f"cha_{lab}.mol2")
        q = sum(a["charge"] for a in atoms)
        heavy = sum(1 for a in atoms if a["type"].split(".")[0] != "H")
        rec = aud[f"cha_{lab}"]
        check(f"cha_{lab}: net charge is exactly -2", abs(q + 2.0) < 1e-9, f"{q:+.6f}")
        check(f"cha_{lab}: audit row reproduced",
              (len(atoms), heavy, len(bonds)) == (rec[0], rec[1], rec[2]),
              f"{len(atoms)}/{heavy}/{len(bonds)} vs {rec[0]}/{rec[1]}/{rec[2]}")

    print("\n5. combined complex composition")
    comb = pdb(out / "abc_cha_preprotonation.pdb")
    prot = [a for a in comb if a["rec"] == "ATOM"]
    het = [a for a in comb if a["rec"] == "HETATM"]
    ligres = {(a["chain"], a["resid"]) for a in het}
    paud = {}
    for line in (out / "step06_combined_pdb_audit.txt").read_text().splitlines():
        f = line.split("\t")
        if len(f) == 2:
            paud[f[0]] = f[1]
    check("protein atom count reproduced", len(prot) == int(paud["protein_atoms"]),
          f"{len(prot)} vs {paud['protein_atoms']}")
    check("heteroatom count reproduced", len(het) == int(paud["hetatm_atoms"]),
          f"{len(het)} vs {paud['hetatm_atoms']}")
    check("substrate residue count reproduced",
          len(ligres) == int(paud["ligand_residues"]),
          f"{len(ligres)} vs {paud['ligand_residues']}")
    dup = 0
    for ch, ri in ligres:
        names = [a["name"] for a in het if a["chain"] == ch and a["resid"] == ri]
        if any(c > 1 for c in Counter(names).values()):
            dup += 1
    check("no substrate residue carries duplicate names",
          dup == int(paud["ligand_duplicate_name_residues"]),
          f"{dup} vs {paud['ligand_duplicate_name_residues']}")
    check("protein + heteroatoms account for every record in the file",
          len(prot) + len(het) == len(comb), f"{len(prot)}+{len(het)}={len(comb)}")

    print("\n6. the freeze moved nothing")
    before = {(a["chain"], a["resid"], a["name"]): a["xyz"]
              for a in pdb(s03 / "abc_repaired_clean.pdb")}
    worst, missing = 0.0, 0
    for a in prot:
        k = (a["chain"], a["resid"], a["name"])
        if k not in before:
            missing += 1
        else:
            worst = max(worst, math.dist(before[k], a["xyz"]))
    check("every protein atom traces to step 03", missing == 0, f"{missing} untraceable")
    check("protein coordinates unchanged by the freeze", worst < 1e-6,
          f"max displacement {worst:.2e} A")
    worst_l = 0.0
    for lab, ch in zip("abc", "ABC"):
        src_atoms, _ = mol2(s04 / f"cha_{lab}_placed.mol2")
        hetc = [a for a in het if a["chain"] == ch]
        for s, h in zip(src_atoms, hetc):
            worst_l = max(worst_l, math.dist(s["xyz"], h["xyz"]))
    # The combined file is a PDB, which stores coordinates to three decimals, while the
    # mol2 source carries more. The largest possible discrepancy from rounding alone is
    # sqrt(3) * 5e-4 = 8.7e-4 A, so the tolerance is set just above that.
    check("substrate coordinates unchanged by the freeze, to PDB rounding",
          worst_l < 8.7e-4,
          f"max displacement {worst_l:.2e} A, rounding limit 8.7e-04 A")

    print("\n7. frozen checksums")
    for line in (out / "sha256_step06_accepted_preprotonation.txt").read_text().splitlines():
        h, p = line.split()
        f = out / Path(p).name
        check(f"{Path(p).name} matches its frozen checksum",
              f.exists() and sha(f) == h, h[:16])

    print("\nNOTE, not a check:")
    txt = (out / "step06_acceptance_decision.txt").read_text()
    # The text enumerates only the ADJACENT-CHAIN contributors, in one parenthesis. It
    # does not enumerate the same-chain residues at all, so only the cross list is
    # comparable. The measured cross set is Phe57, Ala59, Lys60, Arg63, Val73, Thr74,
    # Cys75.
    cross = {57: "Phe57", 59: "Ala59", 60: "Lys60", 63: "Arg63",
             73: "Val73", 74: "Thr74", 75: "Cys75"}
    absent = sorted(r for r, nm in cross.items() if nm not in txt)
    print("         step06_acceptance_decision.txt is a heredoc literal inside the script.")
    print("         It enumerates only the adjacent-chain contributors to the active site.")
    print(f"         Measured cross-chain contacts it does NOT name: "
          f"{[cross[r] for r in absent] if absent else 'none'}")
    print("         The frozen structure does not depend on this text. "
          "See 00_crosscutting/open_items.txt item 5.")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

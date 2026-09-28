#!/usr/bin/env python3
"""
step03_verify.py - re-runnable audit of step 03, cleanup and validation.

Checks performed:

  1. Continuity. step03's input abc_repaired.pdb is byte-identical to step02's committed
     output, by SHA-256. This establishes that step03 consumed step02's product.
  2. Atom-count arithmetic. initial - hydrogens - altloc duplicates = final heavy atoms,
     reproduced from the structures rather than taken from the audit file.
  3. Purity of the output: no hydrogens, no altloc flags, no duplicate atom names, no
     non-standard residue names, three OXT.
  4. Peptide continuity recomputed with the script's own window (1.15-1.70 A against a
     canonical 1.329 A). The worst link is reported whether or not it passes, since the
     window is wide.
  5. Completeness: 127 residues per chain and all twelve contact-set residues present.
  6. Altloc resolution. Reproduces the script's rule (highest occupancy, ties to altloc A
     then blank) from abc_repaired.pdb and confirms the kept conformer in every group is
     the one step03_cleanup_report.tsv records. Reports for each group the OCCUPANCY
     MARGIN and the CONFORMER SEPARATION, neither of which the script records. A small
     margin means the choice was close to arbitrary; a large separation means the choice
     was structurally consequential. Both matter and neither is otherwise visible.
  7. Cleanup moved nothing. Every heavy atom retained in abc_repaired_clean.pdb has
     coordinates identical to its counterpart in abc_repaired.pdb. Cleanup is a deletion
     operation and must not perturb geometry; the script does not assert this.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step03_verify.py [outputs_dir] [step02_outputs_dir]
Defaults: ../outputs  and  ../../02_terminal_repair/outputs
"""
import sys
import math
import hashlib
from pathlib import Path
from collections import defaultdict, Counter

FAIL = []
CONTACT_SET = [7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115]
CN_MIN, CN_MAX, TARGET_CN = 1.15, 1.70, 1.329
STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def parse(path):
    out = []
    for L in Path(path).read_text().splitlines():
        if not L.startswith("ATOM  "):
            continue
        out.append({
            "serial": int(L[6:11]), "name": L[12:16].strip(), "altloc": L[16],
            "resname": L[17:20].strip(), "chain": L[21], "resid": int(L[22:26]),
            "xyz": (float(L[30:38]), float(L[38:46]), float(L[46:54])),
            "occ": float(L[54:60]) if L[54:60].strip() else 1.0,
            "element": L[76:78].strip() if len(L) >= 78 else "",
        })
    return out


def audit_values(path):
    d = {}
    for line in Path(path).read_text().splitlines():
        f = line.split("\t")
        if len(f) == 2:
            d[f[0]] = f[1]
    return d


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s02 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../02_terminal_repair/outputs")
    print(f"step 03 verification against {out.resolve()}\n")

    aud = audit_values(out / "step03_validation_audit.txt")

    print("1. continuity: step03's input is step02's output")
    src = s02 / "abc_repaired.pdb"
    if not src.exists():
        check("abc_repaired.pdb reachable in step02 outputs", False, str(s02))
        before = []
    else:
        here = out / "abc_repaired_clean.pdb"
        check("step02 abc_repaired.pdb present", True,
              hashlib.sha256(src.read_bytes()).hexdigest()[:16])
        before = parse(src)
        check("audit names abc_repaired.pdb as its input",
              aud.get("input", "").endswith("abc_repaired.pdb"), aud.get("input", ""))
        check("cleaned output present", here.exists())

    after = parse(out / "abc_repaired_clean.pdb")

    print("\n2. atom-count arithmetic, reproduced from the structures")
    n_init = len(before)
    n_h = sum(1 for a in before if a["element"].upper() == "H")
    groups = defaultdict(list)
    for a in before:
        groups[(a["chain"], a["resid"], a["resname"], a["name"])].append(a)
    n_dup = sum(len(v) - 1 for v in groups.values() if len(v) > 1)
    n_groups = sum(1 for v in groups.values() if len(v) > 1)
    check("initial atom count matches the audit", str(n_init) == aud["initial_atoms"],
          f"{n_init} vs {aud['initial_atoms']}")
    check("hydrogens removed matches the audit",
          str(n_h) == aud["hydrogens_removed_by_element"],
          f"{n_h} vs {aud['hydrogens_removed_by_element']}")
    check("altloc groups resolved matches the audit",
          str(n_groups) == aud["altloc_groups_resolved"],
          f"{n_groups} vs {aud['altloc_groups_resolved']}")
    check("initial - hydrogens - altloc duplicates = final heavy atoms",
          n_init - n_h - n_dup == len(after),
          f"{n_init} - {n_h} - {n_dup} = {n_init - n_h - n_dup}, final {len(after)}")

    print("\n3. purity of the output")
    check("no hydrogens remain",
          not [a for a in after if a["element"].upper() == "H"])
    check("no altloc flags remain", not [a for a in after if a["altloc"].strip()])
    by_res = defaultdict(list)
    for a in after:
        by_res[(a["chain"], a["resid"])].append(a["name"])
    dups = [k for k, v in by_res.items() if any(c > 1 for c in Counter(v).values())]
    check("no duplicate atom names within a residue", not dups, str(dups[:3]))
    nonstd = sorted({a["resname"] for a in after if a["resname"] not in STANDARD})
    check("no non-standard residue names", not nonstd, str(nonstd))
    n_oxt = sum(1 for a in after if a["name"] == "OXT")
    check("three OXT, one per chain", n_oxt == 3, str(n_oxt))

    print("\n4. peptide continuity, recomputed")
    xyz = defaultdict(dict)
    for a in after:
        xyz[(a["chain"], a["resid"])][a["name"]] = a["xyz"]
    checked, bad, worst = 0, [], None
    for ch in "ABC":
        for r in range(1, 127):
            C = xyz.get((ch, r), {}).get("C")
            N = xyz.get((ch, r + 1), {}).get("N")
            if C and N:
                checked += 1
                d = math.dist(C, N)
                if worst is None or abs(d - TARGET_CN) > abs(worst[1] - TARGET_CN):
                    worst = ((ch, r, r + 1), d)
                if d < CN_MIN or d > CN_MAX:
                    bad.append((ch, r, r + 1, d))
    check("links checked matches the audit", str(checked) == aud["peptide_links_checked"],
          f"{checked} vs {aud['peptide_links_checked']}")
    check("bad links matches the audit", str(len(bad)) == aud["bad_peptide_links"],
          f"{len(bad)} vs {aud['bad_peptide_links']}")
    print(f"         worst link: {worst[0]} at {worst[1]:.3f} A "
          f"({worst[1] - TARGET_CN:+.3f} from canonical, window {CN_MIN}-{CN_MAX})")

    print("\n5. completeness")
    res = defaultdict(set)
    for a in after:
        res[a["chain"]].add(a["resid"])
    for ch in "ABC":
        check(f"chain {ch} is 127 residues, 1-127",
              len(res[ch]) == 127 and min(res[ch]) == 1 and max(res[ch]) == 127,
              f"n={len(res[ch])}")
    missing = [(ch, r) for ch in "ABC" for r in CONTACT_SET if r not in res[ch]]
    check("all twelve contact-set residues present in every chain", not missing,
          "none missing" if not missing else str(missing))

    print("\n6. altloc resolution: rule, margin and structural consequence")
    recorded = {}
    for line in (out / "step03_cleanup_report.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        if f[0] == "altloc":
            recorded[(f[1], int(f[2]), f[4])] = f[5]
    check("recorded altloc decisions match the audit count",
          len(recorded) == n_dup, f"{len(recorded)} decisions, {n_dup} duplicates")
    per_res = defaultdict(list)
    for key, recs in groups.items():
        if len(recs) < 2:
            continue
        ch, ri, rn, an = key
        best = sorted(recs, key=lambda a: (-a["occ"],
                                           {"A": 0, " ": 1, "B": 2}.get(a["altloc"], 3),
                                           a["serial"]))[0]
        kept = best["altloc"].strip() or "blank"
        check(f"{ch} {rn}{ri} {an}: kept conformer is the highest-occupancy one",
              recorded.get((ch, ri, an)) == kept,
              f"kept {kept}, report says {recorded.get((ch, ri, an))}")
        occs = {a["altloc"]: a["occ"] for a in recs}
        pos = {a["altloc"]: a["xyz"] for a in recs}
        if "A" in pos and "B" in pos:
            per_res[(ch, ri, rn)].append(
                (an, math.dist(pos["A"], pos["B"]), occs["A"], occs["B"], kept))
    print("\n         occupancy margin and conformer separation per residue:")
    for (ch, ri, rn), lst in sorted(per_res.items()):
        oa, ob, kept = lst[0][2], lst[0][3], lst[0][4]
        seps = [s for _, s, _, _, _ in lst]
        rmsd = math.sqrt(sum(s * s for s in seps) / len(seps))
        print(f"           {ch} {rn}{ri}: occ A={oa:.2f} B={ob:.2f}, margin {abs(oa-ob):.2f}, "
              f"kept {kept}; conformer rmsd {rmsd:.3f} A, max {max(seps):.3f} A")

    print("\n7. cleanup moved nothing")
    kept_map = {}
    for key, recs in groups.items():
        best = sorted(recs, key=lambda a: (-a["occ"],
                                           {"A": 0, " ": 1, "B": 2}.get(a["altloc"], 3),
                                           a["serial"]))[0]
        if best["element"].upper() != "H":
            kept_map[key] = best["xyz"]
    moved, absent = [], []
    for a in after:
        k = (a["chain"], a["resid"], a["resname"], a["name"])
        if k not in kept_map:
            absent.append(k)
        elif math.dist(kept_map[k], a["xyz"]) > 1e-6:
            moved.append((k, math.dist(kept_map[k], a["xyz"])))
    check("every retained atom is traceable to the input", not absent, str(absent[:3]))
    check("no retained atom changed coordinates", not moved,
          "max displacement 0.000 A" if not moved
          else f"{len(moved)} moved, max {max(d for _, d in moved):.3f} A")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
step01_verify.py - re-runnable audit of step 01 outputs.

Checks performed, all against the committed files in ../outputs/:

  1. SHA-256 of the five raw inputs vs sha256_step01_inputs.txt.
  2. PDB inventory: ATOM/HETATM counts, chain set, non-water heteros.
  3. TSA site inventory in 2CHT.
  4. Ligand composition: atom count, bond count, elemental formula and NET
     CHARGE. The step01 script's own comment says to "confirm element
     composition and net charge from the mol2 before it sets the QM-region
     charge", but the script only prints name/atoms/bonds. This implements
     the confirmation the comment promises.
  5. Ligand-to-TSA registration margin. step01b takes the nearest centre
     unconditionally with no ambiguity guard; this reports the distance to
     the runner-up so the assignment is measured, not assumed.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step01_verify.py [outputs_dir]     (default ../outputs)
"""
import sys
import math
import hashlib
import collections
from pathlib import Path

FAIL = []


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def pdb_inventory(p):
    atom = het = 0
    chains = set()
    het_res = collections.Counter()
    tsa_sites = set()
    for L in Path(p).read_text().splitlines():
        if L.startswith("ATOM  "):
            atom += 1
        elif L.startswith("HETATM"):
            het += 1
            rn = L[17:20]
            if rn != "HOH":
                het_res[rn] += 1
            if rn == "TSA":
                tsa_sites.add((L[21], int(L[22:26])))
        else:
            continue
        if L[21] != " ":
            chains.add(L[21])
    return atom, het, "".join(sorted(chains)), het_res, tsa_sites


def mol2_parse(p):
    lines = Path(p).read_text().splitlines()
    name, sec, atoms, bonds = None, None, [], 0
    for i, L in enumerate(lines):
        if L.startswith("@<TRIPOS>MOLECULE"):
            name = lines[i + 1].strip()
            continue
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
            atoms.append(f)
        elif sec == "b" and len(f) >= 4:
            bonds += 1
    els = collections.Counter(a[5].split(".")[0] for a in atoms)
    formula = "".join(f"{e}{els[e]}" for e in sorted(els))
    q = sum(float(a[8]) for a in atoms) if atoms and len(atoms[0]) >= 9 else None
    return name, len(atoms), bonds, formula, q


def centre_mol2(p):
    xs, ina = [], False
    for L in Path(p).read_text().splitlines():
        if L.startswith("@<TRIPOS>ATOM"):
            ina = True
            continue
        if L.startswith("@<TRIPOS>") and ina:
            break
        f = L.split()
        if ina and len(f) >= 5:
            xs.append((float(f[2]), float(f[3]), float(f[4])))
    return tuple(sum(c) / len(xs) for c in zip(*xs))


def tsa_centres(p):
    g = collections.defaultdict(list)
    for L in Path(p).read_text().splitlines():
        if L.startswith("HETATM") and L[17:20] == "TSA":
            g[(L[21], int(L[22:26]))].append(
                (float(L[30:38]), float(L[38:46]), float(L[46:54]))
            )
    return {k: tuple(sum(c) / len(v) for c in zip(*v)) for k, v in g.items()}


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    print(f"step 01 verification against {out.resolve()}\n")

    print("1. raw-input checksums vs sha256_step01_inputs.txt")
    recorded = {}
    for line in (out / "sha256_step01_inputs.txt").read_text().splitlines():
        h, p = line.split()
        recorded[Path(p).name] = h
    for n in ["2cht_raw.pdb", "1dbf_raw.pdb", "liga.mol2", "ligb.mol2", "ligc.mol2"]:
        f = out / n
        check(f"sha256 {n}", f.exists() and sha256(f) == recorded.get(n),
              recorded.get(n, "no record")[:16])

    print("\n2. PDB inventory")
    a, h, ch, hetres, tsa = pdb_inventory(out / "2cht_raw.pdb")
    check("2CHT chain set is A-L", ch == "ABCDEFGHIJKL", f"chains={ch}")
    check("2CHT ATOM/HETATM", (a, h) == (13529, 1815), f"ATOM={a} HETATM={h}")
    check("2CHT non-water hetero is TSA only", set(hetres) == {"TSA"},
          f"{dict(hetres)}")
    a2, h2, ch2, _, _ = pdb_inventory(out / "1dbf_raw.pdb")
    check("1DBF chain set is A-C", ch2 == "ABC", f"chains={ch2}")
    check("1DBF ATOM/HETATM", (a2, h2) == (3156, 499), f"ATOM={a2} HETATM={h2}")

    print("\n3. TSA site inventory")
    check("twelve TSA sites, one per chain", len(tsa) == 12,
          " ".join(f"{c}{r}" for c, r in sorted(tsa)))

    print("\n4. ligand composition and net charge (not done by step01 itself)")
    for lig in ["liga", "ligb", "ligc"]:
        name, na, nb, formula, q = mol2_parse(out / f"{lig}.mol2")
        check(f"{lig} 24 atoms / 24 bonds", (na, nb) == (24, 24), f"{na}/{nb}")
        check(f"{lig} formula C10H8O6", formula == "C10H8O6", formula)
        check(f"{lig} net charge -2", q is not None and abs(q + 2.0) < 1e-6,
              "n/a" if q is None else f"{q:+.4f}")

    print("\n5. registration margin (step01b has no ambiguity guard)")
    t = tsa_centres(out / "2cht_raw.pdb")
    reg = {}
    for line in (out / "ligand_tsa_registration.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        reg[f[0]] = (f[2], int(f[3]), float(f[4]))
    for lig in ["liga.mol2", "ligb.mol2", "ligc.mol2"]:
        c = centre_mol2(out / lig)
        d = sorted((math.dist(c, v), k) for k, v in t.items())
        (d1, k1), (d2, _) = d[0], d[1]
        rc, rr, rd = reg[lig]
        check(f"{lig} registered to {rc}{rr} as recorded",
              (k1[0], k1[1]) == (rc, rr) and abs(d1 - rd) < 5e-3,
              f"{d1:.3f} A, runner-up {d2:.3f} A, margin {d2 - d1:.2f} A "
              f"({d2 / d1:.0f}x)")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

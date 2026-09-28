#!/usr/bin/env bash
# Step 01 addendum — record which crystallographic TSA site each CP2K chorismate
# ligand is registered to. The CP2K coordinates are pre-aligned to specific 2CHT
# active sites; this mapping is the provenance for the source sites used later
# when superposing chorismate onto the A/B/C TSA pose. Nearest-centre match only.
set -euo pipefail

root="$HOME/system_development"
pdb="$root/01_inputs/structures/2cht_raw.pdb"
ligdir="$root/01_inputs/ligands"
out="$root/00_admin/ligand_tsa_registration.tsv"

python3 - "$pdb" "$ligdir" "$out" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict

pdb, ligdir, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])

def dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))

def tsa_centres(path):
    groups = defaultdict(list)
    for line in path.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20] == "TSA":
            key = (line[21], int(line[22:26]))
            groups[key].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return {k: tuple(sum(c)/len(v) for c in zip(*v)) for k, v in groups.items()}

def mol2_centre(path):
    xs = []
    inatoms = False
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>ATOM"): inatoms = True; continue
        if line.startswith("@<TRIPOS>") and inatoms: break
        if inatoms and len(line.split()) >= 5:
            p = line.split(); xs.append((float(p[2]), float(p[3]), float(p[4])))
    return tuple(sum(c)/len(xs) for c in zip(*xs)), len(xs)

tsa = tsa_centres(pdb)
rows = []
for lig in ["liga.mol2", "ligb.mol2", "ligc.mol2"]:
    centre, n = mol2_centre(ligdir / lig)
    dists = sorted(
        (dist(centre, c), chain, resid)
        for (chain, resid), c in tsa.items()
    )
    d, chain, resid = dists[0]
    rows.append((lig, n, chain, resid, d))
    print(f"{lig:10s} atoms={n:2d}  registered to TSA {chain}{resid}  (centre offset {d:.3f} A)")

with out.open("w") as f:
    f.write("ligand\tatoms\tsource_tsa_chain\tsource_tsa_resid\tcentre_offset_A\n")
    for lig, n, chain, resid, d in rows:
        f.write(f"{lig}\t{n}\t{chain}\t{resid}\t{d:.3f}\n")
print("\nWROTE " + str(out))
PY
echo "STEP 01b DONE"

#!/usr/bin/env bash
# Step 05 - provenance check for the site-C Arg63' rotamer asymmetry.
#
# Question: is the swung-away Arg63' in our placed cha_c site INHERITED from raw
# 2CHT (crystallographic) or INTRODUCED by repair/placement?
#
# Approach: in RAW 2CHT, measure each of the 12 TSA sites against the nearest
# Arg63 across all chains (guanidinium N vs any-atom), recording distance and
# which atom is closest. If the TSA site our cha_c derives from (source L211,
# per step 01b) already shows the guanidinium swung away (nearest approach via a
# non-N atom at long distance) while other sites show short N contacts, the
# asymmetry is crystallographic. We also confirm the placed complex reproduces it.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
raw2cht="$root/01_inputs/structures/2cht_raw.pdb"
placed="$root/02_preparation/ligand_placement/abc_with_chorismate_unprotonated.pdb"
admin="$root/00_admin"

python3 - "$raw2cht" "$placed" "$admin" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict
import numpy as np

raw2cht, placed, admin = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out = admin / "step05_arg63_provenance_report.tsv"

def dist(a, b): return math.sqrt(float(((a - b) ** 2).sum()))
guanidinium_N = {"NE", "NH1", "NH2"}

# --- RAW 2CHT: TSA sites and all Arg63 copies ---
tsa = defaultdict(list)       # (chain,resid) -> list of xyz (TSA heavy atoms)
arg63 = defaultdict(dict)     # chain -> {atomname: xyz}
for line in raw2cht.read_text().splitlines():
    if line.startswith("HETATM") and line[17:20].strip() == "TSA":
        tsa[(line[21], int(line[22:26]))].append(
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    elif line.startswith("ATOM  ") and line[17:20].strip() == "ARG" and int(line[22:26]) == 63:
        arg63[line[21]][line[12:16].strip()] = \
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])

def nearest_arg63(tsa_atoms):
    # across all Arg63 copies: closest guanidinium-N approach, and closest any-atom
    best_N = None; best_any = None
    for ch, atoms in arg63.items():
        for aname, axyz in atoms.items():
            for t in tsa_atoms:
                dd = dist(t, axyz)
                if best_any is None or dd < best_any[0]:
                    best_any = (dd, ch, aname)
                if aname in guanidinium_N and (best_N is None or dd < best_N[0]):
                    best_N = (dd, ch, aname)
    return best_N, best_any

print("RAW 2CHT: each TSA site vs nearest Arg63 (across all chains)")
print(f"  {'TSA':6s} {'N_dist':>7s} {'N_chain':>7s} {'any_dist':>8s} {'any_chain':>9s} {'any_atom':>8s}")
rows = []
for key in sorted(tsa):
    bN, bA = nearest_arg63(tsa[key])
    site = f"{key[0]}{key[1]}"
    nd = f"{bN[0]:.3f}" if bN else "NA"; nc = bN[1] if bN else "NA"
    ad = f"{bA[0]:.3f}"; ac = bA[1]; aa = bA[2]
    swung = "SWUNG" if (bN and bN[0] > 4.5) else ""
    print(f"  {site:6s} {nd:>7s} {nc:>7s} {ad:>8s} {ac:>9s} {aa:>8s}  {swung}")
    rows.append(("raw_2cht", site, nd, nc, ad, ac, aa, swung))

# --- PLACED complex: each CHA vs nearest Arg63 (across all chains) ---
cha = defaultdict(list)
parg63 = defaultdict(dict)
for line in placed.read_text().splitlines():
    if line.startswith("HETATM") and line[17:20].strip() == "CHA":
        cha[(line[21], int(line[22:26]))].append(
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]))
    elif line.startswith("ATOM  ") and line[17:20].strip() == "ARG" and int(line[22:26]) == 63:
        parg63[line[21]][line[12:16].strip()] = \
            np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])

def nearest_parg63(cha_atoms):
    best_N = None; best_any = None
    for ch, atoms in parg63.items():
        for aname, axyz in atoms.items():
            for t in cha_atoms:
                dd = dist(t, axyz)
                if best_any is None or dd < best_any[0]:
                    best_any = (dd, ch, aname)
                if aname in guanidinium_N and (best_N is None or dd < best_N[0]):
                    best_N = (dd, ch, aname)
    return best_N, best_any

print("\nPLACED complex: each CHA site vs nearest Arg63 (across all chains)")
print(f"  {'CHA':6s} {'N_dist':>7s} {'N_chain':>7s} {'any_dist':>8s} {'any_chain':>9s} {'any_atom':>8s}")
for key in sorted(cha):
    bN, bA = nearest_parg63(cha[key])
    site = f"{key[0]}{key[1]}"
    nd = f"{bN[0]:.3f}" if bN else "NA"; nc = bN[1] if bN else "NA"
    ad = f"{bA[0]:.3f}"; ac = bA[1]; aa = bA[2]
    swung = "SWUNG" if (bN and bN[0] > 4.5) else ""
    print(f"  {site:6s} {nd:>7s} {nc:>7s} {ad:>8s} {ac:>9s} {aa:>8s}  {swung}")
    rows.append(("placed", site, nd, nc, ad, ac, aa, swung))

with out.open("w") as f:
    f.write("source\tsite\tnearest_N_dist_A\tN_chain\tnearest_any_dist_A\tany_chain\tany_atom\tflag\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")

# --- verdict ---
# cha_c derives from source TSA L211 (step 01b). Find that raw TSA's Arg63-N dist.
src_site = "L211"
raw_L = next((r for r in rows if r[0]=="raw_2cht" and r[1]==src_site), None)
placed_C = next((r for r in rows if r[0]=="placed" and r[1]=="C201"), None)
print("\nverdict:")
if raw_L:
    print(f"  raw 2CHT source site {src_site}: nearest Arg63-N = {raw_L[2]} A ({raw_L[7] or 'normal'})")
if placed_C:
    print(f"  placed cha_c (C201):        nearest Arg63-N = {placed_C[2]} A ({placed_C[7] or 'normal'})")
if raw_L and placed_C and raw_L[7] == "SWUNG":
    print("  => Arg63 swung-away geometry is PRESENT IN RAW 2CHT at the cha_c source site.")
    print("     The asymmetry is CRYSTALLOGRAPHIC (inherited), not introduced by our workflow.")
elif raw_L and placed_C and raw_L[7] != "SWUNG" and placed_C[7] == "SWUNG":
    print("  => raw source site is NORMAL but placed cha_c is SWUNG: asymmetry may be")
    print("     INTRODUCED by repair/placement - investigate before freezing.")
print(f"\nWROTE {out}")
PY
echo "STEP 05 DONE"

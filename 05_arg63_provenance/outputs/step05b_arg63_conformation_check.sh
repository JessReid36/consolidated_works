#!/usr/bin/env bash
# Step 05b - confirm the mechanism behind the Arg63' contact pattern.
#
# Hypothesis (from step 05): sites A/B show short Arg63' contacts and site C shows
# a swung rotamer because the A/B/C-frame active sites are completed by DIFFERENT
# physical Arg63 copies than the raw J/K/L source sites, and those A/B/C-chain
# Arg63 side chains simply have their own (crystallographic) conformations.
#
# Test: each placed Arg63' is a specific physical residue from raw 2CHT (the same
# chain letter, since our repaired protein keeps 2CHT A/B/C coordinates for core
# residues). Superpose each placed Arg63 on its OWN backbone (N,CA,C) against the
# same residue in raw 2CHT, then RMSD the side-chain atoms (CB,CG,CD,NE,CZ,NH1,NH2).
# ~0 RMSD => side chain unchanged by our workflow (conformation is crystallographic).
# Large RMSD => repair/cleanup moved the side chain (must investigate before freeze).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
raw2cht="$root/01_inputs/structures/2cht_raw.pdb"
placed="$root/02_preparation/protein_only/abc_repaired_clean.pdb"
admin="$root/00_admin"

python3 - "$raw2cht" "$placed" "$admin" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict
import numpy as np

raw2cht, placed, admin = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out = admin / "step05b_arg63_conformation_report.tsv"

backbone = ["N", "CA", "C"]
sidechain = ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"]

def read_arg63(path):
    # chain -> {atomname: xyz} for ARG 63 (heavy atoms only)
    d = defaultdict(dict)
    for line in path.read_text().splitlines():
        if line.startswith("ATOM  ") and line[17:20].strip() == "ARG" and int(line[22:26]) == 63:
            name = line[12:16].strip()
            elem = line[76:78].strip() if len(line) >= 78 else ""
            if elem.upper() == "H" or name.startswith("H"):
                continue
            d[line[21]][name] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return d

raw = read_arg63(raw2cht)
plc = read_arg63(placed)

def kabsch_fit(P, Q):
    # returns rotation+translation that best fits P onto Q, and applies to P
    Pc, Qc = P.mean(0), Q.mean(0)
    C = (P - Pc).T @ (Q - Qc)
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1, 1, d]) @ Wt
    return Pc, Qc, U

def rmsd(A, B):
    return math.sqrt(((A - B) ** 2).sum() / len(A))

print("Arg63 side-chain conformation: placed (repaired A/B/C) vs raw 2CHT, same chain")
print("(superposed on backbone N,CA,C; RMSD over CB,CG,CD,NE,CZ,NH1,NH2)")
print(f"  {'chain':6s} {'bb_atoms':>8s} {'sc_atoms':>8s} {'sc_RMSD_A':>10s}  {'verdict'}")
rows = []
for ch in sorted(plc):
    if ch not in raw:
        print(f"  {ch:6s}  no raw Arg63 on this chain"); continue
    bb = [a for a in backbone if a in plc[ch] and a in raw[ch]]
    sc = [a for a in sidechain if a in plc[ch] and a in raw[ch]]
    if len(bb) < 3:
        print(f"  {ch:6s}  insufficient backbone atoms"); continue
    P = np.array([plc[ch][a] for a in bb]); Q = np.array([raw[ch][a] for a in bb])
    Pc, Qc, U = kabsch_fit(P, Q)
    # apply backbone-derived transform to placed side chain, compare to raw side chain
    Psc = np.array([plc[ch][a] for a in sc])
    Qsc = np.array([raw[ch][a] for a in sc])
    Psc_fit = (Psc - Pc) @ U + Qc
    r = rmsd(Psc_fit, Qsc)
    verdict = "unchanged" if r < 0.5 else ("shifted" if r < 1.5 else "ROTAMER CHANGE")
    print(f"  {ch:6s} {len(bb):>8d} {len(sc):>8d} {r:>10.3f}  {verdict}")
    rows.append((ch, len(bb), len(sc), f"{r:.3f}", verdict))

with out.open("w") as f:
    f.write("chain\tbackbone_atoms\tsidechain_atoms\tsidechain_rmsd_A\tverdict\n")
    for r in rows:
        f.write("\t".join(str(x) for x in r) + "\n")

print("\ninterpretation:")
print("  all 'unchanged' (<0.5 A) => repair/cleanup did NOT move Arg63 side chains;")
print("    the A/B/C contact pattern is entirely crystallographic (re-partnering).")
print("  any 'ROTAMER CHANGE' => our workflow altered that Arg63 - investigate before freeze.")
print(f"\nWROTE {out}")
PY
echo "STEP 05b DONE"

#!/usr/bin/env bash
# Step 07a - prepare the protein-only PDB for submission to the H++ server.
#
# H++ protonates titratable PROTEIN residues; chorismate is excluded (it is
# separately parameterized by AM1-BCC/GAFF in step 08 and rejoins at the tleap
# build in step 09). The input is the cleaned, repaired heavy-atom protein from
# step 03. Clean TER records are written between chains so the returned file can
# be re-chained if H++ strips chain IDs (a known H++ behaviour).
#
# Server settings to use manually (matching physiological choice + original run):
#   pH 7.0, salinity 0.15 M, internal dielectric 10, external dielectric 80.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
protein="$root/02_preparation/protein_only/abc_repaired_clean.pdb"
hppdir="$root/03_amber/protonation_hpp"
admin="$root/00_admin"
mkdir -p "$hppdir" "$admin"

python3 - "$protein" "$hppdir" "$admin" <<'PY'
import sys
from pathlib import Path
from collections import defaultdict

protein, hppdir, admin = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
out = hppdir / "hpp_input_protein_only_with_TER.pdb"
report = admin / "step07a_hpp_input_report.tsv"

atoms = []
for line in protein.read_text().splitlines():
    if line.startswith("ATOM  "):
        atoms.append(line)
if not atoms:
    sys.exit("FAIL no ATOM records in protein input")

# group by chain to write clean TER records between chains
by_chain = defaultdict(list)
order = []
for line in atoms:
    ch = line[21]
    if ch not in by_chain:
        order.append(ch)
    by_chain[ch].append(line)

# sanity: expect A, B, C only, protein residues 1-127 each
counts = {}
for ch in order:
    resids = sorted(set(int(l[22:26]) for l in by_chain[ch]))
    counts[ch] = (len(by_chain[ch]), min(resids), max(resids), len(resids))

serial = 1
with out.open("w") as f:
    for ch in order:
        for line in by_chain[ch]:
            # renumber atom serials cleanly, keep everything else verbatim
            f.write(f"{line[:6]}{serial:5d}{line[11:]}\n")
            serial += 1
        f.write("TER\n")
    f.write("END\n")

with report.open("w") as f:
    f.write("chain\tatoms\tmin_resid\tmax_resid\tn_residues\n")
    for ch in order:
        a, mn, mx, n = counts[ch]
        f.write(f"{ch}\t{a}\t{mn}\t{mx}\t{n}\n")

print("H++ input prepared (protein only):")
for ch in order:
    a, mn, mx, n = counts[ch]
    print(f"  chain {ch}: {a} atoms, residues {mn}-{mx} ({n} residues)")
total = sum(counts[ch][0] for ch in order)
print(f"  total: {total} heavy atoms, {len(order)} chains, {len(order)} TER records")
if order != ["A","B","C"]:
    print(f"  WARNING unexpected chain set: {order}")
print(f"\nWROTE {out}")
PY
echo "STEP 07a DONE"

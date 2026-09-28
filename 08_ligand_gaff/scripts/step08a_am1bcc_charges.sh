#!/usr/bin/env bash
# Step 08a - derive AM1-BCC charges for chorismate on a local Ubuntu workstation.
#
# WHY OFF-HPC: AM1-BCC needs `sqm`, which cannot run on the HPC (its AMBER22 sqm
# is missing libopenblas.so.0 and no BLAS module fixed it). We derive the charges
# once here with a self-contained conda AmberTools (ships a working BLAS), then
# carry only the charge file to the HPC for GAFF typing in step 08b. Atom order is
# preserved, so charges_am1bcc.dat lines up 1:1 with the canonical CHA template.
#
# Charge model: AM1-BCC - the standard Antechamber/GAFF scheme and the method
# implied by the reference protocol ("Antechamber + GAFF").
# ==============================================================================
# RUN LOCATION: LOCAL Ubuntu workstation - NOT the HPC.
# AM1-BCC needs sqm, which is broken on the HPC; this half is done on your PC with
# conda AmberTools. The HPC half of step 08 is step08b_ligand_gaff.sh.
# ==============================================================================
set -euo pipefail

# --- set these to match your one-time install ----------------------------------
MAMBA="${MAMBA:-$HOME/bin/micromamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
ENV="${ENV:-ambertools}"
# -------------------------------------------------------------------------------

# fail-safe: needs the local conda AmberTools; refuse with guidance if absent
# (e.g. if this was launched on the HPC by mistake).
if ! command -v "$MAMBA" >/dev/null 2>&1; then
  echo "STOP: step 08a runs on your LOCAL Ubuntu workstation, not the HPC."
  echo "  It needs conda AmberTools (working sqm) at: $MAMBA"
  echo "  On the HPC?  Wrong machine - run 08a on your PC, scp charges_am1bcc.dat"
  echo "              up, then run step08b_ligand_gaff.sh on the HPC."
  echo "  On your PC?  Run the one-time micromamba/ambertools setup first."
  exit 1
fi

workdir="$HOME/system_dev_offline/step08a_am1bcc"
input="$workdir/cha_a.mol2"          # scp this canonical CHA template from the HPC first
mkdir -p "$workdir"; cd "$workdir"

[[ -s "$input" ]] || { echo "FAIL missing $input"; echo "  scp it first, e.g.:"; \
  echo "  scp 18660916@hpc1.sun.ac.za:~/system_development/02_preparation/accepted_preprotonation/cha_a.mol2 $workdir/"; exit 1; }

echo "=== step 08a: verify canonical CHA input (24 atoms, 24 bonds, unique names, net -2) ==="
python3 - "$input" <<'PY'
import sys
from collections import Counter
atoms=[]; bonds=0; sec=None
for l in open(sys.argv[1]):
    if l.startswith("@<TRIPOS>"): sec=l.strip(); continue
    if sec=="@<TRIPOS>ATOM" and len(l.split())>=9:
        p=l.split(); atoms.append((p[1], float(p[8])))
    elif sec=="@<TRIPOS>BOND" and l.split() and len(l.split())>=4:
        bonds+=1
if len(atoms)!=24: sys.exit(f"FAIL expected 24 atoms, found {len(atoms)}")
if bonds!=24:      sys.exit(f"FAIL expected 24 bonds, found {bonds}")
names=[a[0] for a in atoms]
dup=[n for n,c in Counter(names).items() if c>1]
if dup: sys.exit(f"FAIL duplicate atom names (need unique CHA template): {dup}")
net=sum(a[1] for a in atoms)
if abs(net+2.0)>1e-4: sys.exit(f"FAIL expected net charge -2, found {net:.4f}")
print(f"PASS 24 atoms, 24 bonds, unique names, net charge {net:+.4f}")
PY

echo
echo "=== step 08a: run AM1-BCC via Antechamber (this calls sqm) ==="
rm -f cha_am1bcc.mol2 antechamber_am1bcc.log sqm.in sqm.out sqm.pdb ANTECHAMBER* ATOMTYPE.INF NEWPDB.PDB PREP.INF
if ! "$MAMBA" run -n "$ENV" antechamber \
      -i "$input" -fi mol2 \
      -o cha_am1bcc.mol2 -fo mol2 \
      -at gaff -c bcc -nc -2 -rn CHA \
      > antechamber_am1bcc.log 2>&1; then
  echo "FAIL antechamber -c bcc failed"; echo "--- log tail ---"; tail -20 antechamber_am1bcc.log; exit 1
fi
[[ -s cha_am1bcc.mol2 ]] || { echo "FAIL no cha_am1bcc.mol2 (sqm may have failed)"; tail -20 antechamber_am1bcc.log; exit 1; }
echo "PASS AM1-BCC produced cha_am1bcc.mol2"

echo
echo "=== step 08a: extract charges, renormalise to exactly -2, write charge file ==="
amber_ver="$("$MAMBA" list -n "$ENV" ambertools 2>/dev/null | awk '$1=="ambertools"{print $2}' | head -1)"
amber_ver="${amber_ver:-unknown}"
python3 - "$input" cha_am1bcc.mol2 "$amber_ver" <<'PY'
import sys, hashlib
from pathlib import Path

inp, bcc, amber_ver = sys.argv[1], sys.argv[2], sys.argv[3]

def read(path):
    atoms=[]; sec=None
    for l in open(path):
        if l.startswith("@<TRIPOS>"): sec=l.strip(); continue
        if sec=="@<TRIPOS>ATOM" and len(l.split())>=9:
            p=l.split(); atoms.append({"name":p[1],"type":p[5],"q":float(p[8])})
    return atoms

src = read(inp)
out = read(bcc)
if len(out)!=24: sys.exit(f"FAIL AM1-BCC mol2 has {len(out)} atoms, expected 24")

q = [a["q"] for a in out]
raw_sum = sum(q)

# renormalise to exactly -2 by distributing the tiny rounding residual over all atoms
residual = -2.0 - raw_sum
q = [x + residual/len(q) for x in q]
assert abs(sum(q)+2.0) < 1e-9, "renormalisation failed"

types = sorted(set(a["type"] for a in out))
dot = [t for t in types if "." in t]
if dot: sys.exit(f"FAIL Tripos dot-types remain: {dot}")

Path("charges_am1bcc.dat").write_text("".join(f"{x: .8f}\n" for x in q))

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

Path("step08a_provenance.txt").write_text(
f"""Step 08a - AM1-BCC charge derivation (off-HPC)

Reason for off-HPC derivation:
  The HPC AMBER22 sqm cannot run (missing libopenblas.so.0; no BLAS module fixed
  it). AM1-BCC requires sqm, so charges were derived on a local Ubuntu workstation
  with a self-contained conda AmberTools that ships a working BLAS.

Method:
  antechamber -i cha_a.mol2 -fi mol2 -o cha_am1bcc.mol2 -fo mol2 \\
              -at gaff -c bcc -nc -2 -rn CHA
  AmberTools: {amber_ver}
  Charges extracted in canonical atom order; renormalised to exactly -2 by
  distributing the {residual:+.6f} e mol2 rounding residual evenly over 24 atoms.

Charge summary:
  raw AM1-BCC sum (mol2, 4 dp) = {raw_sum:+.6f}
  renormalised sum             = {sum(q):+.6f}
  GAFF atom types              = {','.join(types)}

Inputs / outputs (sha256):
  cha_a.mol2         {sha(inp)}
  cha_am1bcc.mol2    {sha(bcc)}
  charges_am1bcc.dat {sha('charges_am1bcc.dat')}

Next: scp charges_am1bcc.dat to the HPC; step 08b runs
  antechamber -at gaff -c rc -cf charges_am1bcc.dat -nc -2 (no sqm needed).
""")

print(f"PASS wrote charges_am1bcc.dat (24 charges, sum {sum(q):+.6f})")
print(f"  raw AM1-BCC sum was {raw_sum:+.6f}; residual {residual:+.6f} distributed")
print(f"  GAFF atom types: {','.join(types)}")
print(f"  max |charge|: {max(abs(x) for x in q):.4f} e")
PY

echo
echo "=== step 08a: outputs ==="
for f in cha_am1bcc.mol2 charges_am1bcc.dat step08a_provenance.txt; do
  [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }
  echo "PASS $workdir/$f"
done
sha256sum "$input" cha_am1bcc.mol2 charges_am1bcc.dat > sha256_step08a.txt
echo
echo "NEXT: scp charges_am1bcc.dat to the HPC, then run step 08b there."
echo "STEP 08a DONE"

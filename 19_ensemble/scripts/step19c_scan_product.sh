#!/usr/bin/env bash
#
# step19c_scan_product.sh - prepare the restrained scan and the free product
# optimisation for a frame whose reactant optimisation has converged.
#
# This completes the per-frame chain that frame 820 followed:
#
#   reactant opt (19a)  ->  restrained scan (here)  ->  product opt (here)  ->  NEB (19a)
#
# The scan exists only to produce a product-side starting guess. Frame 820's
# product_start.pdb was verified byte-identical to its final scan window
# (win_20.pdb), so the chain is scan -> last window -> free optimisation, and
# that is what is reproduced here.
#
# WHY TARGETS ARE REGENERATED PER FRAME
# -------------------------------------
# Frame 820's targets.txt began at 1.464 / 3.251 A, which are the distances of
# its OPTIMISED reactant, not of the raw MD frame. Each frame relaxes to its own
# reactant geometry, so reusing frame 820's targets would start every scan with a
# restraint that fights the structure it is applied to. The first window's target
# is therefore taken from this frame's own optimised reactant, and the windows
# interpolate linearly to the same product-side endpoint (2.400 / 1.540), so the
# path spans the same chemistry in every frame.
#
# Colvar atoms are 0-based absolute indices in the 55,680-atom system and are the
# same in every frame because the topology is identical:
#   break  6215 - 6214  = C4 - O3
#   form   6219 - 6207  = C6 - C1
# These match 00_admin/reacting_atoms.tsv offsets against QMAtoms {6207:6230}.
#
# USAGE
#   bash step19c_scan_product.sh 2450          # one frame
#   bash step19c_scan_product.sh 2450 4085     # several
#   bash step19c_scan_product.sh               # every frame with a converged reactant
#
# Prepares only. Submit what it prints.

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

ORCA=/home/apps2/ORCA/6.0.1
MPI=/apps/openmpi/4.1.1
BLAS=/apps/mambaforge/envs/medaka/lib
NPROC=8
EMAIL=18660916@sun.ac.za

root="$HOME/system_development"
ensdir="$root/05_qmmm/19_ensemble"
ref="$root/05_qmmm/16_scan"
active="$root/05_qmmm/15_active_region/active_atoms_R12.txt"
bridge="$root/05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"

# product-side endpoint of the scan, from frame 820's targets.txt final line
END_BRK=2.400
END_FRM=1.540
NWIN=20

for f in "$active" "$bridge" "$ref/run_scan.sh"; do
  [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }
done

if [[ $# -gt 0 ]]; then
  FRAMES=("$@")
else
  FRAMES=()
  for d in "$ensdir"/frame_*; do
    if grep -q "The minimization has converged" "$d/reactant_opt.out" 2>/dev/null; then
      FRAMES+=("$(basename "$d" | sed 's/frame_0*//')")
    fi
  done
fi
[[ ${#FRAMES[@]} -gt 0 ]] || { echo "no frames with a converged reactant"; exit 0; }
echo "frames: ${FRAMES[*]}"
echo

ACT="$(cat "$active")"
SUBMIT=""

for FR in "${FRAMES[@]}"; do
  PAD=$(printf "%05d" "$FR")
  work="$ensdir/frame_${PAD}"
  rpdb="$work/reactant_opt.pdb"

  if ! grep -q "The minimization has converged" "$work/reactant_opt.out" 2>/dev/null; then
    echo "  SKIP frame $FR - reactant not converged"; continue
  fi
  [[ -s "$rpdb" ]] || { echo "  SKIP frame $FR - no reactant_opt.pdb"; continue; }

  scan="$work/scan"; mkdir -p "$scan"
  cp "$rpdb" "$scan/win_00.pdb"
  cp "$bridge" "$scan/complex_solvated.ORCAFF.prms"
  cp "$ref/run_scan.sh" "$scan/run_scan.sh"
  cp "$ref/stage_tsguess.sh" "$scan/stage_tsguess.sh"
  echo "$ACT" > "$scan/active.txt"

  # --- targets from THIS frame's optimised reactant
  python3 - "$scan/win_00.pdb" "$scan/targets.txt" "$END_BRK" "$END_FRM" "$NWIN" <<'PY'
import sys, math
pdb, out, endb, endf, n = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
xyz = []
for ln in open(pdb):
    if ln[:6] in ("ATOM  ", "HETATM"):
        xyz.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
def d(i, j):
    # explicit rather than math.dist: the cluster python3 predates 3.8
    a, b = xyz[i], xyz[j]
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
brk0 = d(6215, 6214)          # C4 - O3
frm0 = d(6219, 6207)          # C6 - C1
print(f"    optimised reactant: break {brk0:.3f}  form {frm0:.3f} A")
if not (1.35 < brk0 < 1.70 and 2.6 < frm0 < 4.0):
    sys.exit(f"FAIL: reactant distances {brk0:.3f}/{frm0:.3f} are outside the "
             f"expected reactant range - check the optimisation before scanning")
with open(out, "w") as fh:
    # Frame 820's targets.txt begins AT the optimised reactant distances and
    # steps to the product-side endpoint in n-1 equal intervals, so the first
    # window restrains the structure where it already sits. Reproduced here.
    for k in range(n):
        t = k / (n - 1)
        fh.write(f"{brk0 + t*(endb-brk0):.3f} {frm0 + t*(endf-frm0):.3f}\n")
print(f"    wrote {n} targets, {brk0:.3f}/{frm0:.3f} -> {endb:.3f}/{endf:.3f}")
PY

  # --- scan job
  {
    echo '#!/bin/bash'
    echo "#PBS -N cm19_s${PAD}"
    echo "#PBS -l select=1:ncpus=$NPROC:mem=16gb"
    echo '#PBS -l walltime=176:00:00'
    echo '#PBS -m ae'; echo "#PBS -M $EMAIL"; echo '#PBS -j oe'
    echo "#PBS -o $scan/scan.pbs.out"
    echo "cd $scan"
    echo '# RUNLOCK: mkdir is atomic - a second job in this directory exits'
    echo 'if ! mkdir .running 2>/dev/null; then echo "FAIL: another job holds .running here"; exit 1; fi'
    echo 'trap "rmdir .running 2>/dev/null" EXIT'
    echo "export PATH=\"$MPI/bin:\$PATH\""
    echo "export LD_LIBRARY_PATH=\"$ORCA/lib:$MPI/lib:$BLAS:\$LD_LIBRARY_PATH\""
    echo 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1'
    echo 'echo "host=$(hostname) start=$(date)"'
    echo 'bash run_scan.sh'
    echo 'RC=$?'
    echo 'echo "end=$(date) rc=$RC"'
    echo '# the scan is only useful if it reached the final window'
    echo "if [ -s win_$(printf '%02d' $NWIN).pdb ]; then echo SCAN_PASS; else echo SCAN_INCOMPLETE; exit 1; fi"
    echo 'bash stage_tsguess.sh || exit 1'
    echo 'tail -3 scan_progress.tsv 2>/dev/null'
  } > "$work/scan.pbs"

  # --- product optimisation, identical to step 17 but from this frame's win_20
  {
    echo "! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX L-Opt"
    echo "%maxcore 3000"
    echo "%pal nprocs $NPROC end"
    echo "%scf MaxIter 200 end"
    echo "%geom MaxIter 2000 end"
    echo "%qmmm"
    echo "  QMAtoms {6207:6230} end"
    echo "  ActiveAtoms {$ACT} end"
    echo '  ORCAFFFilename "complex_solvated.ORCAFF.prms"'
    echo "end"
    echo "*pdbfile -2 1 product_start.pdb"
  } > "$work/product_opt.inp"

  {
    echo '#!/bin/bash'
    echo "#PBS -N cm19_p${PAD}"
    echo "#PBS -l select=1:ncpus=$NPROC:mem=16gb"
    echo '#PBS -l walltime=176:00:00'
    echo '#PBS -m ae'; echo "#PBS -M $EMAIL"; echo '#PBS -j oe'
    echo "#PBS -o $work/product_opt.pbs.out"
    echo "cd $work"
    echo "# the product guess is the final scan window, as it was for frame 820"
    echo "cp scan/win_$(printf '%02d' $NWIN).pdb product_start.pdb || { echo 'FAIL: scan did not reach the final window'; exit 1; }"
    echo "export PATH=\"$MPI/bin:\$PATH\""
    echo "export LD_LIBRARY_PATH=\"$ORCA/lib:$MPI/lib:$BLAS:\$LD_LIBRARY_PATH\""
    echo 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1'
    echo 'echo "host=$(hostname) start=$(date)"'
    echo '( while sleep 1800; do printf "  [%s] steps=%s  lastE=%s\n" "$(date +%H:%M)" "$(wc -l < product_opt-minimize-ener.csv 2>/dev/null)" "$(tail -1 product_opt-minimize-ener.csv 2>/dev/null | cut -d, -f2)"; done ) & HB=$!'
    echo "$ORCA/orca product_opt.inp > product_opt.out 2>&1"
    echo 'RC=$?'
    echo 'kill $HB 2>/dev/null'
    echo 'echo "orca_exit=$RC end=$(date)"'
    echo 'grep -q "ORCA TERMINATED NORMALLY" product_opt.out && grep -q "The minimization has converged" product_opt.out && echo PRODUCT_OPT_PASS || echo PRODUCT_OPT_INCOMPLETE'
    echo '# stage the endpoints under the names the NEB expects'
    echo 'cp reactant_opt.pdb reactant.pdb 2>/dev/null'
    echo 'cp product_opt.pdb  product.pdb  2>/dev/null'
  } > "$work/product_opt.pbs"

  echo "  prepared frame $FR"
  SUBMIT="$SUBMIT  qsub $work/scan.pbs      # then, when SCAN_PASS:\n  qsub $work/product_opt.pbs\n"
done

echo
echo "STEP 19c prepared. Submit the scans first:"
echo -e "$SUBMIT"
echo "product_opt.pbs stages reactant.pdb and product.pdb for the NEB on success,"
echo "so after it reports PRODUCT_OPT_PASS you can submit neb.pbs directly."

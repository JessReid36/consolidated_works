#!/usr/bin/env bash
# Step 13b - minimal QM/MM single-point smoke test (HF-3c EnGrad, NON-production) on one
# selected frame, validating bridge + coords + substrate-only QM region + electrostatic
# embedding terminate normally. frame_00820_CHA2, QM = 24 CHA atoms {6207:6230} (0-based),
# charge -2, mult 1. Single-process. All tools pinned to 1 OpenBLAS thread (login-node cap).
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
BLAS=/apps/mambaforge/envs/medaka/lib
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
bridge="$root/05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"
rst="$root/05_qmmm/12_frame_selection/frames/frame_00820_CHA2.rst7"
work="$root/05_qmmm/13_smoketest"
email="18660916@sun.ac.za"
for f in "$prmtop" "$bridge" "$rst"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done
[[ -e "$BLAS/libopenblas.so.0" ]] || { echo "FAIL openblas missing at $BLAS"; exit 1; }
mkdir -p "$work"; cd "$work"

echo "=== rst7 -> PDB (ambpdb, single-threaded) ==="
set +u; module load app/amber22/22; set -u
export LD_LIBRARY_PATH="$BLAS:${LD_LIBRARY_PATH:-}"
ambpdb -p "$prmtop" -c "$rst" > frame_00820_CHA2.pdb 2>ambpdb.err || { echo "FAIL ambpdb"; cat ambpdb.err; exit 1; }
n=$(grep -cE '^(ATOM|HETATM)' frame_00820_CHA2.pdb); echo "  PDB atoms: $n (expect 55680)"
[[ "$n" -eq 55680 ]] || { echo "FAIL pdb atoms $n"; exit 1; }

echo "=== bridge + input ==="
cp "$bridge" ./complex_solvated.ORCAFF.prms
printf '%s\n' \
  '! QMMM HF-3c EnGrad' \
  '%maxcore 2000' \
  '%pal nprocs 1 end' \
  '%qmmm' \
  '  QMAtoms {6207:6230} end' \
  '  ORCAFFFilename "complex_solvated.ORCAFF.prms"' \
  'end' \
  '*pdbfile -2 1 frame_00820_CHA2.pdb' > smoke.inp
cat smoke.inp

echo "=== PBS ==="
printf '%s\n' \
  '#!/bin/bash' \
  '#PBS -N cm13b_smoke' \
  '#PBS -l select=1:ncpus=1:mem=8gb' \
  '#PBS -l walltime=02:00:00' \
  "#PBS -m ae" \
  "#PBS -M $email" \
  '#PBS -j oe' \
  "#PBS -o $work/smoke.pbs.out" \
  "cd $work" \
  "export LD_LIBRARY_PATH=\"$ORCA/lib:$BLAS:\$LD_LIBRARY_PATH\"" \
  'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1' \
  'echo "host=$(hostname) start=$(date)"' \
  "$ORCA/orca smoke.inp > smoke.out 2>&1" \
  'echo "orca_exit=$? end=$(date)"' \
  'grep -iE "FINAL SINGLE POINT ENERGY|ORCA TERMINATED NORMALLY" smoke.out | tail -5' \
  'grep -q "ORCA TERMINATED NORMALLY" smoke.out && echo SMOKE_TEST_PASS || echo SMOKE_TEST_FAIL' \
  > smoke.pbs

jobid="$(qsub smoke.pbs)"; echo "$jobid" > smoke_jobid.txt
echo "PASS submitted: $jobid"

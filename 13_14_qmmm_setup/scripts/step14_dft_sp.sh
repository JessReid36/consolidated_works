#!/usr/bin/env bash
# Step 14 - production-level QM/MM single-point (validates the DFT level runs on our
# frame; nothing moves, so active-region decision not yet needed). B3LYP-D3BJ/def2-SVP/
# RIJCOSX EnGrad on frame_00820_CHA2. QM = 24 CHA atoms {6207:6230}, charge -2, mult 1.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
BLAS=/apps/mambaforge/envs/medaka/lib
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
bridge="$root/05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"
rst="$root/05_qmmm/12_frame_selection/frames/frame_00820_CHA2.rst7"
work="$root/05_qmmm/14_dft_sp"
email="18660916@sun.ac.za"
for f in "$prmtop" "$bridge" "$rst"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done
[[ -e "$BLAS/libopenblas.so.0" ]] || { echo "FAIL openblas missing"; exit 1; }
mkdir -p "$work"; cd "$work"

echo "=== rst7 -> PDB ==="
set +u; module load app/amber22/22; set -u
export LD_LIBRARY_PATH="$BLAS:${LD_LIBRARY_PATH:-}"
ambpdb -p "$prmtop" -c "$rst" > frame_00820_CHA2.pdb 2>ambpdb.err || { echo "FAIL ambpdb"; cat ambpdb.err; exit 1; }
n=$(grep -cE '^(ATOM|HETATM)' frame_00820_CHA2.pdb); echo "  PDB atoms: $n (expect 55680)"
[[ "$n" -eq 55680 ]] || { echo "FAIL pdb atoms $n"; exit 1; }

cp "$bridge" ./complex_solvated.ORCAFF.prms
printf '%s\n' \
  '! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX EnGrad' \
  '%maxcore 3000' \
  '%pal nprocs 1 end' \
  '%qmmm' \
  '  QMAtoms {6207:6230} end' \
  '  ORCAFFFilename "complex_solvated.ORCAFF.prms"' \
  'end' \
  '*pdbfile -2 1 frame_00820_CHA2.pdb' > dft_sp.inp
echo "--- dft_sp.inp ---"; cat dft_sp.inp

printf '%s\n' \
  '#!/bin/bash' \
  '#PBS -N cm14_dftsp' \
  '#PBS -l select=1:ncpus=1:mem=16gb' \
  '#PBS -l walltime=04:00:00' \
  '#PBS -m ae' \
  "#PBS -M $email" \
  '#PBS -j oe' \
  "#PBS -o $work/dft_sp.pbs.out" \
  "cd $work" \
  "export LD_LIBRARY_PATH=\"$ORCA/lib:$BLAS:\$LD_LIBRARY_PATH\"" \
  'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1' \
  'echo "host=$(hostname) start=$(date)"' \
  "$ORCA/orca dft_sp.inp > dft_sp.out 2>&1" \
  'echo "orca_exit=$? end=$(date)"' \
  'grep -iE "FINAL SINGLE POINT ENERGY|ORCA TERMINATED NORMALLY" dft_sp.out | tail -5' \
  'grep -q "ORCA TERMINATED NORMALLY" dft_sp.out && echo DFT_SP_PASS || echo DFT_SP_FAIL' \
  > dft_sp.pbs

jobid="$(qsub dft_sp.pbs)"; echo "$jobid" > dft_sp_jobid.txt
echo "PASS submitted: $jobid"

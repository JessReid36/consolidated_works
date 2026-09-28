#!/usr/bin/env bash
# Step 15b - QM/MM reactant optimisation. B3LYP-D3BJ/def2-SVP/RIJCOSX with the L-Opt
# (Cartesian L-BFGS via orca_md) optimiser - REQUIRED for the ~1959-atom active region
# (plain internal-coordinate Opt is ~6h/cycle and never converges; L-Opt ~5.5h total).
# QM = 24 CHA atoms {6207:6230}; ActiveAtoms = Claeyssens-style free region from 15a.
# Parallel: ORCA 6.0.1 + OpenMPI 4.1.1. NOTE: on ORCA 6.1+, L-Opt is renamed MD-L-Opt.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1
MPI=/apps/openmpi/4.1.1
BLAS=/apps/mambaforge/envs/medaka/lib
NPROC=8
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
bridge="$root/05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"
rst="$root/05_qmmm/12_frame_selection/frames/frame_00820_CHA2.rst7"
active="$root/05_qmmm/15_active_region/active_atoms_R12.txt"
work="$root/05_qmmm/15_reactant_opt_lopth"
email="18660916@sun.ac.za"
for f in "$prmtop" "$bridge" "$rst" "$active"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done
[[ -x "$MPI/bin/mpirun" ]] || { echo "FAIL mpirun missing"; exit 1; }
mkdir -p "$work"; cd "$work"

echo "=== rst7 -> PDB ==="
set +u; module load app/amber22/22; set -u
export LD_LIBRARY_PATH="$BLAS:${LD_LIBRARY_PATH:-}"
ambpdb -p "$prmtop" -c "$rst" > frame_00820_CHA2.pdb 2>ambpdb.err || { echo "FAIL ambpdb"; cat ambpdb.err; exit 1; }
n=$(grep -cE '^(ATOM|HETATM)' frame_00820_CHA2.pdb); [[ "$n" -eq 55680 ]] || { echo "FAIL pdb atoms $n"; exit 1; }
echo "  PDB atoms: $n"

cp "$bridge" ./complex_solvated.ORCAFF.prms
ACT="$(cat "$active")"; echo "  active atoms: $(wc -w < "$active")"

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
  echo "*pdbfile -2 1 frame_00820_CHA2.pdb"
} > reactant_opt.inp

{
  echo '#!/bin/bash'
  echo '#PBS -N cm15b_ropt'
  echo "#PBS -l select=1:ncpus=$NPROC:mem=120gb"
  echo '#PBS -l walltime=168:00:00'
  echo '#PBS -m ae'; echo "#PBS -M $email"; echo '#PBS -j oe'
  echo "#PBS -o $work/reactant_opt.pbs.out"
  echo "cd $work"
  echo "export PATH=\"$MPI/bin:\$PATH\""
  echo "export LD_LIBRARY_PATH=\"$ORCA/lib:$MPI/lib:$BLAS:\$LD_LIBRARY_PATH\""
  echo 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1'
  echo 'echo "host=$(hostname) start=$(date)"'
  echo "$ORCA/orca reactant_opt.inp > reactant_opt.out 2>&1"
  echo 'echo "orca_exit=$? end=$(date)"'
  echo 'grep -iE "The minimization has converged|ORCA TERMINATED NORMALLY" reactant_opt.out | tail -4'
  echo 'grep -q "ORCA TERMINATED NORMALLY" reactant_opt.out && grep -q "The minimization has converged" reactant_opt.out && echo REACTANT_OPT_PASS || echo REACTANT_OPT_INCOMPLETE'
  echo '# cleanup: keep pdb/mdrestart/out, drop the bulky trajectory'
  echo 'rm -f reactant_opt_trj.dcd 2>/dev/null'
} > reactant_opt.pbs

echo "STEP 15b prepared. Submit: qsub $work/reactant_opt.pbs"

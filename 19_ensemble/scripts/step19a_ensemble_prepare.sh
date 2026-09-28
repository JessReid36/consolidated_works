#!/usr/bin/env bash
#
# step19a_ensemble_prepare.sh - prepare QM/MM barrier calculations for the remaining
# frames of the step-12a selection, so the Phase-1 barrier can be reported as an
# ensemble mean rather than a single path.
#
# WHY
# ---
# The literature value this work is compared against (Claeyssens et al. 2011,
# 11.3 kcal/mol for BsCM) is an AVERAGE over multiple QM/MM pathways. A single
# path is not the observable. Frame 820 gives +16.00 kcal/mol, +3.3 above the
# experimental dH‡ of 12.7, and without a spread there is no way to tell whether
# that is the model's central value or a tail.
#
# The step-12a selection already chose twelve catalytically competent frames for
# exactly this purpose ("Claeyssens-style: N competent reactant frames for
# multiple QM/MM paths -> barrier distribution"). Only frame 820 was carried
# through. This script prepares the other eleven.
#
# WHAT IS AND IS NOT REPEATED
# ---------------------------
# Repeated per frame:   reactant optimisation (L-Opt, large active region)
#                       product optimisation  (L-Opt, same region)
#                       NEB-TS barrier, from the scan maximum as guess
#
# NOT repeated:         reduced-region TS optimisation with Hybrid Hessian
#                       612-displacement endpoint frequencies
#                       IRC in both directions
#
# The omitted steps exist once, on frame 820, and establish that the saddle is a
# genuine first-order transition state connecting the right basins. Their result
# also licenses the shortcut: on frame 820 the converged NEB climbing image gave
# 15.94 kcal/mol against the fully characterised 16.00, agreeing to 0.06. NEB-TS
# is therefore a validated proxy FOR THIS SYSTEM, demonstrated rather than
# assumed, and the ensemble is reported by that proxy.
#
# COST
# ----
# Eleven frames x (two L-Opt optimisations + one 8-image NEB). The frame-820
# reactant optimisation took roughly 5.5 h on 8 cores; budget similarly per
# endpoint and longer for the NEB. Submit in batches rather than all at once.
#
# USAGE
# -----
#   bash step19a_ensemble_prepare.sh            # prepare all remaining frames
#   bash step19a_ensemble_prepare.sh 2450 4085  # prepare specific frames
#
# Then submit what it prints. Nothing is submitted automatically.

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

ORCA=/home/apps2/ORCA/6.0.1
MPI=/apps/openmpi/4.1.1
BLAS=/apps/mambaforge/envs/medaka/lib
NPROC=8
EMAIL=18660916@sun.ac.za

root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
bridge="$root/05_qmmm/13_bridge/complex_solvated.ORCAFF.prms"
seldir="$root/05_qmmm/12_frame_selection"
manifest="$seldir/selection_manifest.tsv"
active="$root/05_qmmm/15_active_region/active_atoms_R12.txt"
ensdir="$root/05_qmmm/19_ensemble"

for f in "$prmtop" "$bridge" "$manifest" "$active"; do
  [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }
done
[[ -x "$MPI/bin/mpirun" ]] || { echo "FAIL mpirun missing"; exit 1; }

# frames to prepare: those given on the command line, else every frame in the
# manifest except 820, which is already done
if [[ $# -gt 0 ]]; then
  FRAMES=("$@")
else
  mapfile -t FRAMES < <(awk 'NR>1 && $2!=820 {print $2}' "$manifest")
fi
echo "frames to prepare: ${FRAMES[*]}"
echo

mkdir -p "$ensdir"
ACT="$(cat "$active")"
NACT=$(wc -w < "$active")
echo "active region: $NACT atoms (same definition as frame 820)"
echo

set +u; module load app/amber22/22; set -u
export LD_LIBRARY_PATH="$BLAS:${LD_LIBRARY_PATH:-}"

SUBMIT=""
for FR in "${FRAMES[@]}"; do
  PAD=$(printf "%05d" "$FR")
  rst="$seldir/frames/frame_${PAD}_CHA2.rst7"
  if [[ ! -s "$rst" ]]; then
    echo "  SKIP frame $FR - $rst not found"
    continue
  fi
  work="$ensdir/frame_${PAD}"
  mkdir -p "$work"; cd "$work"

  # --- rst7 -> pdb, with the same atom-count check as step 15b
  if [[ ! -s "frame_${PAD}_CHA2.pdb" ]]; then
    ambpdb -p "$prmtop" -c "$rst" > "frame_${PAD}_CHA2.pdb" 2>ambpdb.err || {
      echo "  FAIL ambpdb frame $FR"; cat ambpdb.err; continue; }
  fi
  n=$(grep -cE '^(ATOM|HETATM)' "frame_${PAD}_CHA2.pdb")
  [[ "$n" -eq 55680 ]] || { echo "  FAIL frame $FR: pdb has $n atoms, expected 55680"; continue; }

  cp "$bridge" ./complex_solvated.ORCAFF.prms

  # --- reactant optimisation, identical settings to step 15b
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
    echo "*pdbfile -2 1 frame_${PAD}_CHA2.pdb"
  } > reactant_opt.inp

  {
    echo '#!/bin/bash'
    echo "#PBS -N cm19_r${PAD}"
    echo "#PBS -l select=1:ncpus=$NPROC:mem=16gb"
    echo '#PBS -l walltime=176:00:00'
    echo '#PBS -m ae'; echo "#PBS -M $EMAIL"; echo '#PBS -j oe'
    echo "#PBS -o $work/reactant_opt.pbs.out"
    echo "cd $work"
    echo "export PATH=\"$MPI/bin:\$PATH\""
    echo "export LD_LIBRARY_PATH=\"$ORCA/lib:$MPI/lib:$BLAS:\$LD_LIBRARY_PATH\""
    echo 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1'
    echo 'echo "host=$(hostname) start=$(date)"'
    echo "$ORCA/orca reactant_opt.inp > reactant_opt.out 2>&1"
    echo 'echo "orca_exit=$? end=$(date)"'
    echo 'grep -q "ORCA TERMINATED NORMALLY" reactant_opt.out && grep -q "The minimization has converged" reactant_opt.out && echo REACTANT_OPT_PASS || echo REACTANT_OPT_INCOMPLETE'
    echo 'rm -f reactant_opt_trj.dcd 2>/dev/null'
  } > reactant_opt.pbs

  # --- NEB, prepared but only runnable once the endpoints exist
  {
    echo "! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX NEB-TS SlowConv TightSCF"
    echo "%maxcore 3000"
    echo "%pal nprocs $NPROC end"
    echo "%scf MaxIter 250 end"
    echo "%neb"
    echo '  Product_PDBFile "product.pdb"'
    echo '  TS_PDBFile "tsguess.pdb"'
    echo "  NImages 8"
    echo "  Interpolation Linear"
    echo "  Prepare_Frags false"
    echo "end"
    echo "%qmmm"
    echo "  QMAtoms {6207:6230} end"
    echo "  ActiveAtoms {$ACT} end"
    echo '  ORCAFFFilename "complex_solvated.ORCAFF.prms"'
    echo "end"
    echo "*pdbfile -2 1 reactant.pdb"
  } > neb.inp

  {
    echo '#!/bin/bash'
    echo "#PBS -N cm19_n${PAD}"
    echo "#PBS -l select=1:ncpus=$NPROC:mem=250gb"
    echo '#PBS -l walltime=176:00:00'
    echo '#PBS -m ae'; echo "#PBS -M $EMAIL"; echo '#PBS -j oe'
    echo "#PBS -o $work/neb.pbs.out"
    echo "cd $work"
    echo "export PATH=\"$MPI/bin:\$PATH\""
    echo "export LD_LIBRARY_PATH=\"$ORCA/lib:$MPI/lib:$BLAS:\$LD_LIBRARY_PATH\""
    echo 'export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1'
    echo '# the NEB needs both optimised endpoints present under these names'
    echo 'for f in reactant.pdb product.pdb; do [[ -s "$f" ]] || { echo "FAIL missing $f - run the endpoint optimisations first"; exit 1; }; done'
    echo 'echo "host=$(hostname) start=$(date)"'
    echo "$ORCA/orca neb.inp > neb.out 2>&1"
    echo 'echo "orca_exit=$? end=$(date)"'
    echo 'grep -iE "THE NEB OPTIMIZATION HAS CONVERGED|ORCA TERMINATED NORMALLY" neb.out | tail -4'
    echo '# the full-system trajectories are large and the quantum-region'
    echo '# versions of the same paths are kept; the Hessian is kept because'
    echo '# a restart of the transition-state optimisation would need it'
    echo 'rm -f neb_MEP_ALL_trj.xyz neb_MEP_trj.xyz neb_initial_path_trj.xyz neb_MEP_ALL.activeRegion_trj.xyz neb*.tmp neb_im*.gbw 2>/dev/null'
  } > neb.pbs

  echo "  prepared frame $FR -> $work"
  SUBMIT="$SUBMIT  qsub $work/reactant_opt.pbs\n"
done

echo
echo "STEP 19a prepared. Submit the reactant optimisations first:"
echo -e "$SUBMIT"
echo "When each reactant converges, produce its product endpoint the same way"
echo "frame 820 did (step 17), copy the optimised endpoints into the frame"
echo "directory as reactant.pdb and product.pdb, then submit neb.pbs."
echo
echo "Collect with: python3 step19b_collect.py"

#!/usr/bin/env bash
# step18e/18f - reduced-region re-optimisation of the reaction-path minima.
# Brings reactant and product into the SAME 102-atom reduced active region as the
# step18c OptTS, on the shared ts_start.pdb frozen bulk, so R/TS/P form one consistent
# reduced model and the barrier TS-R is a clean same-region number. Minima can afford the
# full active region, the TS could not (Hessian memory) - so consistency is enforced by
# bringing the minima down to the TS region rather than the reverse.
# Usage: step18e_reduced_region_reopt.sh {reactant|product}
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ORCA=/home/apps2/ORCA/6.0.1; MPI=/apps/openmpi/4.1.1; BLAS=/apps/mambaforge/envs/medaka/lib
NPROC=8; EMAIL=18660916@sun.ac.za
state="${1:-}"
case "$state" in
  reactant) src="15_reactant_opt_lopth/reactant_opt.pdb"; workname="18e_reactant_reduced"; stem="reactant_reduced";;
  product)  src="17_product_opt/product_opt.pdb";         workname="18f_product_reduced";  stem="product_reduced";;
  *) echo "usage: $0 {reactant|product}"; exit 2;;
esac
qmmm="$HOME/system_development/05_qmmm"; work="$qmmm/$workname"; mkdir -p "$work"; cd "$work"
active="$qmmm/18c_reduced_region/active_reduced.txt"
bulk="$qmmm/18c_reduced_region/ts_start.pdb"
target="$qmmm/$src"; orcaff="$qmmm/18c_reduced_region/complex_solvated.ORCAFF.prms"
for f in "$active" "$bulk" "$target" "$orcaff"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done
cp "$orcaff" ./
python3 "$(dirname "$0")/step18e_splice_reduced_start.py" "$active" "$bulk" "$target" "${stem}_start.pdb"
ACT="$(cat "$active")"
cat > "${stem}.inp" <<INP
! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX Opt TightSCF
%maxcore 3000
%pal nprocs $NPROC end
%scf MaxIter 250 end
%geom MaxIter 500 end
%qmmm
 QMAtoms {6207:6230} end
 ActiveAtoms {$ACT} end
 ORCAFFFilename "complex_solvated.ORCAFF.prms"
end
*pdbfile -2 1 ${stem}_start.pdb
INP
cat > "${stem}.pbs" <<PBS
#!/bin/bash
#PBS -N cm18e_${state}_red
#PBS -l select=1:ncpus=$NPROC:mem=24gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M $EMAIL
#PBS -j oe
#PBS -o $work/${stem}.pbs.out
cd $work
export PATH="$MPI/bin:\$PATH"
export LD_LIBRARY_PATH="$ORCA/lib:$MPI/lib:$BLAS:\$LD_LIBRARY_PATH"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "host=\$(hostname) start=\$(date)"
$ORCA/orca ${stem}.inp > ${stem}.out 2>&1
echo "orca_exit=\$? end=\$(date)"
grep -q "ORCA TERMINATED NORMALLY" ${stem}.out && grep -q "THE OPTIMIZATION HAS CONVERGED" ${stem}.out && echo ${state}_REDUCED_PASS || echo ${state}_REDUCED_INCOMPLETE
PBS
echo "prepared $work"
echo "submit: qsub -N cm18e_${state}_red $work/${stem}.pbs"

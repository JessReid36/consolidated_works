#!/usr/bin/env bash
# Step 10a - restrained energy minimisation (paper: solute heavy atoms restrained at
# 200 kcal/mol/A^2 while solvent + hydrogens relax). Parallel pmemd.MPI (16 ranks,
# one node) via the cluster's openmpi build. Runs in node-local scratch; copies
# results back. One job, no continuation. Same input/physics as the serial version.
set -euo pipefail

root="$HOME/system_development"
build="$root/03_amber/tleap_build"
prmtop="$build/complex_solvated.prmtop"
inpcrd="$build/complex_solvated.inpcrd"
rundir="$root/04_amber_md/10a_min"
admin="$root/00_admin"
email="18660916@sun.ac.za"
ncpus=16
mkdir -p "$rundir" "$admin"

echo "=== step 10a: input presence ==="
for f in "$prmtop" "$inpcrd"; do
  [[ -e "$f" ]] || { echo "FAIL missing: $f"; exit 1; }
  echo "PASS $f"
done

# Cancel a previously-submitted 10a (e.g. the serial job) so we don't run two.
if [[ -s "$rundir/10a_min_jobid.txt" ]]; then
  old="$(cat "$rundir/10a_min_jobid.txt")"
  echo "=== cancelling previous 10a job $old (if still queued/running) ==="
  qdel "$old" 2>/dev/null && echo "cancelled $old" || echo "(job $old not active - nothing to cancel)"
fi

echo
echo "=== step 10a: write minimisation input ==="
cat > "$rundir/min10a.in" <<'MDIN'
Restrained minimisation: solute heavy atoms at 200 kcal/mol/A^2; solvent + H relax
 &cntrl
  imin=1, ntmin=1, maxcyc=10000, ncyc=5000,
  ntb=1, cut=9.0,
  ntr=1, restraint_wt=200.0, restraintmask='!(:WAT,Na+) & !@H=',
  ntpr=100, ntwx=0, ntxo=1, ioutfm=1,
 /
MDIN
cat "$rundir/min10a.in"

echo
echo "=== step 10a: write PBS job (pmemd.MPI, $ncpus ranks) ==="
pbase="$(basename "$prmtop")"; ibase="$(basename "$inpcrd")"
cat > "$rundir/10a_min.pbs" <<EOF
#!/bin/bash
#PBS -N cm10a_min
#PBS -l select=1:ncpus=$ncpus:mpiprocs=$ncpus:mem=16gb
#PBS -l walltime=02:00:00
#PBS -m ae
#PBS -M $email
#PBS -j oe
#PBS -o $rundir/10a_min.pbs.out

set -uo pipefail
scratch="/scratch-small-local/\${PBS_JOBID}"
mkdir -p "\$scratch" || { echo "FAIL scratch mkdir"; exit 1; }
echo "host=\$(hostname)  jobid=\${PBS_JOBID:-UNSET}  scratch=\$scratch  start=\$(date)"

cp "$rundir/min10a.in" "$prmtop" "$inpcrd" "\$scratch"/
cd "\$scratch"

set +u
export PERL5LIB="\${PERL5LIB:-}" PYTHONPATH="\${PYTHONPATH:-}"
module load app/ambermpi/22mpi
set -u
echo "AMBERHOME=\${AMBERHOME:-UNSET}"
command -v pmemd.MPI >/dev/null || { echo "FAIL pmemd.MPI not found"; exit 1; }
command -v mpirun    >/dev/null || { echo "FAIL mpirun not found"; exit 1; }

np=\$(wc -l < "\${PBS_NODEFILE}")
echo "MPI ranks (np)=\$np  nodefile=\${PBS_NODEFILE}"

echo "=== running restrained minimisation (pmemd.MPI) ==="
mpirun -np \${np} --hostfile "\${PBS_NODEFILE}" pmemd.MPI -O \\
  -i min10a.in -o 10a_min.out \\
  -p "$pbase" -c "$ibase" \\
  -r 10a_min.rst7 -ref "$ibase" -inf 10a_min.mdinfo
rc=\$?

cp -f 10a_min.out 10a_min.mdinfo "$rundir"/ 2>/dev/null || true
[[ \$rc -eq 0 ]] && cp -f 10a_min.rst7 "$rundir"/ 2>/dev/null || true
cd "$rundir"; rm -rf "\$scratch"

if [[ \$rc -ne 0 ]]; then
  echo "FAIL pmemd.MPI rc=\$rc"; tail -30 "$rundir/10a_min.out" 2>/dev/null || true; exit 1
fi
echo "=== minimisation key lines ==="
grep -E "NSTEP|FINAL RESULTS|ENERGY|RMS|RESTRAINT" "$rundir/10a_min.out" | tail -20 || true
[[ -s "$rundir/10a_min.rst7" ]] || { echo "FAIL no rst7 produced"; exit 1; }
echo "end=\$(date)"
echo "FINAL PASS: 10a restrained minimisation complete"
EOF
cat "$rundir/10a_min.pbs"

echo
echo "=== step 10a: submit ==="
jobid="$(qsub "$rundir/10a_min.pbs")"
echo "$jobid" > "$rundir/10a_min_jobid.txt"
cat > "$admin/step10a_submit_summary.txt" <<EOF
step10a restrained minimisation
enginepmemd.MPI ($ncpus ranks, 1 node), AMBER22 (app/ambermpi/22mpi, openmpi/4.1.1)
restraint200 kcal/mol/A^2 on solute heavy atoms (mask '!(:WAT,Na+) & !@H=')
maxcyc/ncyc10000/5000
prmtop$prmtop
inpcrd$inpcrd
rundir$rundir
jobid$jobid
EOF
echo "PASS submitted: $jobid"
echo "monitor: qstat -u \$USER    (or: qstat -f $jobid)"
echo "STEP 10a SUBMITTED (MPI)"

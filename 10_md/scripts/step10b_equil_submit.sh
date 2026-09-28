#!/usr/bin/env bash
# Step 10b - equilibration: five 20 ps NPT Langevin stages (paper), one GPU job.
# Stage 1 heats 0->300 K while solute heavy atoms are restrained; the restraint is
# then relaxed 200 -> 100 -> 50 -> 10 -> 0 kcal/mol/A^2 across the five stages.
# pmemd.cuda, MC barostat, SHAKE dt=0.002, fixed Langevin seed. One submission.
set -euo pipefail
root="$HOME/system_development"
build="$root/03_amber/tleap_build"
prmtop="$build/complex_solvated.prmtop"
minrst="$root/04_amber_md/10a_min/10a_min.rst7"
inpref="$build/complex_solvated.inpcrd"
rundir="$root/04_amber_md/10b_equil"
admin="$root/00_admin"
email="18660916@sun.ac.za"
cuda_lib="/apps/mambaforge/pkgs/cudatoolkit-11.8.0-h37601d7_11/lib"
blas_lib="/apps/mambaforge/envs/medaka/lib"
seed=531984
mkdir -p "$rundir" "$admin"

echo "=== step 10b: input presence ==="
for f in "$prmtop" "$minrst" "$inpref"; do
  [[ -e "$f" ]] || { echo "FAIL missing: $f"; exit 1; }
  echo "PASS $f"
done

echo
echo "=== step 10b: write five NPT mdin stages (200->100->50->10->0) ==="
weights=(200 100 50 10 0)
for i in 1 2 3 4 5; do
  w="${weights[$((i-1))]}"
  if [[ $i -eq 1 ]]; then startln="irest=0, ntx=1, tempi=0.0, temp0=300.0,"; else startln="irest=1, ntx=5, temp0=300.0,"; fi
  if [[ "$w" == "0" ]]; then restrln="ntr=0,"; else restrln="ntr=1, restraint_wt=${w}.0, restraintmask='!(:WAT,Na+) & !@H=',"; fi
  cat > "$rundir/npt${i}.in" <<MDIN
NPT equilibration stage $i: restraint ${w} kcal/mol/A^2, 20 ps, 300 K, 1 atm
 &cntrl
  imin=0, $startln
  ntb=2, ntp=1, barostat=2, pres0=1.0, cut=9.0,
  nstlim=10000, dt=0.002,
  ntc=2, ntf=2,
  ntt=3, gamma_ln=2.0, ig=$seed,
  $restrln
  ntpr=1000, ntwx=0, ntwr=10000, ntxo=1, ioutfm=1, iwrap=1,
 /
MDIN
done
echo "wrote npt1.in .. npt5.in"; echo "--- npt1.in ---"; cat "$rundir/npt1.in"; echo "--- npt5.in ---"; cat "$rundir/npt5.in"

echo
echo "=== step 10b: write GPU PBS (runs all five stages chained) ==="
pbase="$(basename "$prmtop")"; mbase="$(basename "$minrst")"; rbase="$(basename "$inpref")"
cat > "$rundir/10b_equil.pbs" <<EOF
#!/bin/bash
#PBS -N cm10b_equil
#PBS -l select=1:ncpus=1:ngpus=1:mem=16gb
#PBS -l walltime=01:00:00
#PBS -m ae
#PBS -M $email
#PBS -j oe
#PBS -o $rundir/10b_equil.pbs.out

set -uo pipefail
scratch="/scratch-small-local/\${PBS_JOBID}"
mkdir -p "\$scratch" || { echo "FAIL scratch mkdir"; exit 1; }
echo "host=\$(hostname)  jobid=\${PBS_JOBID:-UNSET}  start=\$(date)"
nvidia-smi -L 2>/dev/null || { echo "FAIL no GPU visible"; exit 1; }

cp "$rundir"/npt1.in "$rundir"/npt2.in "$rundir"/npt3.in "$rundir"/npt4.in "$rundir"/npt5.in \\
   "$prmtop" "$minrst" "$inpref" "\$scratch"/
cd "\$scratch"

set +u
export PERL5LIB="\${PERL5LIB:-}" PYTHONPATH="\${PYTHONPATH:-}"
module load app/amber22/22
set -u
command -v pmemd.cuda >/dev/null || { echo "FAIL pmemd.cuda not found"; exit 1; }
[[ -e "$cuda_lib/libcufft.so.10" ]]   || { echo "FAIL CUDA libs missing at $cuda_lib"; exit 1; }
[[ -e "$blas_lib/libopenblas.so.0" ]] || { echo "FAIL openblas missing at $blas_lib"; exit 1; }
export LD_LIBRARY_PATH="$cuda_lib:$blas_lib:\${LD_LIBRARY_PATH:-}"

run () {  # tag  input_coords  ref(optional)
  local tag="\$1" cin="\$2" ref="\${3:-}"
  local refflag=""; [[ -n "\$ref" ]] && refflag="-ref \$ref"
  echo "=== \$tag (\$(date +%T)) ==="
  pmemd.cuda -O -i \${tag}.in -o \${tag}.out -p "$pbase" -c "\$cin" -r \${tag}.rst7 -inf \${tag}.mdinfo \$refflag
  local rc=\$?
  cp -f \${tag}.out \${tag}.mdinfo "$rundir"/ 2>/dev/null || true
  if [[ \$rc -ne 0 || ! -s \${tag}.rst7 ]]; then
    echo "FAIL stage \$tag rc=\$rc"; tail -30 \${tag}.out 2>/dev/null || true
    cd "$rundir"; rm -rf "\$scratch"; exit 1
  fi
  cp -f \${tag}.rst7 "$rundir"/ 2>/dev/null || true
  grep -E "TIME\\(PS\\) =|Density|ns/day" \${tag}.out | tail -3 || true
}

run npt1 "$mbase" "$rbase"
run npt2 npt1.rst7 "$rbase"
run npt3 npt2.rst7 "$rbase"
run npt4 npt3.rst7 "$rbase"
run npt5 npt4.rst7

cd "$rundir"; rm -rf "\$scratch"
[[ -s "$rundir/npt5.rst7" ]] || { echo "FAIL no npt5.rst7 (final equilibrated restart)"; exit 1; }
echo "end=\$(date)"
echo "FINAL PASS: 10b equilibration complete (final restart: npt5.rst7)"
EOF
cat "$rundir/10b_equil.pbs"

echo
echo "=== step 10b: submit ==="
jobid="$(qsub "$rundir/10b_equil.pbs")"
echo "$jobid" > "$rundir/10b_equil_jobid.txt"
echo "PASS submitted: $jobid"
echo "STEP 10b SUBMITTED"

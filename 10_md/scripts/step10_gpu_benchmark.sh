#!/usr/bin/env bash
# Step 10 GPU benchmark - short restrained NPT on pmemd.cuda from the 10a restart.
# Locks the GPU recipe (incl. the CUDA-11 runtime the amber22 module does not set)
# and measures ns/day to size the 20 ns production into one job.
set -euo pipefail
root="$HOME/system_development"
build="$root/03_amber/tleap_build"
prmtop="$build/complex_solvated.prmtop"
minrst="$root/04_amber_md/10a_min/10a_min.rst7"
inpref="$build/complex_solvated.inpcrd"
rundir="$root/04_amber_md/10b_gpu_benchmark"
admin="$root/00_admin"
email="18660916@sun.ac.za"
# CUDA-11 runtime + OpenBLAS that pmemd.cuda needs (no cuda module exists; validated by ldd)
cuda_lib="/apps/mambaforge/pkgs/cudatoolkit-11.8.0-h37601d7_11/lib"
blas_lib="/apps/mambaforge/envs/medaka/lib"
mkdir -p "$rundir" "$admin"

echo "=== GPU benchmark: input presence ==="
for f in "$prmtop" "$minrst" "$inpref"; do
  [[ -e "$f" ]] || { echo "FAIL missing: $f"; exit 1; }
  echo "PASS $f"
done

echo
echo "=== write benchmark mdin (20 ps NPT, solute restrained 200, heat 0->300) ==="
cat > "$rundir/bench.in" <<'MDIN'
GPU benchmark: 20 ps NPT, solute heavy restrained 200 kcal/mol/A^2, heat 0->300 K
 &cntrl
  imin=0, irest=0, ntx=1,
  ntb=2, ntp=1, barostat=2, pres0=1.0, cut=9.0,
  nstlim=10000, dt=0.002,
  ntc=2, ntf=2,
  tempi=0.0, temp0=300.0, ntt=3, gamma_ln=2.0, ig=531984,
  ntr=1, restraint_wt=200.0, restraintmask='!(:WAT,Na+) & !@H=',
  ntpr=500, ntwx=0, ntwr=10000, ntxo=1, ioutfm=1, iwrap=1,
 /
MDIN
cat "$rundir/bench.in"

echo
echo "=== write GPU PBS ==="
pbase="$(basename "$prmtop")"; mbase="$(basename "$minrst")"; rbase="$(basename "$inpref")"
cat > "$rundir/bench.pbs" <<EOF
#!/bin/bash
#PBS -N cm10_gpubench
#PBS -l select=1:ncpus=1:ngpus=1:mem=16gb
#PBS -l walltime=00:30:00
#PBS -m ae
#PBS -M $email
#PBS -j oe
#PBS -o $rundir/bench.pbs.out

set -uo pipefail
scratch="/scratch-small-local/\${PBS_JOBID}"
mkdir -p "\$scratch" || { echo "FAIL scratch mkdir"; exit 1; }
echo "host=\$(hostname)  jobid=\${PBS_JOBID:-UNSET}  start=\$(date)"

echo "=== GPU visible? ==="
nvidia-smi -L 2>/dev/null || { echo "FAIL no GPU visible (did not land on a GPU)"; exit 1; }

cp "$rundir/bench.in" "$prmtop" "$minrst" "$inpref" "\$scratch"/
cd "\$scratch"

set +u
export PERL5LIB="\${PERL5LIB:-}" PYTHONPATH="\${PYTHONPATH:-}"
module load app/amber22/22
set -u
echo "AMBERHOME=\${AMBERHOME:-UNSET}"
command -v pmemd.cuda >/dev/null || { echo "FAIL pmemd.cuda not found"; exit 1; }

# CUDA-11 runtime for pmemd.cuda (amber22 module does not provide it; no cuda module exists)
[[ -e "$cuda_lib/libcufft.so.10" ]]     || { echo "FAIL CUDA libs missing at $cuda_lib"; exit 1; }
[[ -e "$blas_lib/libopenblas.so.0" ]]   || { echo "FAIL openblas missing at $blas_lib"; exit 1; }
export LD_LIBRARY_PATH="$cuda_lib:$blas_lib:\${LD_LIBRARY_PATH:-}"

echo "=== running benchmark (pmemd.cuda) ==="
pmemd.cuda -O -i bench.in -o bench.out \\
  -p "$pbase" -c "$mbase" -ref "$rbase" \\
  -r bench.rst7 -inf bench.mdinfo
rc=\$?
cp -f bench.out bench.mdinfo "$rundir"/ 2>/dev/null || true
cd "$rundir"; rm -rf "\$scratch"
[[ \$rc -eq 0 ]] || { echo "FAIL pmemd.cuda rc=\$rc"; tail -30 "$rundir/bench.out" 2>/dev/null; exit 1; }

echo "=== throughput ==="
grep -iE "ns/day|seconds/ns|CUDA Device|Device Name|GPU" "$rundir/bench.out" | head -20 || true
echo "end=\$(date)"
echo "FINAL PASS: GPU benchmark complete"
EOF
cat "$rundir/bench.pbs"

echo
echo "=== submit ==="
jobid="$(qsub "$rundir/bench.pbs")"
echo "$jobid" > "$rundir/bench_jobid.txt"
echo "PASS submitted: $jobid"
echo "monitor: qstat -u \$USER"
echo "GPU BENCHMARK SUBMITTED"

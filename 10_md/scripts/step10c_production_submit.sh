#!/usr/bin/env bash
# Step 10c - production: 20 ns unrestrained NPT (300 K, 1 atm), one GPU job.
# From the equilibrated 10b restart (npt5.rst7). dt=0.002 + SHAKE, ntwx=500 ->
# 1 ps/frame -> 20,000 frames (paper).
# Walltime = 168 h (week queue): the GPU nodes are shared and sustained throughput
# swings from ~150 down to ~18 ns/day under load. A week ceiling makes contention
# irrelevant - it completes in ONE job whenever it completes. ncpus=8 feeds the GPU
# on a busy host; scratch self-selected (>=25 GB); trap cleans scratch on any exit.
set -euo pipefail
root="$HOME/system_development"
build="$root/03_amber/tleap_build"
prmtop="$build/complex_solvated.prmtop"
startrst="$root/04_amber_md/10b_equil/npt5.rst7"
rundir="$root/04_amber_md/10c_production"
admin="$root/00_admin"
email="18660916@sun.ac.za"
cuda_lib="/apps/mambaforge/pkgs/cudatoolkit-11.8.0-h37601d7_11/lib"
blas_lib="/apps/mambaforge/envs/medaka/lib"
seed=531984
mkdir -p "$rundir" "$admin"

echo "=== step 10c: input presence ==="
for f in "$prmtop" "$startrst"; do
  [[ -e "$f" ]] || { echo "FAIL missing: $f"; exit 1; }
  echo "PASS $f"
done

echo
echo "=== step 10c: write production mdin (20 ns, unrestrained NPT) ==="
cat > "$rundir/prod.in" <<MDIN
Production: 20 ns unrestrained NPT, 300 K, 1 atm, 1 ps/frame (20000 frames)
 &cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, pres0=1.0, cut=9.0,
  nstlim=10000000, dt=0.002,
  ntc=2, ntf=2,
  temp0=300.0, ntt=3, gamma_ln=2.0, ig=$seed,
  ntpr=5000, ntwx=500, ntwr=500000, ntxo=1, ioutfm=1, iwrap=1,
 /
MDIN
cat "$rundir/prod.in"

echo
echo "=== step 10c: write GPU PBS ==="
pbase="$(basename "$prmtop")"; sbase="$(basename "$startrst")"
cat > "$rundir/10c_prod.pbs" <<EOF
#!/bin/bash
#PBS -N cm10c_prod
#PBS -l select=1:ncpus=8:ngpus=1:mem=16gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M $email
#PBS -j oe
#PBS -o $rundir/10c_prod.pbs.out

set -uo pipefail
echo "host=\$(hostname)  jobid=\${PBS_JOBID:-UNSET}  start=\$(date)"
nvidia-smi -L 2>/dev/null || { echo "FAIL no GPU visible"; exit 1; }

# choose a scratch dir with >= 25 GB free (trajectory is ~13 GB); fail fast if none
need_kb=26214400
scratch=""
for base in /scratch-small-local /scratch-large-network; do
  d="\$base/\${PBS_JOBID}"
  mkdir -p "\$d" 2>/dev/null || continue
  avail=\$(df -Pk "\$d" 2>/dev/null | awk 'END{print \$4}')
  if [[ "\${avail:-0}" -ge "\$need_kb" ]]; then scratch="\$d"; echo "scratch=\$scratch (avail \${avail} KB)"; break; fi
  rmdir "\$d" 2>/dev/null || true
done
[[ -n "\$scratch" ]] || { echo "FAIL no scratch dir with >=25GB free for the trajectory"; exit 1; }
# clean our own scratch on ANY exit (normal end, qdel, or walltime kill)
trap 'cd "$rundir" 2>/dev/null; rm -rf "\$scratch" 2>/dev/null' EXIT TERM

cp "$rundir/prod.in" "$prmtop" "$startrst" "\$scratch"/
cd "\$scratch"

set +u
export PERL5LIB="\${PERL5LIB:-}" PYTHONPATH="\${PYTHONPATH:-}"
module load app/amber22/22
set -u
command -v pmemd.cuda >/dev/null || { echo "FAIL pmemd.cuda not found"; exit 1; }
[[ -e "$cuda_lib/libcufft.so.10" ]]   || { echo "FAIL CUDA libs missing at $cuda_lib"; exit 1; }
[[ -e "$blas_lib/libopenblas.so.0" ]] || { echo "FAIL openblas missing at $blas_lib"; exit 1; }
export LD_LIBRARY_PATH="$cuda_lib:$blas_lib:\${LD_LIBRARY_PATH:-}"

echo "=== running 20 ns production (pmemd.cuda) ==="
pmemd.cuda -O -i prod.in -o prod.out -p "$pbase" -c "$sbase" \\
  -r prod.rst7 -x prod.nc -inf prod.mdinfo
rc=\$?

echo "=== copying results to \$HOME (1 TB tier) ==="
cp -f prod.out prod.mdinfo "$rundir"/ 2>/dev/null || true
if [[ \$rc -eq 0 ]]; then
  cp -f prod.rst7 prod.nc "$rundir"/ 2>/dev/null || true
fi

[[ \$rc -eq 0 ]] || { echo "FAIL pmemd.cuda rc=\$rc"; tail -30 "$rundir/prod.out" 2>/dev/null; exit 1; }
[[ -s "$rundir/prod.nc" && -s "$rundir/prod.rst7" ]] || { echo "FAIL missing trajectory or restart"; exit 1; }
echo "=== summary ==="
grep -E "ns/day|TIME\\(PS\\) =|Density" "$rundir/prod.out" | tail -5 || true
ls -lh "$rundir/prod.nc" "$rundir/prod.rst7"
echo "end=\$(date)"
echo "FINAL PASS: 10c production complete (trajectory prod.nc, restart prod.rst7)"
EOF
cat "$rundir/10c_prod.pbs"

echo
echo "=== step 10c: submit ==="
jobid="$(qsub "$rundir/10c_prod.pbs")"
echo "$jobid" > "$rundir/10c_prod_jobid.txt"
echo "PASS submitted: $jobid"
echo "STEP 10c SUBMITTED"

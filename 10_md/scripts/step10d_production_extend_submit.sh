#!/usr/bin/env bash
# Step 10d - production extension: 40 ns unrestrained NPT (300 K, 1 atm), one GPU
# job, continuing from the 10c restart (prod.rst7).
#
# This is step10c_production_submit.sh with four changes and nothing else:
#   nstlim   10000000 -> 20000000   (20 ns -> 40 ns)
#   startrst npt5.rst7 -> prod.rst7 (continue the production run, not the equil)
#   rundir   10c_production -> 10d_production_extend
#   scratch  25 GB -> 55 GB         (the trajectory roughly doubles)
# Every other setting, including the random seed, the module handling and the
# scratch discipline, is carried over verbatim so the two runs are one trajectory.
#
# WHY EXTEND
# The backbone RMSD of the 20 ns run rises monotonically over the first 10 ns and
# is stable thereafter, so only the second half is equilibrated. Within that half
# the backbone autocorrelation decays with a time constant of about 800-860 ps,
# so the 9 ns of equilibrated trajectory holds only ~5.6 independent protein
# configurations however many frames are drawn from it. Thirty independent
# configurations need ~48 ns equilibrated; 9 exist, so 40 more are run.
#
# dt=0.002 + SHAKE, ntwx=500 -> 1 ps/frame -> 40,000 frames.
set -euo pipefail
root="$HOME/system_development"
build="$root/03_amber/tleap_build"
prmtop="$build/complex_solvated.prmtop"
startrst="$root/04_amber_md/10c_production/prod.rst7"
rundir="$root/04_amber_md/10d_production_extend"
admin="$root/00_admin"
email="18660916@sun.ac.za"
cuda_lib="/apps/mambaforge/pkgs/cudatoolkit-11.8.0-h37601d7_11/lib"
blas_lib="/apps/mambaforge/envs/medaka/lib"
seed=531984
mkdir -p "$rundir" "$admin"

echo "=== step 10d: input presence ==="
for f in "$prmtop" "$startrst"; do
  [[ -e "$f" ]] || { echo "FAIL missing: $f"; exit 1; }
  echo "PASS $f"
done

echo
echo "=== step 10d: write extension mdin (40 ns, unrestrained NPT) ==="
cat > "$rundir/prod_ext.in" <<MDIN
Production extension: 40 ns unrestrained NPT, 300 K, 1 atm, 1 ps/frame (40000 frames)
 &cntrl
  imin=0, irest=1, ntx=5,
  ntb=2, ntp=1, barostat=2, pres0=1.0, cut=9.0,
  nstlim=20000000, dt=0.002,
  ntc=2, ntf=2,
  temp0=300.0, ntt=3, gamma_ln=2.0, ig=$seed,
  ntpr=5000, ntwx=500, ntwr=500000, ntxo=1, ioutfm=1, iwrap=1,
 /
MDIN
cat "$rundir/prod_ext.in"

echo
echo "=== step 10d: write GPU PBS ==="
pbase="$(basename "$prmtop")"; sbase="$(basename "$startrst")"
cat > "$rundir/10d_prod_ext.pbs" <<EOF
#!/bin/bash
#PBS -N cm10d_prodext
#PBS -l select=1:ncpus=8:ngpus=1:mem=16gb
#PBS -l walltime=168:00:00
#PBS -m ae
#PBS -M $email
#PBS -j oe
#PBS -o $rundir/10d_prod_ext.pbs.out

set -uo pipefail
echo "host=\$(hostname)  jobid=\${PBS_JOBID:-UNSET}  start=\$(date)"
nvidia-smi -L 2>/dev/null || { echo "FAIL no GPU visible"; exit 1; }

# choose a scratch dir with >= 55 GB free (the 40 ns trajectory is ~26 GB and the
# restart and output add to it); fail fast if none
need_kb=57671680
scratch=""
for base in /scratch-small-local /scratch-large-network; do
  d="\$base/\${PBS_JOBID}"
  mkdir -p "\$d" 2>/dev/null || continue
  avail=\$(df -Pk "\$d" 2>/dev/null | awk 'END{print \$4}')
  if [[ "\${avail:-0}" -ge "\$need_kb" ]]; then scratch="\$d"; echo "scratch=\$scratch (avail \${avail} KB)"; break; fi
  rmdir "\$d" 2>/dev/null || true
done
[[ -n "\$scratch" ]] || { echo "FAIL no scratch dir with >=55GB free for the trajectory"; exit 1; }
# clean our own scratch on ANY exit (normal end, qdel, or walltime kill)
trap 'cd "$rundir" 2>/dev/null; rm -rf "\$scratch" 2>/dev/null' EXIT TERM

cp "$rundir/prod_ext.in" "$prmtop" "$startrst" "\$scratch"/
cd "\$scratch"

set +u
export PERL5LIB="\${PERL5LIB:-}" PYTHONPATH="\${PYTHONPATH:-}"
module load app/amber22/22
set -u
command -v pmemd.cuda >/dev/null || { echo "FAIL pmemd.cuda not found"; exit 1; }
[[ -e "$cuda_lib/libcufft.so.10" ]]   || { echo "FAIL CUDA libs missing at $cuda_lib"; exit 1; }
[[ -e "$blas_lib/libopenblas.so.0" ]] || { echo "FAIL openblas missing at $blas_lib"; exit 1; }
export LD_LIBRARY_PATH="$cuda_lib:$blas_lib:\${LD_LIBRARY_PATH:-}"

echo "=== running 40 ns production extension (pmemd.cuda) ==="
pmemd.cuda -O -i prod_ext.in -o prod_ext.out -p "$pbase" -c "$sbase" \\
  -r prod_ext.rst7 -x prod_ext.nc -inf prod_ext.mdinfo
rc=\$?

echo "=== copying results to \$HOME (1 TB tier) ==="
cp -f prod_ext.out prod_ext.mdinfo "$rundir"/ 2>/dev/null || true
if [[ \$rc -eq 0 ]]; then
  cp -f prod_ext.rst7 prod_ext.nc "$rundir"/ 2>/dev/null || true
fi

[[ \$rc -eq 0 ]] || { echo "FAIL pmemd.cuda rc=\$rc"; tail -30 "$rundir/prod_ext.out" 2>/dev/null; exit 1; }
[[ -s "$rundir/prod_ext.nc" && -s "$rundir/prod_ext.rst7" ]] || { echo "FAIL missing trajectory or restart"; exit 1; }
echo "=== summary ==="
grep -E "ns/day|TIME\\(PS\\) =|Density" "$rundir/prod_ext.out" | tail -5 || true
ls -lh "$rundir/prod_ext.nc" "$rundir/prod_ext.rst7"
echo "end=\$(date)"
echo "FINAL PASS: 10d production extension complete (trajectory prod_ext.nc, restart prod_ext.rst7)"
echo
echo "The extension continues 10c, so the two are one trajectory: prod.nc then"
echo "prod_ext.nc, 60 ns total, of which the first 10 ns are pre-equilibration."
EOF
cat "$rundir/10d_prod_ext.pbs"

echo
echo "=== step 10d: submit ==="
jobid="$(qsub "$rundir/10d_prod_ext.pbs")"
echo "$jobid" > "$rundir/10d_prod_ext_jobid.txt"
echo "PASS submitted: $jobid"
echo "STEP 10d SUBMITTED"

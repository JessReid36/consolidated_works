#!/usr/bin/env bash
# Step 11 plots - render RMSD-vs-time and per-residue RMSF from the .dat files.
# Ensures matplotlib via whatever pip is available (pip3 / python3 -m pip / module).
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
d="$HOME/system_development/04_amber_md/11_analysis"

ensure_mpl() {
  python3 -c "import matplotlib" 2>/dev/null && { echo "matplotlib present"; return 0; }
  echo "matplotlib missing - trying to install (user)..."
  if   python3 -m pip --version >/dev/null 2>&1; then python3 -m pip install --user matplotlib 2>&1 | tail -3
  elif command -v pip3 >/dev/null 2>&1;          then pip3 install --user matplotlib 2>&1 | tail -3
  elif command -v pip  >/dev/null 2>&1;          then pip  install --user matplotlib 2>&1 | tail -3
  else echo "no pip found on this node"; fi
  python3 -c "import matplotlib" 2>/dev/null
}

if ! ensure_mpl; then
  echo "COULD NOT GET matplotlib on the HPC - the .dat files are complete; plot them locally instead."
  echo "  scp 18660916@hpc1.sun.ac.za:'~/system_development/04_amber_md/11_analysis/*.dat' ."
  exit 0
fi

python3 - "$d" <<'PY'
import sys, os
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=sys.argv[1]
p=os.path.join(d,"rmsd_vs_time.dat")
if os.path.exists(p):
    a=np.loadtxt(p); ns=a[:,1]/1000.0; r=a[:,2]
    plt.figure(figsize=(8,4)); plt.plot(ns,r,lw=0.6,color="#1f77b4")
    plt.xlabel("production time (ns)"); plt.ylabel("backbone RMSD (Å)")
    plt.title("BsCM trimer backbone RMSD vs frame 0"); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig(os.path.join(d,"rmsd_vs_time.png"),dpi=150); plt.close(); print("wrote rmsd_vs_time.png")
p=os.path.join(d,"rmsf_per_residue.dat")
if os.path.exists(p):
    chains=np.genfromtxt(p,dtype=str,usecols=0); resid=np.genfromtxt(p,usecols=1); rmsf=np.genfromtxt(p,usecols=3)
    plt.figure(figsize=(10,4))
    for ch,col in zip(["A","B","C"],["#1f77b4","#ff7f0e","#2ca02c"]):
        m=chains==ch
        if m.any(): plt.plot(resid[m],rmsf[m],lw=0.8,color=col,label="chain "+ch)
    plt.xlabel("residue (global index)"); plt.ylabel("Cα RMSF (Å)")
    plt.title("Per-residue Cα RMSF"); plt.legend(); plt.grid(alpha=.3); plt.tight_layout()
    plt.savefig(os.path.join(d,"rmsf_per_residue.png"),dpi=150); plt.close(); print("wrote rmsf_per_residue.png")
PY
echo "PLOTS DONE (in $d)"

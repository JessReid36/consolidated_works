#!/usr/bin/env bash
set -uo pipefail
ORCA=/home/apps2/ORCA/6.0.1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "$PWD"
ACT="$(cat active.txt)"
SPRING=400.0
prev="win_00.pdb"
n=0
: > scan_progress.tsv
: > scan_energies.tsv
while read -r BRK FRM <&3; do
  n=$((n+1)); w=$(printf "%02d" $n); base="scan_${w}"
  cat > "${base}.inp" <<INP
! QMMM B3LYP D3BJ def2-SVP def2/J RIJCOSX MD
%maxcore 3000
%pal nprocs 8 end
%scf MaxIter 200 end
%md
  Manage_Colvar Define 1 Distance Atom 6215 Atom 6214
  Manage_Colvar Define 2 Distance Atom 6219 Atom 6207
  Restraint Add Colvar 1 Harmonic Target ${BRK}_A Spring ${SPRING}
  Restraint Add Colvar 2 Harmonic Target ${FRM}_A Spring ${SPRING}
  Minimize LBFGS MaxGrad 5.0 RMSGrad 1.0 History 20
end
%qmmm
  QMAtoms {6207:6230} end
  ActiveAtoms {$ACT} end
  ORCAFFFilename "complex_solvated.ORCAFF.prms"
end
*pdbfile -2 1 ${prev}
INP
  echo "=== window $w  target break=$BRK form=$FRM  from $prev  $(date +%H:%M) ==="
  "$ORCA/orca" "${base}.inp" > "${base}.out" 2>&1 </dev/null

  traj="${base}-postrj-all.xyz"
  if [[ ! -s "$traj" ]]; then echo "  NO TRAJ for $w - stopping"; break; fi

  # build win_NN.pdb: reactant PDB records + last-frame XYZ coords; also record distances + energy
  python3 - "$traj" win_00.pdb "win_${w}.pdb" "$w" <<'PY'
import sys, numpy as np
traj, tmpl, outpdb, w = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
# read last frame of the multi-frame xyz
lines=open(traj).read().splitlines()
nat=int(lines[0].split()[0]); blk=nat+2
nfr=len(lines)//blk
s=(nfr-1)*blk
# energy from that frame's comment line
import re
m=re.search(r'E_Pot=(-?\d+\.\d+)', lines[s+1]); epot=float(m.group(1)) if m else float('nan')
coords=[]
for ln in lines[s+2:s+2+nat]:
    p=ln.split(); coords.append((float(p[1]),float(p[2]),float(p[3])))
# splice into template PDB (preserve columns, swap x/y/z)
out=[]; ai=0
for ln in open(tmpl):
    if ln[:6] in ("ATOM  ","HETATM"):
        x,y,z=coords[ai]; ai+=1
        out.append(f"{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:]}")
    else:
        out.append(ln.rstrip("\n"))
open(outpdb,"w").write("\n".join(out)+"\n")
xyz=lambda i: np.array(coords[i])
brk=np.linalg.norm(xyz(6215)-xyz(6214)); frm=np.linalg.norm(xyz(6219)-xyz(6207))
open("scan_progress.tsv","a").write(f"{w}\t{brk:.3f}\t{frm:.3f}\t{brk-frm:+.3f}\n")
open("scan_energies.tsv","a").write(f"{w}\t{brk-frm:+.3f}\t{epot:.6f}\n")
print(f"  achieved break={brk:.3f} form={frm:.3f} r={brk-frm:+.3f}  E={epot:.5f}")
PY

  prev="win_${w}.pdb"
  rm -f "${base}-postrj-all.xyz" "${base}.mdrestart" "${base}.densities" "${base}.gbw" 2>/dev/null
done 3< targets.txt
echo "SCAN COMPLETE through window $w"

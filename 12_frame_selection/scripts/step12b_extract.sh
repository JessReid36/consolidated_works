#!/usr/bin/env bash
# Step 12b - extract each selected frame (from 12a manifest) to QM/MM-ready inputs:
# a full-system AMBER rst7 (matches complex_solvated.prmtop) + the QM atom list
# (1-based indices of the chosen CHA residue = the substrate-only QM region).
# Folds in the old protocol's separate frame-extraction (29c) and QM-region (30) steps.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
traj="$root/04_amber_md/10c_production/prod.nc"
# the production run was extended; the two files are one trajectory
traj2="$root/04_amber_md/10d_production_extend/prod_ext.nc"
seldir="$root/05_qmmm/12_frame_selection"
manifest="$seldir/selection_manifest.tsv"
outdir="$seldir/frames"; mkdir -p "$outdir"
for f in "$prmtop" "$traj" "$manifest"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done
[[ -s "$traj2" ]] || { echo "note: no extension trajectory at $traj2; "\
  "only frames of the first run can be extracted"; traj2=""; }

python3 - "$prmtop" "$traj" "$manifest" "$outdir" "$traj2" <<'PY'
import sys, os, re, struct
import numpy as np
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
prmtop, traj, manifest, outdir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
traj2 = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] else None
def parse_prmtop(path):
    flags,fmts,cur={},{},None
    for ln in open(path,errors="replace"):
        ln=ln.rstrip("\n")
        if ln.startswith("%FLAG"): cur=ln.split()[1]; flags[cur]=[]
        elif ln.startswith("%FORMAT"): fmts[cur]=ln[ln.find("(")+1:ln.rfind(")")]
        elif ln.startswith("%"): continue
        elif cur is not None: flags[cur].append(ln)
    def width(f):
        m=re.match(r"\s*\d*[aAiIeEfFgG](\d+)",f); return int(m.group(1)) if m else None
    def toks(name,cast):
        w=width(fmts.get(name,"")); out=[]
        for ln in flags.get(name,[]):
            if w: out+=[ln[i:i+w].strip() for i in range(0,len(ln.rstrip()),w) if ln[i:i+w].strip()]
            else: out+=ln.split()
        return [cast(x) for x in out]
    return toks("POINTERS",int)[0],toks("ATOM_NAME",str),toks("RESIDUE_LABEL",str),toks("RESIDUE_POINTER",int)
natom,names,labels,resptr=parse_prmtop(prmtop)
nres=len(labels); starts=[p-1 for p in resptr]+[natom]
def nc_open(path):
    data=open(path,"rb").read(1<<20); fsz=os.path.getsize(path); ver=data[3]; off64=(ver==2); pos=[4]
    def u32():
        v=struct.unpack_from(">I",data,pos[0])[0]; pos[0]+=4; return v
    def off():
        v=struct.unpack_from(">Q" if off64 else ">I",data,pos[0])[0]; pos[0]+=8 if off64 else 4; return v
    def nm():
        n=u32(); s=data[pos[0]:pos[0]+n].decode("ascii","replace"); pos[0]+=n+((4-(n%4))%4); return s
    numrecs=u32(); dims=[]; tag=u32(); ne=u32()
    if tag==0x0A:
        for _ in range(ne): dn=nm(); dl=u32(); dims.append((dn,dl))
    def skip_att():
        t=u32(); n=u32()
        if t==0x0C:
            for _ in range(n):
                nm(); tp=u32(); k=u32(); nb=k*{1:1,2:1,3:2,4:4,5:4,6:8}.get(tp,1); pos[0]+=nb+((4-(nb%4))%4)
    skip_att(); rec=[]; tag=u32(); nv=u32()
    if tag==0x0B:
        for _ in range(nv):
            vn=nm(); nd=u32(); dimids=[u32() for _ in range(nd)]; skip_att(); tp=u32(); vs=u32(); bg=off()
            if nd>0 and dims[dimids[0]][1]==0: rec.append({"name":vn,"vsize":vs,"begin":bg})
    dimd={n:l for n,l in dims}; recsize=sum(v["vsize"] for v in rec)
    frames=numrecs if numrecs!=0xFFFFFFFF else (fsz-min(v["begin"] for v in rec))//recsize
    return {"frames":frames,"recsize":recsize,"rec":rec,"atom":dimd.get("atom")}
def open_traj(path):
    """Header and the three record variables of one trajectory file."""
    H=nc_open(path); assert H["atom"]==natom, "atom mismatch in %s"%path
    return {"path":path, "H":H,
            "c":next(v for v in H["rec"] if v["name"]=="coordinates"),
            "l":next((v for v in H["rec"] if v["name"]=="cell_lengths"),None),
            "a":next((v for v in H["rec"] if v["name"]=="cell_angles"),None)}

T1=open_traj(traj)
T2=open_traj(traj2) if traj2 else None
N1=T1["H"]["frames"]
print("trajectory 1: %s, %d frames"%(os.path.basename(traj),N1))
if T2:
    print("trajectory 2: %s, %d frames"%(os.path.basename(traj2),T2["H"]["frames"]))
    print("frames 0 to %d are read from the first, %d and above from the second"
          %(N1-1,N1))

def frame(i):
    """Read production frame i, spanning the two files.

    The two trajectories are one run: the extension restarted from the first
    file's restart with velocities. Production frame i therefore lives in the
    first file if i is below its length and in the second at i-N1 otherwise.
    The boundary comes from the first file's header, not from an assumed
    length."""
    if i < N1 or T2 is None:
        t, j = T1, i
        if i >= N1:
            sys.exit("FAIL frame %d is beyond %s (%d frames) and no extension "
                     "trajectory was given"%(i,os.path.basename(traj),N1))
    else:
        t, j = T2, i - N1
        if j >= t["H"]["frames"]:
            sys.exit("FAIL frame %d is beyond the end of both trajectories"%i)
    path, H, cvar, lvar, avar = t["path"], t["H"], t["c"], t["l"], t["a"]
    o=cvar["begin"]+j*H["recsize"]
    xyz=np.array(np.memmap(path,dtype=">f4",mode="r",offset=o,shape=(natom,3)),float)
    box=[0.0,0.0,0.0,90.0,90.0,90.0]
    if lvar is not None:
        bo=lvar["begin"]+j*H["recsize"]; box[:3]=list(np.array(np.memmap(path,dtype=">f8",mode="r",offset=bo,shape=(3,)),float))
    if avar is not None:
        ao=avar["begin"]+j*H["recsize"]; box[3:]=list(np.array(np.memmap(path,dtype=">f8",mode="r",offset=ao,shape=(3,)),float))
    return xyz,box
def write_rst7(path, xyz, box, title):
    with open(path,"w") as f:
        f.write(title[:80]+"\n"); f.write("%6d\n"%len(xyz))
        flat=xyz.reshape(-1)
        for k in range(0,len(flat),6):
            f.write("".join("%12.7f"%v for v in flat[k:k+6])+"\n")
        f.write("".join("%12.7f"%v for v in box)+"\n")
rows=[ln.rstrip("\n").split("\t") for ln in open(manifest) if ln.strip() and not ln.startswith("#")]
out_manifest=os.path.join(outdir,"frames_manifest.tsv")
mf=open(out_manifest,"w"); mf.write("# idx\tframe\tsite\tCHA_resid\trst7\tqm_atoms_file\tn_qm\n")
print("STEP 12b - extract %d selected frames"%len(rows))
for r in rows:
    idx=int(r[0]); fr=int(r[1]); site=r[3]; resid=int(r[4])
    xyz,box=frame(fr)
    if not np.isfinite(xyz).all(): sys.exit("FAIL non-finite coords in frame %d"%fr)
    tag="frame_%05d_%s"%(fr,site.replace("#",""))
    rst=os.path.join(outdir,tag+".rst7"); write_rst7(rst, xyz, box, "QMMM reactant %s frame %d"%(site,fr))
    ri=resid-1; qm=list(range(starts[ri]+1, starts[ri+1]+1))
    qmf=os.path.join(outdir,tag+".qmatoms"); open(qmf,"w").write(" ".join(map(str,qm))+"\n")
    mf.write("%d\t%d\t%s\t%d\t%s\t%s\t%d\n"%(idx,fr,site,resid,os.path.basename(rst),os.path.basename(qmf),len(qm)))
    print("  #%2d  %s  rst7=%d atoms  QM=%d atoms (res %d: %d..%d)"%(idx,tag,len(xyz),len(qm),resid,qm[0],qm[-1]))
mf.close()
print("wrote", out_manifest)
PY
echo "STEP 12b DONE"

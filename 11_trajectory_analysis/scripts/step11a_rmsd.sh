#!/usr/bin/env bash
# Step 11a - protein backbone RMSD of the 20 ns production vs frame 0 (equilibration/
# stability QC; Agbaglo-style backbone C,O,CA,N,H). cpptraj-free: stdlib NetCDF header
# reader + numpy, with molecule-wise imaging (undoes iwrap chain splits) + Kabsch align.
# Single-threaded to avoid the login-node OpenBLAS RLIMIT_NPROC segfault.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
traj="$root/04_amber_md/10c_production/prod.nc"
outdir="$root/04_amber_md/11_analysis"; mkdir -p "$outdir"
outdat="$outdir/rmsd_vs_time.dat"

for f in "$prmtop" "$traj"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done

python3 - "$prmtop" "$traj" "$outdat" <<'PY'
import sys, os, re, struct
import numpy as np
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
prmtop, traj, outdat = sys.argv[1], sys.argv[2], sys.argv[3]
PROT_RES={"ALA","ARG","ASN","ASP","ASH","CYS","CYX","CYM","GLN","GLU","GLH","GLY","HID","HIE","HIP",
          "ILE","LEU","LYS","LYN","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
BB={"N","CA","C","O","H"}
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
    return (toks("POINTERS",int)[0],toks("ATOM_NAME",str),toks("RESIDUE_LABEL",str),
            toks("RESIDUE_POINTER",int),toks("ATOMS_PER_MOLECULE",int))
natom,names,labels,resptr,apm=parse_prmtop(prmtop)
nres=len(labels); starts=[p-1 for p in resptr]+[natom]
res_of=np.empty(natom,int)
for ri in range(nres): res_of[starts[ri]:starts[ri+1]]=ri
bb_idx=np.array([a for a in range(natom) if labels[res_of[a]] in PROT_RES and names[a] in BB],int)
mol_of=np.empty(natom,int); a0=0
for mi,n in enumerate(apm): mol_of[a0:a0+n]=mi; a0+=n
prot_mols=sorted(set(mol_of[bb_idx].tolist())); ref_mol=prot_mols[0]
mol_ranges={}; a0=0
for mi,n in enumerate(apm): mol_ranges[mi]=(a0,a0+n); a0+=n
n_solute=max(mol_ranges[m][1] for m in prot_mols)
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
H=nc_open(traj); assert H["atom"]==natom,"atom mismatch %d vs %d"%(H["atom"],natom)
cvar=next(v for v in H["rec"] if v["name"]=="coordinates")
lvar=next((v for v in H["rec"] if v["name"]=="cell_lengths"),None)
def frame(i):
    off=cvar["begin"]+i*H["recsize"]
    xyz=np.array(np.memmap(traj,dtype=">f4",mode="r",offset=off,shape=(natom,3))[:n_solute],float)
    box=None
    if lvar is not None:
        bo=lvar["begin"]+i*H["recsize"]; box=np.array(np.memmap(traj,dtype=">f8",mode="r",offset=bo,shape=(3,)),float)
    return xyz,box
def image(xyz,box):
    if box is None: return xyz
    x=xyz.copy(); rc=x[mol_ranges[ref_mol][0]:mol_ranges[ref_mol][1]].mean(0)
    for mi in prot_mols:
        s,e=mol_ranges[mi]; x[s:e]-=np.round((x[s:e].mean(0)-rc)/box)*box
    return x
def kabsch(P,Q):
    Pc=P-P.mean(0); Qc=Q-Q.mean(0); M=Pc.T@Qc; U,S,Vt=np.linalg.svd(M)
    d=np.sign(np.linalg.det(Vt.T@U.T)); R=Vt.T@np.diag([1,1,d])@U.T
    return float(np.sqrt(np.mean(np.sum(((Pc@R.T)-Qc)**2,1))))
x0,b0=frame(0); ref=image(x0,b0)[bb_idx]; rms=np.empty(H["frames"])
for i in range(H["frames"]):
    xi,bi=frame(i); rms[i]=kabsch(image(xi,bi)[bb_idx],ref)
with open(outdat,"w") as f:
    f.write("# frame\tprod_ps\trmsd_A  (backbone N,CA,C,O,H vs frame 0, imaged+aligned)\n")
    for i,v in enumerate(rms): f.write("%d\t%d\t%.4f\n"%(i,i+1,v))
w=max(1,H["frames"]//10); first=float(rms[:w].mean()); last=float(rms[-w:].mean())
print("STEP 11a - protein backbone RMSD vs production frame 0")
print("  backbone atoms %d   protein molecules %s   frames %d"%(len(bb_idx),prot_mols,H["frames"]))
print("  mean %.3f A   std %.3f   min %.3f   max %.3f"%(rms.mean(),rms.std(),rms.min(),rms.max()))
print("  first 10%% mean %.3f A   last 10%% mean %.3f A   drift %.3f A"%(first,last,last-first))
print("  wrote",outdat)
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    ps=np.arange(1,H["frames"]+1)
    plt.figure(figsize=(8,4)); plt.plot(ps/1000.0,rms,lw=0.6)
    plt.xlabel("production time (ns)"); plt.ylabel("backbone RMSD (A)")
    plt.title("BsCM trimer backbone RMSD vs frame 0"); plt.tight_layout()
    png=outdat.replace(".dat",".png"); plt.savefig(png,dpi=150); print("  wrote",png)
except Exception as e:
    print("  (no plot: %s - install matplotlib or plot the .dat locally)"%str(e)[:50])
print("  RESULT:", "PASS_STABLE" if abs(last-first)<1.0 and rms.max()<6.0 else "REVIEW_DRIFT_OR_HIGH_RMSD")
PY
echo "STEP 11a DONE"

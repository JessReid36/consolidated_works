#!/usr/bin/env bash
# Step 11d - catalytic-contact validation: does each chorismate make the defining
# Claeyssens contacts? ether O3(=O13)->Arg (Arg90), carboxylates->Arg (Arg7/Arg63),
# hydroxyl O4->acidic (Glu78). Validates each site as a genuine Michaelis complex and
# diagnoses the site-to-site asymmetry. Min-image distances; numpy-only.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
traj="$root/04_amber_md/10c_production/prod.nc"
outdir="$root/04_amber_md/11_analysis"
for f in "$prmtop" "$traj"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done

python3 - "$prmtop" "$traj" "$outdir" <<'PY'
import sys, os, re, struct
from collections import Counter
import numpy as np
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
prmtop, traj, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
ARG_N={"NE","NH1","NH2"}; ACID_O={"OE1","OE2","OD1","OD2"}
CHA_ETHER="O3"; CHA_CARBOX=["O1","O2","O5","O6"]; CHA_OH="O4"
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
res_of=np.empty(natom,int)
for ri in range(nres): res_of[starts[ri]:starts[ri+1]]=ri
argN=np.array([a for a in range(natom) if labels[res_of[a]]=="ARG" and names[a] in ARG_N])
acidO=np.array([a for a in range(natom) if labels[res_of[a]] in ("GLU","ASP","GLH","ASH") and names[a] in ACID_O])
cha=[]
for ri in range(nres):
    if labels[ri]=="CHA":
        amap={names[a]:a for a in range(starts[ri],starts[ri+1])}; cha.append((ri,amap))
if not cha: sys.exit("FAIL no CHA")
chatoms=[a for _,amap in cha for a in amap.values()]
maxidx=int(max(argN.max(), acidO.max(), max(chatoms)))+1     # only read protein+substrate
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
H=nc_open(traj); assert H["atom"]==natom
cvar=next(v for v in H["rec"] if v["name"]=="coordinates")
lvar=next((v for v in H["rec"] if v["name"]=="cell_lengths"),None)
def frame(i):
    off=cvar["begin"]+i*H["recsize"]
    xyz=np.array(np.memmap(traj,dtype=">f4",mode="r",offset=off,shape=(natom,3))[:maxidx],float)
    box=np.array([1e9,1e9,1e9])
    if lvar is not None:
        bo=lvar["begin"]+i*H["recsize"]; box=np.array(np.memmap(traj,dtype=">f8",mode="r",offset=bo,shape=(3,)),float)
    return xyz,box
def mindist(a, pts, box):
    d=pts-a; d-=np.round(d/box)*box; r=np.sqrt((d*d).sum(1)); k=int(r.argmin()); return float(r[k]),k
def rid(ai): ri=res_of[ai]; return "%s%d"%(labels[ri],ri+1)
report={}
for j,(ri,amap) in enumerate(cha):
    acc={"ether":[], "carbox":[], "oh":[]}; who={"ether":Counter(),"carbox":Counter(),"oh":Counter()}
    for i in range(H["frames"]):
        x,box=frame(i)
        d,k=mindist(x[amap[CHA_ETHER]], x[argN], box); acc["ether"].append(d); who["ether"][rid(argN[k])]+=1
        best=(1e9,None)
        for co in CHA_CARBOX:
            if co in amap:
                d,k=mindist(x[amap[co]], x[argN], box)
                if d<best[0]: best=(d, argN[k])
        acc["carbox"].append(best[0]); who["carbox"][rid(best[1])]+=1
        if len(acidO):
            d,k=mindist(x[amap[CHA_OH]], x[acidO], box); acc["oh"].append(d); who["oh"][rid(acidO[k])]+=1
    report[j]=(ri,acc,who)
print("STEP 11d - catalytic-contact validation (Claeyssens key contacts)")
print("  Arg guanidinium N: %d   acidic O: %d   CHA copies: %d   frames: %d"%(len(argN),len(acidO),len(cha),H["frames"]))
for j,(ri,acc,who) in report.items():
    print("  CHA#%d (res %d):"%(j+1,ri+1))
    for key,lbl in [("ether","ether O3->Arg  (Arg90/O13, the key TS contact)"),
                    ("carbox","carboxylate->Arg (Arg7/Arg63)"),
                    ("oh","hydroxyl O4->acidic (Glu78)")]:
        a=np.array(acc[key]); part=who[key].most_common(1)[0] if who[key] else ("none",0)
        frac=100.0*part[1]/len(a) if len(a) else 0
        print("    %-44s mean %.2f A  min %.2f | %s (%.0f%% of frames)"%(lbl,a.mean(),a.min(),part[0],frac))
PY
echo "STEP 11d DONE"

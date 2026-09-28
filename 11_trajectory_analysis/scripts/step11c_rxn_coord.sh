#!/usr/bin/env bash
# Step 11c - per-chorismate reaction-coordinate tracking over the 20 ns production.
# Claeyssens r = d(C2-O13) - d(C4-C14), mapped to cha_gaff names:
#   breaking C2-O13 -> C4-O3 ;  forming C4-C14 -> C6-C1 ;  carboxylate C16-C17 -> C3-C10.
# Confirms all three substrates stay bound AND produces the near-attack metric that
# step 12 (frame selection for QM/MM) keys on. numpy-only, intramolecular distances.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
traj="$root/04_amber_md/10c_production/prod.nc"
outdir="$root/04_amber_md/11_analysis"
for f in "$prmtop" "$traj"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done

python3 - "$prmtop" "$traj" "$outdir" <<'PY'
import sys, os, re, struct
import numpy as np
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
prmtop, traj, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
PAIRS={"break_C4_O3":("C4","O3"), "form_C6_C1":("C6","C1"), "carbox_C3_C10":("C3","C10")}
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
cha=[]
for ri in range(nres):
    if labels[ri]=="CHA":
        amap={names[a]:a for a in range(starts[ri],starts[ri+1])}; cha.append((ri+1,amap))
if not cha: sys.exit("FAIL no CHA residues found")
need=set(a for pr in PAIRS.values() for a in pr)
for rid,amap in cha:
    miss=[n for n in need if n not in amap]
    if miss: sys.exit("FAIL CHA %d missing atoms %s"%(rid,miss))
maxidx=max(a for _,amap in cha for a in amap.values())+1
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
def frame(i):
    off=cvar["begin"]+i*H["recsize"]
    return np.array(np.memmap(traj,dtype=">f4",mode="r",offset=off,shape=(natom,3))[:maxidx],float)
os.makedirs(outdir,exist_ok=True)
per={rid:{k:np.empty(H["frames"]) for k in list(PAIRS)+["rc"]} for rid,_ in cha}
for i in range(H["frames"]):
    x=frame(i)
    for rid,amap in cha:
        d={}
        for k,(a,b) in PAIRS.items():
            d[k]=float(np.linalg.norm(x[amap[a]]-x[amap[b]])); per[rid][k][i]=d[k]
        per[rid]["rc"][i]=d["break_C4_O3"]-d["form_C6_C1"]
with open(os.path.join(outdir,"rxn_coord_per_frame.dat"),"w") as f:
    f.write("# frame"+"".join("\tcha%d_break\tcha%d_form\tcha%d_rc"%(j+1,j+1,j+1) for j in range(len(cha)))+"\n")
    for i in range(H["frames"]):
        row="%d"%i
        for rid,_ in cha: row+="\t%.3f\t%.3f\t%.3f"%(per[rid]["break_C4_O3"][i],per[rid]["form_C6_C1"][i],per[rid]["rc"][i])
        f.write(row+"\n")
print("STEP 11c - substrate reaction-coordinate tracking (Claeyssens r)")
print("  CHA copies %d   frames %d"%(len(cha),H["frames"]))
for j,(rid,_) in enumerate(cha):
    b=per[rid]["break_C4_O3"]; fo=per[rid]["form_C6_C1"]; rc=per[rid]["rc"]; cx=per[rid]["carbox_C3_C10"]
    print("  CHA#%d (res %d):"%(j+1,rid))
    print("    break C4-O3   mean %.3f A  [%.3f..%.3f]  (bonded, should stay ~1.4)"%(b.mean(),b.min(),b.max()))
    print("    form  C6-C1   mean %.3f A  min %.3f (most near-attack)  max %.3f"%(fo.mean(),fo.min(),fo.max()))
    print("    r=break-form  mean %.3f A  min %.3f (closest to TS ~-0.5)  max %.3f"%(rc.mean(),rc.min(),rc.max()))
    print("    carbox C3-C10 mean %.3f A"%cx.mean())
print("  wrote", os.path.join(outdir,"rxn_coord_per_frame.dat"))
PY
echo "STEP 11c DONE"

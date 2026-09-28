#!/usr/bin/env bash
# Step 12a - QM/MM frame selection. Claeyssens-style: N competent reactant frames for
# multiple QM/MM paths (-> barrier distribution). competent = Arg90-O13 intact
# (nearest-Arg-to-O3 < CONTACT_CUT) AND near-attack (form C6-C1 < NAC_CUT). Only sites
# that are competent in >= SITE_MIN of frames are used as reactant sources (excludes
# sites whose substrate is only transiently bound). Bin-centre equal-interval per site.
# Args: [stride=5] [N=12] [contact_cut=3.2] [nac_cut=3.7] [site_min=0.5]
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
root="$HOME/system_development"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"
traj="$root/04_amber_md/10c_production/prod.nc"
outdir="$root/05_qmmm/12_frame_selection"; mkdir -p "$outdir"
for f in "$prmtop" "$traj"; do [[ -s "$f" ]] || { echo "FAIL missing $f"; exit 1; }; done

python3 - "$prmtop" "$traj" "$outdir" "${1:-5}" "${2:-12}" "${3:-3.2}" "${4:-3.7}" "${5:-0.5}" <<'PY'
import sys, os, re, struct
import numpy as np
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
prmtop, traj, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
stride=int(sys.argv[4]); N=int(sys.argv[5]); CONTACT_CUT=float(sys.argv[6]); NAC_CUT=float(sys.argv[7]); SITE_MIN=float(sys.argv[8])
ARG_N={"NE","NH1","NH2"}
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
cha=[]
for ri in range(nres):
    if labels[ri]=="CHA":
        amap={names[a]:a for a in range(starts[ri],starts[ri+1])}; cha.append((ri,amap))
if not cha: sys.exit("FAIL no CHA")
maxidx=int(max(argN.max(), max(a for _,m in cha for a in m.values())))+1
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
def mind(a,pts,box):
    d=pts-a; d-=np.round(d/box)*box; return float(np.sqrt((d*d).sum(1)).min())
scan=list(range(0,H["frames"],stride))
comp={j:[] for j in range(len(cha))}; total={j:0 for j in range(len(cha))}
for i in scan:
    x,box=frame(i)
    for j,(ri,amap) in enumerate(cha):
        form=float(np.linalg.norm(x[amap["C6"]]-x[amap["C1"]]))
        brk =float(np.linalg.norm(x[amap["C4"]]-x[amap["O3"]]))
        contact=mind(x[amap["O3"]], x[argN], box); total[j]+=1
        if contact<CONTACT_CUT and form<NAC_CUT: comp[j].append((i,form,brk-form,contact))
frac={j:len(comp[j])/max(1,total[j]) for j in comp}
qualifying=[j for j in comp if frac[j]>=SITE_MIN]
counts={j:len(comp[j]) for j in qualifying}; tot=sum(counts.values()); sel=[]
if tot>0:
    alloc={j:max(1,round(N*counts[j]/tot)) for j in qualifying}
    while sum(alloc.values())>N: alloc[max(alloc,key=lambda k:alloc[k])]-=1
    while sum(alloc.values())<N and any(counts[j]>alloc[j] for j in qualifying):
        alloc[max(qualifying,key=lambda k:counts[k]-alloc[k])]+=1
    for j in qualifying:
        pool=comp[j]; k=alloc[j]
        if k<=0 or not pool: continue
        idx=(((np.arange(k)+0.5)/k)*len(pool)).astype(int).clip(0,len(pool)-1)   # bin centres
        for t in idx: sel.append((j,)+pool[t])
sel.sort(key=lambda s:s[1])
mani=os.path.join(outdir,"selection_manifest.tsv")
with open(mani,"w") as f:
    f.write("# idx\tframe\tprod_ps\tsite\tCHA_resid\tform_C6_C1\tr\tArg90_O13\n")
    for n,(j,fr,form,r,contact) in enumerate(sel):
        f.write("%d\t%d\t%.1f\tCHA#%d\t%d\t%.3f\t%.3f\t%.3f\n"%(n+1,fr,float(fr),j+1,cha[j][0]+1,form,r,contact))
print("STEP 12a - frame selection  [stride %d, %d scanned, N %d, site_min %.0f%%]"%(stride,len(scan),N,SITE_MIN*100))
print("  filter: Arg90-O13 < %.1f A AND form C6-C1 < %.1f A"%(CONTACT_CUT,NAC_CUT))
for j in comp:
    q="QUALIFIES" if j in qualifying else "excluded (< site_min)"
    print("  CHA#%d (res %d): competent %d / %d (%.0f%%)  %s"%(j+1,cha[j][0]+1,len(comp[j]),total[j],100*frac[j],q))
print("  selected %d frames from %d qualifying site(s) -> %s"%(len(sel),len(qualifying),mani))
for n,(j,fr,form,r,contact) in enumerate(sel):
    print("    #%2d  frame %5d (%.0f ps)  CHA#%d  form %.2f  r %+.2f  Arg90-O13 %.2f"%(n+1,fr,float(fr),j+1,form,r,contact))
PY
echo "STEP 12a DONE"

#!/usr/bin/env bash
# Step 10c audit - verify the 20 ns production: reached nstlim, prod.nc holds the
# expected 20,000 frames (all atoms, finite), T/density stable, valid restart.
# Pure stdlib + numpy (no scipy/netCDF4/cpptraj). Single-threaded to avoid the
# login-node OpenBLAS RLIMIT_NPROC segfault.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
root="$HOME/system_development"
rundir="$root/04_amber_md/10c_production"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"

python3 - "$rundir" "$prmtop" <<'PY'
import sys, os, re, struct
import numpy as np
os.environ.setdefault("OPENBLAS_NUM_THREADS","1")
rundir, prmtop = sys.argv[1], sys.argv[2]
def die(m): print("FAIL "+m); sys.exit(1)
J=os.path.join

# ---- expected atoms from prmtop; expected frames from prod.in ----
def prmtop_natom(p):
    flags,cur={},None
    for ln in open(p,errors="replace"):
        ln=ln.rstrip("\n")
        if ln.startswith("%FLAG"): cur=ln.split()[1]; flags[cur]=[]
        elif ln.startswith(("%FORMAT","%")): continue
        elif cur is not None: flags[cur].append(ln)
    t=[]
    for ln in flags.get("POINTERS",[]):
        t+=[ln[i:i+8].strip() for i in range(0,len(ln.rstrip()),8) if ln[i:i+8].strip()]
    return int(t[0])
expected_atoms=prmtop_natom(prmtop)
pin=open(J(rundir,"prod.in")).read()
nstlim=int(re.search(r"nstlim\s*=\s*(\d+)",pin).group(1))
ntwx=int(re.search(r"ntwx\s*=\s*(\d+)",pin).group(1))
expected_frames=nstlim//ntwx

# ---- NetCDF3 header reader (stdlib struct) ----
STREAMING=0xFFFFFFFF; NC_DIM,NC_VAR,NC_ATT=0x0A,0x0B,0x0C
TS={1:1,2:1,3:2,4:4,5:4,6:8}
def read_header(path):
    data=open(path,"rb").read(1<<20); fsz=os.path.getsize(path)
    if data[:3]!=b"CDF": die("prod.nc not NetCDF classic (magic %r)"%data[:3])
    ver=data[3]; off64=(ver==2); pos=[4]
    def u32():
        v=struct.unpack_from(">I",data,pos[0])[0]; pos[0]+=4; return v
    def offset():
        if off64: v=struct.unpack_from(">Q",data,pos[0])[0]; pos[0]+=8
        else: v=struct.unpack_from(">I",data,pos[0])[0]; pos[0]+=4
        return v
    def nm():
        n=u32(); s=data[pos[0]:pos[0]+n].decode("ascii","replace"); pos[0]+=n+((4-(n%4))%4); return s
    numrecs=u32()
    dims=[]; tag=u32(); ne=u32()
    if tag==NC_DIM:
        for _ in range(ne): dn=nm(); dl=u32(); dims.append((dn,dl))
    def skip_att():
        t=u32(); n=u32()
        if t==NC_ATT:
            for _ in range(n):
                nm(); tp=u32(); k=u32(); nb=k*TS.get(tp,1); pos[0]+=nb+((4-(nb%4))%4)
    skip_att()
    recvars=[]; tag=u32(); nv=u32()
    if tag==NC_VAR:
        for _ in range(nv):
            vn=nm(); nd=u32(); dimids=[u32() for _ in range(nd)]; skip_att()
            tp=u32(); vsize=u32(); begin=offset()
            isrec=(nd>0 and dims[dimids[0]][1]==0)
            if isrec: recvars.append({"name":vn,"vsize":vsize,"begin":begin})
    dimd={n:l for n,l in dims}
    recsize=sum(v["vsize"] for v in recvars)
    if numrecs!=STREAMING: frames=numrecs
    else:
        first=min((v["begin"] for v in recvars),default=fsz)
        frames=(fsz-first)//recsize if recsize else 0
    return {"ver":ver,"fsz":fsz,"numrecs":numrecs,"dims":dimd,"recvars":recvars,"recsize":recsize,"frames":frames}

nc=J(rundir,"prod.nc")
if not os.path.exists(nc): die("missing prod.nc")
h=read_header(nc)
atoms=h["dims"].get("atom")
frames=h["frames"]
# filesize cross-check
first=min((v["begin"] for v in h["recvars"]),default=h["fsz"])
implied=(h["fsz"]-first)//h["recsize"] if h["recsize"] else 0

# finite spot-check of first & last coordinate frames
def coord_finite(fidx):
    cv=next((v for v in h["recvars"] if v["name"]=="coordinates"),None)
    if cv is None or atoms is None: return None
    nsp=h["dims"].get("spatial",3); off=cv["begin"]+fidx*h["recsize"]
    mm=np.memmap(nc,dtype=">f4",mode="r",offset=off,shape=(atoms*nsp,))
    return bool(np.isfinite(np.asarray(mm)).all())
f0=coord_finite(0); fN=coord_finite(frames-1)

# ---- prod.out: last NSTEP + final AVERAGES temp/density ----
out=open(J(rundir,"prod.out"),errors="replace").read()
if "nan" in out.lower(): die("prod.out contains NaN")
nsteps=[int(x) for x in re.findall(r"NSTEP\s*=\s*(\d+)",out)]
last_nstep=max(nsteps) if nsteps else -1
am=re.search(r"A V E R A G E S.*?(?=R M S|\Z)",out,re.S)
avgT=avgD=None
if am:
    mt=re.search(r"TEMP\(K\)\s*=\s*([\d.]+)",am.group(0)); avgT=float(mt.group(1)) if mt else None
    md=re.search(r"Density\s*=\s*([\d.]+)",am.group(0)); avgD=float(md.group(1)) if md else None

# ---- prod.rst7 ----
rl=open(J(rundir,"prod.rst7"),errors="replace").read().splitlines()
rn=int(rl[1].split()[0]); fl=[]
for ln in rl[2:]:
    for t in ln.split():
        try: fl.append(float(t.replace("D","E")))
        except: pass
nc_f=len(fl); has_box=nc_f in (3*rn+6,6*rn+6); has_vel=nc_f in (6*rn,6*rn+6)
nan_r=bool(re.search(r"nan|inf|\*\*\*","\n".join(rl[2:]),re.I))

print("STEP 10c audit - 20 ns production")
print("  prod.nc: version %d, frames %d (expected %d), atoms %s (expected %d)" % (h["ver"],frames,expected_frames,atoms,expected_atoms))
print("           filesize %.2f GB, recsize %d B/frame, frames-by-filesize %d" % (h["fsz"]/1e9,h["recsize"],implied))
print("           coords finite: frame0 %s, frame%d %s" % (f0,frames-1,fN))
print("  prod.out: last NSTEP %d (nstlim %d), final avg TEMP %s K, Density %s g/cm^3" % (last_nstep,nstlim,avgT,avgD))
print("  prod.rst7: NATOM %d (== prmtop %d), floats %d, velocities %s, box %s, finite %s" % (rn,expected_atoms,nc_f,has_vel,has_box,not nan_r))

if atoms!=expected_atoms: die("prod.nc atom dim %s != prmtop %d"%(atoms,expected_atoms))
if frames!=expected_frames: die("prod.nc frames %d != expected %d"%(frames,expected_frames))
if implied not in (frames, frames):  # filesize must corroborate frame count within 1
    if abs(implied-frames)>1: die("filesize implies %d frames, header says %d"%(implied,frames))
if f0 is False or fN is False: die("prod.nc has non-finite coordinates")
if last_nstep!=nstlim: die("prod.out last NSTEP %d != nstlim %d (run did not complete)"%(last_nstep,nstlim))
if avgT is None or not (295<=avgT<=305): die("final avg temp %s outside 295-305 K"%avgT)
if avgD is None or not (0.98<=avgD<=1.05): die("final avg density %s outside 0.98-1.05"%avgD)
if rn!=expected_atoms: die("prod.rst7 NATOM %d != prmtop %d"%(rn,expected_atoms))
if not has_vel: die("prod.rst7 has no velocities")
if not has_box: die("prod.rst7 missing box")
if nan_r: die("prod.rst7 contains NaN/Inf")
print("\n  RESULT: PASS - 20 ns complete, 20,000 frames, T/density stable, restart valid")
PY
echo "STEP 10c AUDIT DONE"

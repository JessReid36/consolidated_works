#!/usr/bin/env bash
# Step 10b audit - verify the five-stage NPT equilibration: density climbs to ~1.0,
# T~300 K, restraint ramped 200->0, final restart live (velocities + box, all atoms).
set -euo pipefail
root="$HOME/system_development"
rundir="$root/04_amber_md/10b_equil"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"

python3 - "$rundir" "$prmtop" <<'PY'
import sys, re
from pathlib import Path
rundir=Path(sys.argv[1]); prmtop=Path(sys.argv[2])
def die(m): print("FAIL "+m); sys.exit(1)
def fnum(s):
    s=s.strip()
    if "*" in s or s.lower() in ("nan","inf","-inf"): return None
    try: return float(s.replace("D","E"))
    except: return None
def prmtop_natom(p):
    flags,cur={},None
    for ln in p.read_text(errors="replace").splitlines():
        if ln.startswith("%FLAG"): cur=ln.split()[1]; flags[cur]=[]
        elif ln.startswith(("%FORMAT","%")): continue
        elif cur is not None: flags[cur].append(ln)
    t=[]
    for ln in flags.get("POINTERS",[]):
        t+=[ln[i:i+8].strip() for i in range(0,len(ln.rstrip()),8) if ln[i:i+8].strip()]
    return int(t[0])
expected=prmtop_natom(prmtop)

def avg_block(outfile):
    txt=outfile.read_text(errors="replace")
    if "nan" in txt.lower(): die("%s contains NaN" % outfile.name)
    m=re.search(r"A V E R A G E S.*?(?=R M S|\Z)", txt, re.S)
    if not m: die("%s: no AVERAGES block (stage did not finish)" % outfile.name)
    b=m.group(0)
    temp=re.search(r"TEMP\(K\)\s*=\s*(-?[\d.]+)", b)
    dens=re.search(r"Density\s*=\s*(-?[\d.]+)", b)
    press=re.search(r"PRESS\s*=\s*(-?[\d.]+)", b)
    return (fnum(temp.group(1)) if temp else None,
            fnum(press.group(1)) if press else None,
            fnum(dens.group(1)) if dens else None)

rows=[]
for i in range(1,6):
    f=rundir/("npt%d.out"%i)
    if not f.is_file(): die("missing npt%d.out"%i)
    t,p,d=avg_block(f); rows.append((i,t,p,d))

ladder=[]
for i in range(1,6):
    mdin=rundir/("npt%d.in"%i); w=None
    if mdin.is_file():
        s=mdin.read_text()
        mt=re.search(r"restraint_wt\s*=\s*([\d.]+)", s)
        w=0.0 if re.search(r"ntr\s*=\s*0", s) else (fnum(mt.group(1)) if mt else None)
    ladder.append(w)

dens=[d for _,_,_,d in rows]
if any(d is None for d in dens): die("could not read density from every stage")
mono=all(dens[i] <= dens[i+1]+0.01 for i in range(len(dens)-1))
final_d=dens[-1]; final_t=rows[-1][1]

rst=rundir/"npt5.rst7"
if not rst.is_file(): die("missing npt5.rst7")
rl=rst.read_text(errors="replace").splitlines()
rst_nat=int(rl[1].split()[0]); floats=[]
for ln in rl[2:]:
    for tok in ln.split():
        try: floats.append(float(tok.replace("D","E")))
        except: pass
nc=len(floats)
has_box=nc in (3*rst_nat+6, 6*rst_nat+6)
has_vel=nc in (6*rst_nat, 6*rst_nat+6)
ok_count=nc in (3*rst_nat,3*rst_nat+6,6*rst_nat,6*rst_nat+6)
nan_c=bool(re.search(r"nan|inf|\*\*\*","\n".join(rl[2:]),re.I))

print("STEP 10b audit - five-stage NPT equilibration")
print("  stage averages (restraint -> TEMP / PRESS / Density):")
for (i,t,p,d),w in zip(rows,ladder):
    print("    npt%d  restraint=%s  TEMP=%s K  PRESS=%s  Density=%.4f g/cm^3" % (i,("%.1f"%w if w is not None else "?"),t,p,d))
print("  density monotonic up: %s   final density: %.4f g/cm^3" % (mono, final_d))
print("  final avg temperature: %s K" % final_t)
print("  restraint ladder from mdin: %s" % ladder)
print("  npt5.rst7: NATOM %d (== prmtop %d), floats %d, velocities %s, box %s, finite %s" % (
    rst_nat, expected, nc, has_vel, has_box, not nan_c))

if rst_nat!=expected: die("npt5.rst7 NATOM %d != prmtop %d" % (rst_nat,expected))
if not ok_count: die("npt5.rst7 float count %d inconsistent with NATOM %d" % (nc,rst_nat))
if nan_c: die("npt5.rst7 contains NaN/Inf")
if not has_box: die("npt5.rst7 missing box")
if not has_vel: die("npt5.rst7 has no velocities")
if ladder!=[200.0,100.0,50.0,10.0,0.0]: die("restraint ladder %s != 200/100/50/10/0" % ladder)
if not (0.97<=final_d<=1.03): die("final density %.4f outside 0.97-1.03" % final_d)
if final_t is None or not (290<=final_t<=310): die("final avg temp %s outside 290-310 K" % final_t)
if not mono: die("density not monotonically increasing")
print("\n  RESULT: PASS - density equilibrated to ~1.0, T~300 K, restraint relaxed to 0, restart live")
PY
echo "STEP 10b AUDIT DONE"

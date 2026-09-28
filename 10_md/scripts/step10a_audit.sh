#!/usr/bin/env bash
# Step 10a audit - verify the restrained minimisation: finished, energy finite and
# lowered, no atomic overlap/NaN, restart intact (all atoms, box present).
set -euo pipefail
root="$HOME/system_development"
rundir="$root/04_amber_md/10a_min"
prmtop="$root/03_amber/tleap_build/complex_solvated.prmtop"

python3 - "$rundir" "$prmtop" <<'PY'
import sys, re
from pathlib import Path

rundir = Path(sys.argv[1]); prmtop = Path(sys.argv[2])
out = rundir/"10a_min.out"; rst = rundir/"10a_min.rst7"
def die(m): print("FAIL " + m); sys.exit(1)

def prmtop_natom(p):
    flags, cur = {}, None
    for ln in p.read_text(errors="replace").splitlines():
        if ln.startswith("%FLAG"): cur=ln.split()[1]; flags[cur]=[]
        elif ln.startswith(("%FORMAT","%")): continue
        elif cur is not None: flags[cur].append(ln)
    toks=[]
    for ln in flags.get("POINTERS",[]):
        toks += [ln[i:i+8].strip() for i in range(0,len(ln.rstrip()),8) if ln[i:i+8].strip()]
    return int(toks[0])
expected_natom = prmtop_natom(prmtop)

if not out.is_file(): die("missing 10a_min.out")
lines = out.read_text(errors="replace").splitlines()
text = "\n".join(lines)
hdr = re.compile(r"NSTEP\s+ENERGY\s+RMS\s+GMAX")
def fnum(s):
    s=s.strip()
    if "*" in s or s.lower() in ("nan","inf","-inf","infinity"): return None
    try: return float(s.replace("D","E"))
    except ValueError: return None

records=[]; final_next=False; seen_final = "FINAL RESULTS" in text; i=0
while i < len(lines):
    if "FINAL RESULTS" in lines[i]: final_next=True; i+=1; continue
    if hdr.search(lines[i]):
        j=i+1
        while j<len(lines) and not lines[j].strip(): j+=1
        if j>=len(lines): break
        parts=lines[j].split()
        rec={"nstep":parts[0] if parts else "","energy":fnum(parts[1]) if len(parts)>1 else None,
             "rms":fnum(parts[2]) if len(parts)>2 else None,"gmax":fnum(parts[3]) if len(parts)>3 else None,
             "final":final_next,"comp":{}}
        window="\n".join(lines[j+1:j+12])
        for key in ["BOND","ANGLE","DIHED","UB","IMP","CMAP","VDWAALS","EEL","EGB",
                    "1-4 VDW","1-4 EEL","RESTRAINT","HBOND","EELEC","EKtot","EPtot"]:
            patt=key.replace("-",r"\-").replace(" ",r"\s+")
            m=re.search(patt+r"\s*=\s*(-?\d+\.?\d*(?:[EeDd][+\-]?\d+)?|\*+)", window)
            if m: rec["comp"][key]=fnum(m.group(1))
        records.append(rec); final_next=False; i=j+1; continue
    i+=1

if not seen_final: die("no FINAL RESULTS block - minimisation did not finish")
finals=[r for r in records if r["final"]]
if not finals: die("FINAL RESULTS present but no parseable final energy line")
first=records[0]; last=finals[-1]
if last["energy"] is None: die("final ENERGY is non-numeric (NaN/overflow)")
for key,v in last["comp"].items():
    if v is None: die("final energy component %s is non-numeric (NaN/overflow)" % key)
vdw=last["comp"].get("VDWAALS")
if vdw is not None and abs(vdw) > 1e6: die("VDWAALS %.3e implausibly large (atomic overlap?)" % vdw)
decreased=(first["energy"] is not None) and (last["energy"] < first["energy"])

if not rst.is_file(): die("missing 10a_min.rst7")
rl=rst.read_text(errors="replace").splitlines()
rst_natom=int(rl[1].split()[0]); floats=[]
for ln in rl[2:]:
    for t in ln.split():
        try: floats.append(float(t.replace("D","E")))
        except ValueError: pass
nan_coords=bool(re.search(r"nan|inf|\*\*\*","\n".join(rl[2:]),re.I))
ncoord=len(floats); has_box=ncoord==3*rst_natom+6; coord_ok=has_box or ncoord==3*rst_natom
if rst_natom!=expected_natom: die("rst7 NATOM %d != prmtop NATOM %d" % (rst_natom,expected_natom))
if not coord_ok: die("rst7 coord count %d != 3*%d (+6 box)" % (ncoord,rst_natom))
if nan_coords: die("rst7 contains NaN/Inf coordinates")
if not has_box: die("rst7 missing box (expected periodic system)")

print("STEP 10a audit - restrained minimisation")
print("  minimisation finished: FINAL RESULTS present")
print("  energy: initial %s -> final %.4f kcal/mol  (decreased: %s)" % (
    ("%.4f"%first["energy"]) if first["energy"] is not None else "?", last["energy"], decreased))
print("  final RMS gradient %s  GMAX %s (kcal/mol/A)" % (last["rms"], last["gmax"]))
c=last["comp"]
print("  components: BOND %s ANGLE %s DIHED %s VDWAALS %s EEL %s 1-4VDW %s RESTRAINT %s" % (
    c.get("BOND"),c.get("ANGLE"),c.get("DIHED"),c.get("VDWAALS"),c.get("EEL"),c.get("1-4 VDW"),c.get("RESTRAINT")))
print("  rst7: NATOM %d (== prmtop %d), coords %d, box %s, finite %s" % (
    rst_natom, expected_natom, ncoord, has_box, not nan_coords))
if not decreased: die("energy did not decrease (initial %.4f -> final %.4f)" % (first["energy"], last["energy"]))
print("\n  RESULT: PASS - minimisation completed, energy lowered, no overlap/NaN, restart intact")
PY
echo "STEP 10a AUDIT DONE"

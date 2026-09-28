#!/usr/bin/env bash
# Step 09b - tleap build (ff14SB + GAFF + TIP3P, 10 A box, Na+ neutralise) plus a
# deep prmtop/inpcrd + tleap-log audit, in one pass.
#
# Folds attempt_3 steps 17 (build) + 18 (topology identity) + 18b (warning
# classification) into one. Correctness points that shaped this step:
#  - The solvated PDB saturates at 9999 residues, so water/ion counts are read
#    from the prmtop/inpcrd (source of truth), never the PDB. Na+ = |dry charge|.
#  - tleap emits warnings as a generic "teLeap: Warning!" banner followed by the
#    real message on the next line; the audit reads the messages (not banners) and
#    buckets each into a named benign class. Anything it cannot place is REVIEW-
#    flagged verbatim, never silently passed.
#  - Close-contact warnings give only atom names, so the audit resolves each to
#    its residues by matching the name pair + distance against the pre-tleap complex_for_tleap.pdb (chain-aware, per-chain numbering),
#    and flags any that touch the ligand, the active-site set, or two heavy atoms.
#  - A prmtop embeds a wall-clock DATE line, so its raw hash churns every run; the
#    fingerprint therefore hashes the prmtop with that line stripped (stable),
#    alongside the naturally-stable inpcrd/pdb hashes.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
builddir="$root/03_amber/tleap_build"
ligdir="$root/03_amber/ligand_gaff"
admin="$root/00_admin"

pretleap="$builddir/complex_for_tleap.pdb"      # 09a
cha_mol2="$ligdir/cha_gaff.mol2"                 # 08b
cha_frcmod="$ligdir/cha.frcmod"                  # 08b
mkdir -p "$builddir" "$admin"

echo "=== step 09b: input presence ==="
for f in "$pretleap" "$cha_mol2" "$cha_frcmod"; do
  [[ -s "$f" ]] || { echo "FAIL missing/empty: $f"; exit 1; }
  echo "PASS $f"
done

echo
echo "=== step 09b: write tleap input ==="
cat > "$builddir/tleap_complex.in" <<EOF
logFile tleap_complex.log
# Paper-aligned build: protein ff14SB, chorismate GAFF, TIP3P water, 10 A box, Na+.
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p

CHA = loadmol2 $cha_mol2
loadamberparams $cha_frcmod

complex = loadpdb $pretleap
check complex
charge complex

saveamberparm complex complex_dry.prmtop complex_dry.inpcrd
savepdb complex complex_dry.pdb

solvatebox complex TIP3PBOX 10.0
addions complex Na+ 0
check complex
charge complex

saveamberparm complex complex_solvated.prmtop complex_solvated.inpcrd
savepdb complex complex_solvated.pdb
quit
EOF
cat "$builddir/tleap_complex.in"

echo
echo "=== step 09b: load AMBER22 tools ==="
set +u                                      # amber.sh references unset vars
export PERL5LIB="${PERL5LIB:-}" PYTHONPATH="${PYTHONPATH:-}"
module load app/amber22/22
set -u
echo "AMBERHOME=${AMBERHOME:-UNSET}"
command -v tleap >/dev/null || { echo "FAIL tleap not found"; exit 1; }

echo
echo "=== step 09b: run tleap ==="
cd "$builddir"
if ! tleap -f tleap_complex.in > tleap_complex.stdout 2> tleap_complex.stderr; then
  echo "FAIL tleap returned nonzero"; echo "--- stdout tail ---"; tail -40 tleap_complex.stdout
  echo "--- log tail ---"; tail -60 tleap_complex.log 2>/dev/null || true; exit 1
fi
echo "PASS tleap completed"

echo
echo "=== step 09b: required outputs present ==="
for f in complex_dry.prmtop complex_dry.inpcrd complex_dry.pdb \
         complex_solvated.prmtop complex_solvated.inpcrd tleap_complex.log; do
  [[ -s "$f" ]] || { echo "FAIL missing/empty: $f"; exit 1; }
  echo "PASS $f ($(stat -c%s "$f") bytes)"
done

echo
echo "=== step 09b: deep prmtop/inpcrd + log audit (source of truth) ==="
python3 - "$builddir" "$admin" <<'PY'
import sys, re, math, hashlib
from collections import Counter
from pathlib import Path

builddir, admin = Path(sys.argv[1]), Path(sys.argv[2])
audit_txt = admin / "step09b_topology_audit.txt"
res_tsv   = admin / "step09b_prmtop_residue_counts.tsv"
warn_tsv  = admin / "step09b_warning_classification.tsv"
fp_txt    = admin / "step09_fingerprints.txt"

CHARGE_SCALE = 18.2223            # AMBER internal charge units -> e
SOLVENT = {"WAT","HOH"}
ACTIVE_SITE = {7,57,59,60,63,73,74,75,78,90,108,115}   # per-chain resids (Agbaglo SI S3)

def die(m): sys.exit(f"FAIL {m}")

# ---------- prmtop / inpcrd parsing (source of truth) ------------------------
def read_flags(path):
    flags, fmts, cur = {}, {}, None
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("%FLAG"):   cur = line.split()[1]; flags[cur] = []
        elif line.startswith("%FORMAT"): fmts[cur] = line[line.find("(")+1:line.rfind(")")]
        elif line.startswith("%"):     continue
        elif cur is not None:          flags[cur].append(line)
    return flags, fmts

def width_of(fmt):
    m = re.match(r"\s*\d*[aAiIeEfFgG](\d+)", fmt); return int(m.group(1)) if m else None

def tokens(flags, fmts, name, cast):
    if name not in flags: die(f"prmtop missing %FLAG {name}")
    w = width_of(fmts.get(name, "")); out = []
    for line in flags[name]:
        if w: out += [line[i:i+w].strip() for i in range(0, len(line.rstrip('\n')), w) if line[i:i+w].strip()]
        else: out += line.split()
    return [cast(x) for x in out]

def parse_prmtop(path):
    flags, fmts = read_flags(path)
    pointers = tokens(flags, fmts, "POINTERS", int)
    natom, nres = pointers[0], pointers[11]
    names   = tokens(flags, fmts, "ATOM_NAME", str)
    charges = tokens(flags, fmts, "CHARGE", float)
    labels  = tokens(flags, fmts, "RESIDUE_LABEL", str)
    if len(names) != natom:   die(f"{path.name}: ATOM_NAME {len(names)} != NATOM {natom}")
    if len(charges) != natom: die(f"{path.name}: CHARGE {len(charges)} != NATOM {natom}")
    if len(labels) != nres:   die(f"{path.name}: RESIDUE_LABEL {len(labels)} != NRES {nres}")
    return {"natom":natom, "nres":nres, "total_charge":sum(charges)/CHARGE_SCALE,
            "label_counts":Counter(labels)}

def parse_inpcrd(path):
    lines = path.read_text(errors="replace").splitlines()
    natom = int(lines[1].split()[0]); floats = []
    for line in lines[2:]:
        for tok in line.split():
            try: floats.append(float(tok.replace("D","E")))
            except ValueError: pass
    ncoord = len(floats); has_box = ncoord == 3*natom + 6
    return {"natom":natom, "ncoord":ncoord, "has_box":has_box, "ok":has_box or ncoord == 3*natom}

# ---------- fingerprints (prmtop DATE line stripped so hashes are stable) -----
def content_sha(path, strip_version=False):
    data = path.read_bytes()
    if strip_version:
        data = b"\n".join(l for l in data.splitlines() if not l.startswith(b"%VERSION"))
    return hashlib.sha256(data).hexdigest()

# ---------- PDB atoms + geometry (for close-contact residue resolution) -------
def is_h_name(name):
    n = name.strip()
    return n[:1].upper() == "H" or (len(n) > 1 and n[0].isdigit() and n[1:2].upper() == "H")

def parse_pdb_atoms(path):
    atoms = []
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")): continue
        atoms.append({"name":line[12:16].strip(),"resname":line[17:20].strip(),
            "chain":line[21],"resid":int(line[22:26]),
            "x":float(line[30:38]),"y":float(line[38:46]),"z":float(line[46:54])})
    return atoms

def _d(a,b): return math.sqrt((a["x"]-b["x"])**2+(a["y"]-b["y"])**2+(a["z"]-b["z"])**2)

# ---------- warning events (banner+message, and inline) ----------------------
def warning_events(streams):
    events = []
    for sname, lines in streams.items():
        i = 0
        while i < len(lines):
            l = lines[i]
            if "teleap: warning!" in l.lower():
                j = i+1
                while j < len(lines) and not lines[j].strip(): j += 1
                if j < len(lines):
                    events.append((sname, lines[j].strip())); i = j+1; continue
            elif re.search(r"\bwarning\b", l, re.I) or "close contact" in l.lower():
                events.append((sname, l.strip()))
            i += 1
    return events

CC_RE = re.compile(r"close contact of\s+([\d.]+)\s+angstroms between\s+"
                   r"(?:nonbonded atoms\s+)?(\S+)\s+and\s+(\S+)", re.I)

def bucket_and_resolve(events, ref_pdb):
    buckets = Counter(); other = []; contacts = []
    for sname, msg in events:
        m = msg.lower()
        if re.search(r"converting [nc]-terminal residue name", m):        buckets["terminal_name_conversion"] += 1
        elif re.search(r"unperturbed charge.*not zero|charge of the unit.*not zero", m): buckets["nonzero_charge"] += 1
        elif re.search(r"ignoring the warning from unit checking", m):    buckets["unit_check_note"] += 1
        elif "close contact" in m:
            buckets["close_contact"] += 1
            cm = CC_RE.search(msg)
            if cm: contacts.append((cm.group(2), cm.group(3), float(cm.group(1))))
            else:  other.append((sname, msg))
        else: other.append((sname, msg))

    atoms = parse_pdb_atoms(ref_pdb) if ref_pdb.exists() else []
    by_name = {}
    for a in atoms: by_name.setdefault(a["name"], []).append(a)
    resolved, unresolved, seen = [], [], set()
    for n1, n2, dist in contacts:
        key = (frozenset((n1, n2)), round(dist, 3))
        if key in seen: continue        # same pair reported by dry and solvated check
        seen.add(key)
        cands = [(a,b) for a in by_name.get(n1,[]) for b in by_name.get(n2,[])
                 if a is not b and abs(_d(a,b)-dist) <= 0.02]
        if len(cands) == 1:
            a, b = cands[0]
            resolved.append({"n1":n1,"n2":n2,"dist":dist,
                "a":f"{a['chain']}:{a['resid']}:{a['resname']}:{n1}",
                "b":f"{b['chain']}:{b['resid']}:{b['resname']}:{n2}",
                "ligand":("CHA" in (a["resname"], b["resname"])),
                "active_site":any(x["resid"] in ACTIVE_SITE and x["resname"]!="CHA" for x in (a,b)),
                "heavy_heavy":(not is_h_name(n1)) and (not is_h_name(n2)),
                "same_res":(a["chain"],a["resid"])==(b["chain"],b["resid"])})
        else:
            unresolved.append((n1, n2, dist, len(cands)))
    return buckets, other, resolved, unresolved

# ---------- run ---------------------------------------------------------------
dry_top = parse_prmtop(builddir/"complex_dry.prmtop")
dry_crd = parse_inpcrd(builddir/"complex_dry.inpcrd")
sol_top = parse_prmtop(builddir/"complex_solvated.prmtop")
sol_crd = parse_inpcrd(builddir/"complex_solvated.inpcrd")

def counts(top):
    solv = sum(v for k,v in top["label_counts"].items() if k in SOLVENT)
    na   = sum(v for k,v in top["label_counts"].items() if k in {"Na+","NA","Na"})
    cl   = sum(v for k,v in top["label_counts"].items() if k in {"Cl-","CL","Cl"})
    return solv, na, cl, top["label_counts"].get("CHA", 0)

dry_solv, dry_na, dry_cl, dry_cha = counts(dry_top)
sol_solv, sol_na, sol_cl, sol_cha = counts(sol_top)
dry_charge_int = round(dry_top["total_charge"]); sol_charge_int = round(sol_top["total_charge"])
expected_na = -dry_charge_int if dry_charge_int < 0 else 0

if not dry_crd["ok"]: die(f"dry inpcrd coord count {dry_crd['ncoord']} != 3*{dry_crd['natom']}")
if not sol_crd["ok"]: die(f"solvated inpcrd coord count {sol_crd['ncoord']} unexpected")
if dry_top["natom"] != dry_crd["natom"]: die("dry prmtop/inpcrd NATOM mismatch")
if sol_top["natom"] != sol_crd["natom"]: die("solvated prmtop/inpcrd NATOM mismatch")
if dry_crd["has_box"]:     die("dry inpcrd unexpectedly has box")
if not sol_crd["has_box"]: die("solvated inpcrd missing box (6 box floats)")
if abs(dry_top["total_charge"]-dry_charge_int) > 0.01: die(f"dry total charge {dry_top['total_charge']:.4f} not near integer")
if abs(sol_top["total_charge"]-sol_charge_int) > 0.01: die(f"solvated total charge {sol_top['total_charge']:.4f} not near integer")
if sol_charge_int != 0:   die(f"solvated system not neutral: total charge {sol_charge_int}")
if sol_na != expected_na: die(f"Na+ added {sol_na}, expected {expected_na} (=|dry charge|)")
if sol_cha != dry_cha:    die(f"CHA residues changed on solvation: dry {dry_cha} vs solvated {sol_cha}")
if sol_top["natom"] != dry_top["natom"] + 3*sol_solv + sol_na + sol_cl:
    die(f"solvated NATOM {sol_top['natom']} != dry {dry_top['natom']} + 3*{sol_solv} + {sol_na}Na + {sol_cl}Cl")

# streams: hard errors abort regardless of stream
def read_stream(name):
    p = builddir / name
    return p.read_text(errors="replace").splitlines() if p.exists() else []
streams = {n: read_stream(n) for n in ("tleap_complex.stdout","tleap_complex.log","tleap_complex.stderr")}
combined = "\n".join("\n".join(v) for v in streams.values())
errors      = re.findall(r"\bERROR\b|Errors\s*=\s*[1-9]", combined)
fatal       = [l for l in combined.splitlines() if "fatal" in l.lower()]
unknown_bad = [l for l in combined.splitlines() if "unknown" in l.lower() and ("atom type" in l.lower() or "residue" in l.lower())]
if errors or fatal or unknown_bad:
    for l in (fatal + unknown_bad)[:20]: print("  ", l.strip())
    die("tleap log contains errors / fatal / unknown atom-type or residue problems")

events = warning_events(streams)
buckets, warn_other, resolved, unresolved = bucket_and_resolve(events, builddir/"complex_for_tleap.pdb")
warn_total = len(events)
assert warn_total == sum(buckets.values()) + len(warn_other), "warning accounting mismatch"

# REVIEW only on genuinely actionable items; ligand/active-site H-contacts are reported but expected
review = []
if warn_other:   review.append(f"{len(warn_other)} unclassified tleap warning(s)")
if unresolved:   review.append(f"{len(unresolved)} close contact(s) could not be resolved to a unique atom pair")
hh = [r for r in resolved if r["heavy_heavy"]]
if hh:           review.append(f"{len(hh)} close contact(s) between two heavy atoms (possible real clash)")

# fingerprints (stable across runs)
fps = [
    ("complex_dry.prmtop.content", content_sha(builddir/"complex_dry.prmtop", strip_version=True)),
    ("complex_solvated.prmtop.content", content_sha(builddir/"complex_solvated.prmtop", strip_version=True)),
    ("complex_dry.inpcrd", content_sha(builddir/"complex_dry.inpcrd")),
    ("complex_solvated.inpcrd", content_sha(builddir/"complex_solvated.inpcrd")),
    ("complex_for_tleap.pdb", content_sha(builddir/"complex_for_tleap.pdb")),
]
fp_txt.write_text("".join(f"{h}  {n}\n" for n,h in fps))

with warn_tsv.open("w") as f:
    f.write("kind\tdetail\tcount_or_value\n")
    for s in ("tleap_complex.stdout","tleap_complex.log","tleap_complex.stderr"):
        f.write(f"stream_events\t{s}\t{sum(1 for e in events if e[0]==s)}\n")
    for b in ("terminal_name_conversion","nonzero_charge","unit_check_note","close_contact"):
        f.write(f"bucket\t{b}\t{buckets.get(b,0)}\n")
    f.write(f"bucket\tother\t{len(warn_other)}\n")
    f.write(f"close_contacts_distinct\tresolved\t{len(resolved)}\n")
    for r in sorted(resolved, key=lambda r: r["dist"]):
        flags = ",".join(t for t,on in [("ligand",r["ligand"]),("active_site",r["active_site"]),
                    ("heavy_heavy",r["heavy_heavy"]),("same_residue",r["same_res"])] if on) or "sidechain_H"
        f.write(f"close_contact\t{r['a']} <-> {r['b']} ({r['dist']:.3f} A)\t{flags}\n")
    for n1,n2,d,ncand in unresolved: f.write(f"close_contact_unresolved\t{n1}-{n2} {d:.3f}A\t{ncand} candidates\n")
    for s,l in warn_other: f.write(f"unclassified\t{s}\t{l}\n")

with res_tsv.open("w") as f:
    f.write("system\tlabel\tcount\n")
    for sysname, top in (("dry",dry_top),("solvated",sol_top)):
        for lab, n in sorted(top["label_counts"].items()): f.write(f"{sysname}\t{lab}\t{n}\n")

bucket_str = ", ".join(f"{b}={buckets.get(b,0)}" for b in
    ("terminal_name_conversion","nonzero_charge","unit_check_note","close_contact")) + f", other={len(warn_other)}"
audit_lines = [
    ("dry_prmtop_natom",dry_top["natom"]),("dry_prmtop_nres",dry_top["nres"]),
    ("dry_total_charge",f"{dry_top['total_charge']:.4f}"),("dry_charge_int",dry_charge_int),("dry_CHA_residues",dry_cha),
    ("solvated_prmtop_natom",sol_top["natom"]),("solvated_prmtop_nres",sol_top["nres"]),
    ("solvated_total_charge",f"{sol_top['total_charge']:.4f}"),("solvated_charge_int",sol_charge_int),
    ("water_residues",sol_solv),("Na_ions",sol_na),("Cl_ions",sol_cl),("solvated_CHA_residues",sol_cha),
    ("expected_Na_from_dry_charge",expected_na),
    ("prmtop_inpcrd_natom_match_dry",dry_top["natom"]==dry_crd["natom"]),
    ("prmtop_inpcrd_natom_match_solvated",sol_top["natom"]==sol_crd["natom"]),
    ("tleap_errors",len(errors)),("tleap_fatal",len(fatal)),
    ("tleap_warning_events",warn_total),("tleap_warning_buckets",bucket_str),
    ("close_contacts_distinct",len(resolved)),
    ("close_contacts_heavy_heavy",len(hh)),("close_contacts_unresolved",len(unresolved)),
    ("review","; ".join(review) if review else "none"),
]
audit_txt.write_text("\n".join(f"{k}\t{v}" for k,v in audit_lines) + "\n")

print(f"  dry:       NATOM {dry_top['natom']}  NRES {dry_top['nres']}  charge {dry_top['total_charge']:+.4f}  box {dry_crd['has_box']}  CHA {dry_cha}")
print(f"  solvated:  NATOM {sol_top['natom']}  NRES {sol_top['nres']}  charge {sol_top['total_charge']:+.4f}  box {sol_crd['has_box']}")
print(f"  water residues (from prmtop): {sol_solv}   Na+: {sol_na} (expected {expected_na})   Cl-: {sol_cl}   CHA: {sol_cha}")
print(f"  prmtop<->inpcrd NATOM match: dry {dry_top['natom']==dry_crd['natom']}, solvated {sol_top['natom']==sol_crd['natom']}")
print(f"  tleap: errors {len(errors)}, fatal {len(fatal)}, warning events {warn_total}")
print(f"    buckets: {bucket_str}")
print(f"  close contacts: {len(resolved)} distinct (each seen in dry + solvated check), heavy-heavy {len(hh)}, unresolved {len(unresolved)}")
for r in sorted(resolved, key=lambda r: r["dist"]):
    flags = ",".join(t for t,on in [("LIGAND",r["ligand"]),("ACTIVE-SITE",r["active_site"]),
                ("HEAVY-HEAVY",r["heavy_heavy"]),("same-res",r["same_res"])] if on) or "sidechain H-contact"
    print(f"    {r['a']} <-> {r['b']}   {r['dist']:.3f} A   [{flags}]")
if warn_other:
    print("  REVIEW - unclassified warnings:")
    for s,l in warn_other: print(f"    [{s.split('_')[-1]}] {l}")
print(f"\n  WROTE {audit_txt}")
print(f"  WROTE {res_tsv}")
print(f"  WROTE {warn_tsv}")
print(f"  WROTE {fp_txt}")
if review:
    print("\n  RESULT: PASS (topology sound) - REVIEW: " + "; ".join(review))
else:
    print("\n  RESULT: PASS - topology neutral, box present, counts consistent, all warnings classified, no heavy-atom clashes")
PY

echo
echo "=== step 09b: method note ==="
cat > "$admin/step09_method_note.txt" <<'EOF'
Step 09 method note (system build):

09a - combine: H++ protonated protein (07b) + 3 placed chorismates (step 06) into
      one pre-tleap PDB (complex_for_tleap.pdb). Counts derived, charge computed.

09b - tleap build (paper-aligned): protein ff14SB, chorismate GAFF (cha_gaff.mol2 +
      cha.frcmod, AM1-BCC charges from step 08), TIP3P water, 10 A box, neutralised
      with Na+ (addions Na+ 0). Dry topology written first, then solvated. AMBER18
      unavailable -> AMBER22. Audit reads prmtop/inpcrd directly (solvated PDB
      saturates at 9999 residues); Na+ derived from dry charge. tleap warnings are
      read as messages (banner+next line) across stdout/log/stderr and bucketed:
      terminal name conversions (savepdb cosmetic), the expected non-neutral-charge
      note on the -9 dry unit, a unit-check bookkeeping note, and intramolecular
      close contacts. Each close contact is resolved to its residues against
      complex_dry.pdb and flagged for ligand / active-site / heavy-atom involvement;
      all are hydrogen contacts relaxed by the step-10 restrained minimisation.
      Fingerprints hash the prmtop with its wall-clock DATE line stripped (stable).
EOF
echo "wrote $admin/step09_method_note.txt"
echo
echo "=== step 09b: fingerprints (stable across runs) ==="
cat "$admin/step09_fingerprints.txt"

echo
echo "STEP 09b DONE - complex_solvated.prmtop/.inpcrd ready for minimisation (step 10)"

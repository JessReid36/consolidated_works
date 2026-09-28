#!/usr/bin/env bash
# Step 07b - audit, repair, and ACCEPT the H++ server output for tleap.
#
# One pass replacing attempt_3's four scripts (15 audit, 15b chain-restore+rename,
# 15c drift-classify, 15d corrected-drift+accept). Two things carried over
# deliberately from 15d (the corrected version):
#   * drift is classified on THREE independent axes by EXACT set membership -
#     backbone/side-chain, flippable/non-flippable (ASN/GLN/HIS reorient under
#     H++/Reduce), active-site/not - so a benign flip is never mistaken for
#     scaffold damage, and active-site integrity is asserted directly. (15c had a
#     substring-matching bug that mislabelled every atom active-site; fixed here.)
#   * an explicit ACCEPT / DO_NOT_ACCEPT decision with reasons is written.
# Expected counts are derived from the step-07a template, not hard-coded.
#
# Active-site set = published Agbaglo QM-cluster residues (matches step 04):
#   {7,57,59,60,63,73,74,75,78,90,108,115}.
#
# H++ run settings for the record: pH 7.0, salinity 0.15 M, internal dielectric
# 10, external dielectric 80.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
       NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

root="$HOME/system_development"
hppdir="$root/03_amber/protonation_hpp"
admin="$root/00_admin"
template="$hppdir/hpp_input_protein_only_with_TER.pdb"   # 07a output: chains + 1-127 numbering
hppout="$hppdir/hpp_output_from_server.pdb"              # raw H++ return
accepted="$hppdir/abc_protonated_hpp_accepted.pdb"       # <- written for tleap
mkdir -p "$admin"

python3 - "$template" "$hppout" "$accepted" "$admin" <<'PY'
import sys, math
from collections import defaultdict, Counter
from pathlib import Path

template, hppout, accepted, admin = (Path(a) for a in sys.argv[1:5])
his_tsv     = admin / "step07b_his_assignment.tsv"
drift_atom  = admin / "step07b_drift_by_atom.tsv"
drift_res   = admin / "step07b_drift_by_residue.tsv"
summ_tsv    = admin / "step07b_audit_summary.tsv"
decision_txt= admin / "step07b_acceptance_decision.txt"

BACKBONE     = {"N", "CA", "C", "O", "OXT"}
ACTIVE_SITE  = {7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115}   # Agbaglo QM set (step 04)
FLIPPABLE    = {"ASN", "GLN", "HIS", "HID", "HIE", "HIP"}
HIS_RING     = {"ND1", "NE2", "CE1", "CD2", "CG"}
FF14SB_STD   = {
    "ALA","ARG","ASN","ASP","ASH","CYS","CYX","CYM","GLN","GLU","GLH","GLY",
    "HID","HIE","HIP","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER","THR",
    "TRP","TYR","VAL",
}
# acceptance thresholds (from attempt_3 step 15d)
BB_MAX, NONFLIP_MAX, ACTIVE_BB_MAX, ACTIVE_BIG = 0.05, 0.20, 0.05, 1.0

def die(msg): sys.exit(f"FAIL {msg}")
for p in (template, hppout):
    if not p.is_file(): die(f"missing input: {p}")

# ---- fixed-column PDB parsing ------------------------------------------------
def is_atom(l): return l.startswith("ATOM  ") or l.startswith("HETATM")
def rec(l):     return l[0:6]
def aname(l):   return l[12:16].strip()
def resname(l): return l[17:20].strip()
def chain(l):   return l[21]
def resseq(l):  return l[22:26].strip()
def icode(l):   return l[26]
def xyz(l):     return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
def element(l): return l[76:78].strip()
def is_h(l):
    e = element(l)
    if e: return e == "H"
    return aname(l).lstrip("0123456789")[:1] == "H"
def dist(p, q): return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2)

def read_residues(path, use_chain):
    residues, cur, cur_key, ter = [], None, None, 0
    for l in path.read_text().splitlines():
        if l.startswith("TER"): ter += 1; cur = None; cur_key = None; continue
        if not is_atom(l): continue
        key = (chain(l), resseq(l), icode(l)) if use_chain else (resseq(l), icode(l), resname(l))
        if key != cur_key:
            cur = {"chain": chain(l), "resseq": resseq(l), "resname": resname(l), "lines": []}
            residues.append(cur); cur_key = key
        cur["lines"].append(l)
    return residues, ter

tmpl_res, _      = read_residues(template, use_chain=True)
hpp_res, hpp_ter = read_residues(hppout,  use_chain=False)
if not tmpl_res: die("template has no residues")
if not hpp_res:  die("H++ output has no residues")
if len(tmpl_res) != len(hpp_res):
    die(f"residue-count mismatch: template {len(tmpl_res)} vs H++ {len(hpp_res)}")

def is_histidine(reslines): return HIS_RING.issubset({aname(l) for l in reslines})

# ---- lockstep: restore chain+number, rename His -----------------------------
restored, his_rows = [], []
serial = 0
for t, h in zip(tmpl_res, hpp_res):
    ch = t["chain"]; rs = int(t["resseq"])
    rn = h["resname"] or t["resname"]
    if is_histidine(h["lines"]) or rn in {"HIS","HID","HIE","HIP","HDP"}:
        names = {aname(l) for l in h["lines"]}
        hd1, he2 = "HD1" in names, "HE2" in names
        if rn in {"HID","HIE","HIP"}: suggested = rn
        elif hd1 and he2: suggested = "HIP"
        elif hd1:         suggested = "HID"
        elif he2:         suggested = "HIE"
        else:             suggested = "HIE"     # neither placed: default, flagged
        his_rows.append((ch, rs, h["resname"], "yes" if hd1 else "no",
                         "yes" if he2 else "no", suggested))
        rn = suggested
    for l in h["lines"]:
        serial += 1
        restored.append((ch, rs, rn, t["resname"],
                         f"{rec(l)}{serial:5d} {l[12:16]}{l[16]}{rn:>3s} {ch}{rs:4d}{l[26]}{l[27:]}"))

# ---- accepted PDB (TER between chains) --------------------------------------
with accepted.open("w") as f:
    prev = None
    for ch, rs, rn, orn, line in restored:
        if prev is not None and ch != prev: f.write("TER\n")
        f.write(line + "\n"); prev = ch
    f.write("TER\nEND\n")

# ---- indexes ----------------------------------------------------------------
tmpl_heavy = {}
tmpl_resname = {}
for t in tmpl_res:
    for l in t["lines"]:
        if not is_h(l):
            tmpl_heavy[(t["chain"], int(t["resseq"]), aname(l))] = xyz(l)
    tmpl_resname[(t["chain"], int(t["resseq"]))] = t["resname"]

restored_heavy, atoms_by_res, xyz_by_res = {}, defaultdict(list), defaultdict(dict)
resname_by_res, res_ids = {}, defaultdict(set)
for ch, rs, rn, orn, line in restored:
    a = aname(line)
    atoms_by_res[(ch, rs)].append(a); resname_by_res[(ch, rs)] = rn; res_ids[ch].add(rs)
    if not is_h(line):
        restored_heavy[(ch, rs, a)] = xyz(line); xyz_by_res[(ch, rs)][a] = xyz(line)

# ---- chain ranges (expected from template) ----------------------------------
tmpl_span, tmpl_nres = {}, Counter(t["chain"] for t in tmpl_res)
for t in tmpl_res:
    ch = t["chain"]; r = int(t["resseq"])
    lo, hi = tmpl_span.get(ch, (10**9, -1)); tmpl_span[ch] = (min(lo, r), max(hi, r))
chain_rows, chain_fail = [], []
for ch in sorted(tmpl_span):
    ids = sorted(res_ids.get(ch, []))
    if not ids: chain_rows.append((ch,"MISSING","MISSING",0,"FAIL")); chain_fail.append(f"chain {ch} missing"); continue
    lo, hi, n = min(ids), max(ids), len(ids)
    elo, ehi = tmpl_span[ch]; en = tmpl_nres[ch]
    ok = (lo==elo and hi==ehi and n==en)
    chain_rows.append((ch, lo, hi, n, "PASS" if ok else "FAIL"))
    if not ok: chain_fail.append(f"chain {ch}: {lo}-{hi} n={n} (expected {elo}-{ehi}, {en})")

# ---- three-axis drift classification (exact membership; from step 15d) ------
common = sorted(set(tmpl_heavy) & set(restored_heavy))
missing_heavy = sorted(set(tmpl_heavy) - set(restored_heavy))
extra_heavy   = sorted(set(restored_heavy) - set(tmpl_heavy))
rows = []
for (ch, rs, at) in common:
    d = dist(tmpl_heavy[(ch,rs,at)], restored_heavy[(ch,rs,at)])
    orn = tmpl_resname.get((ch,rs), "?"); nrn = resname_by_res.get((ch,rs), "?")
    rows.append({
        "chain": ch, "resid": rs, "atom": at, "orig": orn, "hpp": nrn, "drift": d,
        "atom_class": "backbone" if at in BACKBONE else "sidechain",
        "flip_class": "flippable_residue" if (orn in FLIPPABLE or nrn in FLIPPABLE) else "nonflippable_residue",
        "site_class": "active_site_set" if rs in ACTIVE_SITE else "not_active_site_set",
    })

def stat(label, filt):
    vals = [r["drift"] for r in rows if filt(r)]
    if not vals: return (label, 0, 0.0, 0.0, 0, 0, 0)
    return (label, len(vals), max(vals), math.sqrt(sum(v*v for v in vals)/len(vals)),
            sum(v>0.1 for v in vals), sum(v>0.5 for v in vals), sum(v>1.0 for v in vals))

stat_rows = [
    stat("all_heavy",              lambda r: True),
    stat("backbone",               lambda r: r["atom_class"]=="backbone"),
    stat("sidechain",              lambda r: r["atom_class"]=="sidechain"),
    stat("flippable_residue",      lambda r: r["flip_class"]=="flippable_residue"),
    stat("nonflippable_residue",   lambda r: r["flip_class"]=="nonflippable_residue"),
    stat("active_site_set",        lambda r: r["site_class"]=="active_site_set"),
    stat("not_active_site_set",    lambda r: r["site_class"]=="not_active_site_set"),
    stat("active_site_backbone",   lambda r: r["site_class"]=="active_site_set" and r["atom_class"]=="backbone"),
    stat("active_site_sidechain",  lambda r: r["site_class"]=="active_site_set" and r["atom_class"]=="sidechain"),
]
S = {s[0]: s for s in stat_rows}   # index: [1]=count [2]=max [3]=rms [4]=>0.1 [5]=>0.5 [6]=>1.0

# "bad" large shifts: a >1.0 A move that is non-flippable OR active-site (i.e. NOT a benign flip)
large     = [r for r in rows if r["drift"] > ACTIVE_BIG]
large_bad = [r for r in large if r["flip_class"] != "flippable_residue" or r["site_class"] == "active_site_set"]

# ---- other structural audits ------------------------------------------------
bad_links, checked = [], 0
for ch in sorted(res_ids):
    ids = res_ids[ch]
    for r1 in sorted(ids):
        if r1+1 not in ids: continue
        C = xyz_by_res[(ch,r1)].get("C"); N = xyz_by_res[(ch,r1+1)].get("N")
        if C is None or N is None: bad_links.append((ch,r1,r1+1,"missing_C_or_N",None)); continue
        checked += 1; d = dist(C,N)
        if d < 1.15 or d > 1.70: bad_links.append((ch,r1,r1+1,"abnormal_C_N",d))
missing_CN  = [b for b in bad_links if b[3]=="missing_C_or_N"]
abnormal_CN = [b for b in bad_links if b[3]=="abnormal_C_N"]

dup_rows = [(k, [n for n,c in Counter(v).items() if c>1]) for k,v in atoms_by_res.items()
            if any(c>1 for c in Counter(v).values())]

sg = [(ch,rs,xyz_by_res[(ch,rs)]["SG"]) for (ch,rs) in xyz_by_res
      if "SG" in xyz_by_res[(ch,rs)] and resname_by_res.get((ch,rs)) in {"CYS","CYX","CYM"}]
ss_pairs = [(sg[i][0],sg[i][1],sg[j][0],sg[j][1],dist(sg[i][2],sg[j][2]))
            for i in range(len(sg)) for j in range(i+1,len(sg)) if dist(sg[i][2],sg[j][2]) < 2.5]

restored_h_n = sum(1 for ch,rs,rn,orn,line in restored if is_h(line))
resname_counts = Counter(rn for (ch,rs),rn in resname_by_res.items())
non_std = sorted(n for n in resname_counts if n not in FF14SB_STD)
his_default = [r for r in his_rows if r[2] not in {"HID","HIE","HIP"} and r[3]=="no" and r[4]=="no"]

# ---- ACCEPTANCE RULE (step 15d) ---------------------------------------------
hard = []
if chain_fail:      hard.append("chain range/identity failure")
if restored_h_n==0: hard.append("no hydrogens present")
if dup_rows:        hard.append("duplicate atom names in a residue")
if missing_CN:      hard.append(f"{len(missing_CN)} peptide link(s) missing C/N")
if missing_heavy:   hard.append(f"{len(missing_heavy)} original heavy atom(s) missing")

accept = (not hard
          and len(missing_heavy)==0 and len(extra_heavy)==0
          and S["backbone"][2]           <= BB_MAX
          and S["nonflippable_residue"][2] <= NONFLIP_MAX
          and S["active_site_backbone"][2] <= ACTIVE_BB_MAX
          and S["active_site_set"][6]    == 0
          and len(large_bad)==0)

reasons = []
if len(missing_heavy)==0 and len(extra_heavy)==0: reasons.append("all template heavy atoms preserved, none added")
if S["backbone"][2] <= BB_MAX:                    reasons.append(f"backbone drift {S['backbone'][2]:.3f} A <= {BB_MAX}")
if S["nonflippable_residue"][2] <= NONFLIP_MAX:   reasons.append(f"nonflippable drift {S['nonflippable_residue'][2]:.3f} A <= {NONFLIP_MAX}")
if S["active_site_backbone"][2] <= ACTIVE_BB_MAX: reasons.append(f"active-site backbone drift {S['active_site_backbone'][2]:.3f} A <= {ACTIVE_BB_MAX}")
if S["active_site_set"][6] == 0:                  reasons.append("no active-site heavy atom moved > 1.0 A")
if large and not large_bad:                       reasons.append("large shifts confined to non-active-site flippable ASN/GLN/HIS residues")

# ---- reports ----------------------------------------------------------------
with his_tsv.open("w") as f:
    f.write("chain\tresid\tinput_resname\tHD1\tHE2\tsuggested\n")
    for r in his_rows: f.write("\t".join(map(str,r))+"\n")

with drift_atom.open("w") as f:
    f.write("chain\tresid\tatom\torig_resname\thpp_resname\tdrift_A\tatom_class\tflip_class\tsite_class\n")
    for r in sorted(rows, key=lambda x:x["drift"], reverse=True):
        f.write(f"{r['chain']}\t{r['resid']}\t{r['atom']}\t{r['orig']}\t{r['hpp']}\t"
                f"{r['drift']:.6f}\t{r['atom_class']}\t{r['flip_class']}\t{r['site_class']}\n")

by_res = defaultdict(list)
for r in rows: by_res[(r["chain"],r["resid"],r["orig"],r["hpp"])].append(r)
with drift_res.open("w") as f:
    f.write("chain\tresid\torig_resname\thpp_resname\theavy\tmax_drift_A\trms_drift_A\tgt1.0\tis_flippable\tis_active_site\n")
    recs = []
    for (ch,rs,orn,nrn), rr in by_res.items():
        recs.append((ch,rs,orn,nrn,len(rr),max(x["drift"] for x in rr),
                     math.sqrt(sum(x["drift"]**2 for x in rr)/len(rr)),
                     sum(x["drift"]>1.0 for x in rr),
                     any(x["flip_class"]=="flippable_residue" for x in rr),
                     any(x["site_class"]=="active_site_set" for x in rr)))
    for ch,rs,orn,nrn,nh,mx,rm,g1,fl,asf in sorted(recs, key=lambda x:x[5], reverse=True):
        f.write(f"{ch}\t{rs}\t{orn}\t{nrn}\t{nh}\t{mx:.6f}\t{rm:.6f}\t{g1}\t{fl}\t{asf}\n")

with summ_tsv.open("w") as f:
    f.write("metric\tvalue\n")
    for k,v in [("template_residues",len(tmpl_res)),("hpp_residues",len(hpp_res)),
                ("template_heavy_atoms",len(tmpl_heavy)),("restored_heavy_atoms",len(restored_heavy)),
                ("hydrogens_added",restored_h_n),("hpp_TER_count",hpp_ter),
                ("histidines",len(his_rows)),
                ("suggested_his",",".join(f"{a}:{b}" for a,b in sorted(Counter(r[5] for r in his_rows).items())) or "none"),
                ("resname_counts",",".join(f"{a}:{b}" for a,b in sorted(resname_counts.items()))),
                ("non_ff14sb_resnames",",".join(non_std) or "none"),
                ("missing_heavy_atoms",len(missing_heavy)),("extra_heavy_atoms",len(extra_heavy)),
                ("peptide_links_checked",checked),("missing_CN_links",len(missing_CN)),
                ("abnormal_CN_links",len(abnormal_CN)),("duplicate_atom_residues",len(dup_rows)),
                ("possible_disulfides",len(ss_pairs))]:
        f.write(f"{k}\t{v}\n")
    for label,count,mx,rm,g1,g5,g10 in stat_rows:
        f.write(f"drift[{label}]\tcount={count};max={mx:.4f};rms={rm:.4f};>1.0A={g10}\n")

with decision_txt.open("w") as f:
    f.write("STEP 07b - H++ protonated protein acceptance decision\n\n")
    f.write(("ACCEPT" if accept else "DO_NOT_ACCEPT") + " for tleap after HIS renaming\n\n")
    if accept:
        f.write("Reasons:\n")
        for r in reasons: f.write(f"  - {r}\n")
        f.write(f"\nAccepted structure:\n  {accepted}\n")
        f.write("\nNext: combine with the accepted CHA ligands and build with tleap (ff14SB + GAFF + TIP3P).\n")
        f.write("Note: H++'s own topology/coordinates are NOT used - only its protonation/orientation decisions.\n")
    else:
        f.write("Blocking issues:\n")
        for x in hard: f.write(f"  - {x}\n")
        if len(extra_heavy): f.write(f"  - {len(extra_heavy)} extra heavy atom(s)\n")
        if S["backbone"][2] > BB_MAX: f.write(f"  - backbone drift {S['backbone'][2]:.3f} A > {BB_MAX}\n")
        if S["nonflippable_residue"][2] > NONFLIP_MAX: f.write(f"  - nonflippable drift {S['nonflippable_residue'][2]:.3f} A > {NONFLIP_MAX}\n")
        if S["active_site_backbone"][2] > ACTIVE_BB_MAX: f.write(f"  - active-site backbone drift {S['active_site_backbone'][2]:.3f} A > {ACTIVE_BB_MAX}\n")
        if S["active_site_set"][6] > 0: f.write(f"  - {S['active_site_set'][6]} active-site heavy atom(s) moved > 1.0 A\n")
        if large_bad: f.write(f"  - {len(large_bad)} large shift(s) are non-flippable or active-site\n")

# ---- stdout -----------------------------------------------------------------
print("STEP 07b - H++ output audit + repair + acceptance")
print(f"  residues: template {len(tmpl_res)}  H++ {len(hpp_res)}  (matched)")
print(f"  heavy atoms: template {len(tmpl_heavy)}  restored {len(restored_heavy)}  + {restored_h_n} H  (TER={hpp_ter})")
for ch,lo,hi,n,st in chain_rows: print(f"  chain {ch}: {lo}-{hi} ({n} res) {st}")
print(f"  histidines ({len(his_rows)}): " + ", ".join(f"{c}/{r}->{s}" for c,r,_,_,_,s in his_rows))
print("  drift by class (max A | rms | >1.0A):")
for label in ("backbone","active_site_backbone","active_site_set","nonflippable_residue","sidechain","flippable_residue"):
    s = S[label]; print(f"    {label:22s} {s[2]:7.3f} | {s[3]:6.3f} | {s[6]}")
if large:
    print(f"  large (>1.0 A) shifts: {len(large)} total, {len(large_bad)} 'bad' (non-flip or active-site)")
    for r in sorted(large, key=lambda x:x['drift'], reverse=True)[:5]:
        print(f"    {r['drift']:.3f} A  {r['chain']}/{r['resid']} {r['hpp']} {r['atom']}  ({r['flip_class']}, {r['site_class']})")
print(f"  peptide links: {checked} checked, {len(missing_CN)} missing C/N, {len(abnormal_CN)} abnormal")
print(f"  cysteines: {resname_counts.get('CYS',0)} CYS / {resname_counts.get('CYX',0)} CYX; SG<2.5A pairs: {len(ss_pairs)}")
if non_std: print(f"  non-ff14SB resnames: {','.join(non_std)}")
if his_default: print(f"  ATTENTION {len(his_default)} histidine(s) had neither HD1 nor HE2 (defaulted HIE)")

print()
if accept:
    print("  RESULT: ACCEPT - protonated protein ready for tleap (step 08/09)")
    for r in reasons: print(f"    + {r}")
else:
    print("  RESULT: DO_NOT_ACCEPT")
    for x in hard: print(f"    FAIL {x}")

print(f"\n  WROTE {accepted}")
for p in (his_tsv, drift_atom, drift_res, summ_tsv, decision_txt): print(f"  WROTE {p}")

if not accept:
    sys.exit("\nRESULT: DO_NOT_ACCEPT - resolve issues before tleap")
PY
echo "STEP 07b DONE"

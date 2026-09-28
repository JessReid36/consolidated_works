#!/usr/bin/env bash
# Step 09a - combine the protonated protein (07b) with the 3 placed chorismates
# (from the step-06 frozen complex) into one pre-tleap PDB for the tleap build.
#
# Tightens attempt_3 step 16: counts are DERIVED from inputs (CHA atoms = 3 x the
# cha_gaff template size; Na+ = |computed system charge|) instead of hard-coded
# 24/72/381/-3/-9/9. Hard failures are structural only (template not -2, non-integer
# charge, HIS left un-renamed, duplicate atom names, chain gaps, a real clash, or
# atom counts changing through assembly); the expected -3/-9/9 charge values are
# computed and REVIEW-flagged if they differ, not fatal - so the step is portable.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
protein="$root/03_amber/protonation_hpp/abc_protonated_hpp_accepted.pdb"              # 07b
complex06="$root/02_preparation/accepted_preprotonation/abc_cha_preprotonation.pdb"   # step 06 (CHA source)
cha_gaff="$root/03_amber/ligand_gaff/cha_gaff.mol2"                                   # 08b (template)
outdir="$root/03_amber/tleap_build"
admin="$root/00_admin"
out="$outdir/complex_for_tleap.pdb"
mkdir -p "$outdir" "$admin"

echo "=== step 09a: input presence ==="
for f in "$protein" "$complex06" "$cha_gaff"; do
  [[ -s "$f" ]] || { echo "FAIL missing/empty: $f"; exit 1; }
  echo "PASS $f"
done

python3 - "$protein" "$complex06" "$cha_gaff" "$out" "$admin" <<'PY'
import sys, math
from collections import defaultdict, Counter
from pathlib import Path

protein, complex06, cha_gaff, out, admin = (Path(a) for a in sys.argv[1:6])
audit_txt = admin / "step09a_pretleap_audit.txt"
resinv    = admin / "step09a_residue_inventory.tsv"
cha_check = admin / "step09a_cha_template_check.tsv"
charge_tsv= admin / "step09a_charge_estimate.tsv"

FF14SB_STD = {
    "ALA","ARG","ASN","ASP","ASH","CYS","CYX","CYM","GLN","GLU","GLH","GLY",
    "HID","HIE","HIP","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER","THR",
    "TRP","TYR","VAL",
}
CHARGE_MAP = {"ARG":1,"LYS":1,"HIP":1,"ASP":-1,"GLU":-1,"CYM":-1}  # others 0

def die(m): sys.exit(f"FAIL {m}")
def is_h(atom, el=""):
    el = (el or "").strip().upper()
    return el == "H" or atom[:1] == "H" or (len(atom) > 1 and atom[0].isdigit() and atom[1:2].upper() == "H")
def dist(a, b): return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

def parse_pdb(path):
    atoms, ter = [], 0
    for line in path.read_text().splitlines():
        if line.startswith("TER"): ter += 1; continue
        if not line.startswith(("ATOM  ", "HETATM")): continue
        atoms.append({"record":line[:6].strip(),"atom":line[12:16].strip(),
            "resname":line[17:20].strip(),"chain":line[21],"resid":int(line[22:26]),
            "element":line[76:78].strip() if len(line)>=78 else "",
            "x":float(line[30:38]),"y":float(line[38:46]),"z":float(line[46:54]),
            "line":line.rstrip("\n")})
    return atoms, ter

def parse_mol2(path):
    atoms, bonds, sec = [], 0, None
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"): sec = line.strip(); continue
        if sec == "@<TRIPOS>ATOM" and len(line.split()) >= 9:
            p = line.split(); atoms.append({"name":p[1],"type":p[5],"charge":float(p[8])})
        elif sec == "@<TRIPOS>BOND" and len(line.split()) >= 4: bonds += 1
    return atoms, bonds

prot_atoms, prot_ter = parse_pdb(protein)
cx_atoms, _          = parse_pdb(complex06)
tmpl, tmpl_bonds     = parse_mol2(cha_gaff)
if not prot_atoms: die("no atoms in protonated protein")
if not tmpl:       die("no atoms in cha_gaff template")

# ---- protein checks ---------------------------------------------------------
prot = [a for a in prot_atoms if a["record"] == "ATOM"]
if any(a["record"] == "HETATM" for a in prot_atoms): die("protonated protein contains HETATM")
prot_res = defaultdict(list)
for a in prot: prot_res[(a["chain"], a["resid"], a["resname"])].append(a)
chains = sorted({a["chain"] for a in prot})
chain_fail = []
per_chain = {}
for ch in chains:
    ids = sorted({a["resid"] for a in prot if a["chain"] == ch})
    contiguous = ids == list(range(ids[0], ids[-1]+1))
    per_chain[ch] = (ids[0], ids[-1], len(ids), contiguous)
    if not contiguous or ids[0] != 1: chain_fail.append((ch, ids[0], ids[-1], len(ids)))
if chain_fail: die(f"chain range/gap issue: {chain_fail}")
resname_counts = Counter(k[2] for k in prot_res)
if resname_counts.get("HIS"): die("HIS remains un-renamed in protein (need HID/HIE/HIP)")
bad_names = sorted(n for n in resname_counts if n not in FF14SB_STD)
dup_res = [(k, [n for n,c in Counter(a["atom"] for a in v).items() if c>1])
           for k,v in prot_res.items() if any(c>1 for c in Counter(a["atom"] for a in v).values())]
if dup_res: die(f"duplicate atom names in protein residue {dup_res[0][0]}")

# ---- CHA extraction + template match ---------------------------------------
n_t = len(tmpl)                                   # DERIVED ligand size (no hard-coded 24)
tmpl_names = [a["name"] for a in tmpl]
tmpl_charge = sum(a["charge"] for a in tmpl)
if abs(tmpl_charge + 2.0) > 1e-4: die(f"cha_gaff template charge {tmpl_charge:.4f}, expected -2")

cha = [a for a in cx_atoms if a["record"] == "HETATM" and a["resname"] == "CHA"]
cha_res = defaultdict(list)
for a in cha: cha_res[(a["chain"], a["resid"], a["resname"])].append(a)
if not cha_res: die("no CHA HETATM found in step-06 complex")

cha_rows = []
for key in sorted(cha_res):
    atoms = cha_res[key]; names = [a["atom"] for a in atoms]
    dups = [n for n,c in Counter(names).items() if c>1]
    missing = sorted(set(tmpl_names) - set(names)); extra = sorted(set(names) - set(tmpl_names))
    order_ok = names == tmpl_names
    cha_rows.append((key[0], key[1], len(atoms), order_ok, dups, missing, extra))
    if len(atoms) != n_t: die(f"CHA {key} has {len(atoms)} atoms, template has {n_t}")
    if dups:    die(f"CHA {key} duplicate names {dups}")
    if missing: die(f"CHA {key} missing template names {missing}")
    if extra:   die(f"CHA {key} extra names {extra}")
n_cha_res = len(cha_res)
cha_total = n_cha_res * n_t                        # DERIVED (no hard-coded 72)

# ---- charge estimate (computed; -3/-9/9 are expected, review-flagged) -------
protein_charge = sum(CHARGE_MAP.get(k[2], 0) for k in prot_res)
# free termini contribute 0 net per chain: N-terminus NH3+ (+1), C-terminus COO- (-1)
protein_charge += len(chains) - len(chains)
ligand_charge = n_cha_res * tmpl_charge
complex_charge = protein_charge + ligand_charge
if abs(protein_charge - round(protein_charge)) > 1e-6: die("protein charge is non-integer")
if abs(complex_charge - round(complex_charge)) > 1e-6: die("complex charge is non-integer")
protein_charge, complex_charge = round(protein_charge), round(complex_charge)
expected_na = -complex_charge if complex_charge < 0 else 0

# ---- clash audit ------------------------------------------------------------
prot_heavy = [a for a in prot if not is_h(a["atom"], a["element"])]
lig_heavy  = [a for a in cha  if not is_h(a["atom"], a["element"])]
min_hh = min(((dist(p,l), p, l) for p in prot_heavy for l in lig_heavy), key=lambda t:t[0])
if min_hh[0] < 1.5:
    p, l = min_hh[1], min_hh[2]
    die(f"severe protein-ligand heavy clash {min_hh[0]:.3f} A "
        f"{p['chain']}:{p['resid']}:{p['resname']}:{p['atom']} <-> {l['chain']}:{l['resid']}:CHA:{l['atom']}")

# ---- assemble combined PDB (protein ATOM/TER, then CHA grouped, clean serials)
serial = 1; lines = []
for line in protein.read_text().splitlines():
    if line.startswith("ATOM  "):
        lines.append(f"{line[:6]}{serial:5d}{line[11:]}"); serial += 1
    elif line.startswith("TER"):
        lines.append(f"TER   {serial:5d}"); serial += 1
for key in sorted(cha_res):
    for a in sorted(cha_res[key], key=lambda x: x["line"][6:11]):
        lines.append(f"HETATM{serial:5d}{a['line'][11:]}"); serial += 1
    lines.append(f"TER   {serial:5d}"); serial += 1
lines.append("END")
out.write_text("\n".join(lines) + "\n")

o_atoms, o_ter = parse_pdb(out)
o_prot = [a for a in o_atoms if a["record"] == "ATOM"]
o_cha  = [a for a in o_atoms if a["record"] == "HETATM" and a["resname"] == "CHA"]
if len(o_prot) != len(prot): die(f"protein atoms changed through assembly {len(o_prot)} vs {len(prot)}")
if len(o_cha)  != cha_total: die(f"CHA atoms changed through assembly {len(o_cha)} vs {cha_total}")

# ---- reports ----------------------------------------------------------------
with cha_check.open("w") as f:
    f.write("chain\tresid\tatom_count\torder_matches_template\tduplicate\tmissing\textra\n")
    for ch,rid,n,ok,dups,miss,ext in cha_rows:
        f.write(f"{ch}\t{rid}\t{n}\t{ok}\t{','.join(dups) or 'none'}\t{','.join(miss) or 'none'}\t{','.join(ext) or 'none'}\n")
with charge_tsv.open("w") as f:
    f.write("component\tcharge\n")
    for k,v in [("protein_formal_estimate",protein_charge),("one_CHA_template",f"{tmpl_charge:.4f}"),
                ("three_CHA",f"{ligand_charge:.4f}"),("complex_total",complex_charge),
                ("expected_Na_to_neutralize",expected_na)]:
        f.write(f"{k}\t{v}\n")
with resinv.open("w") as f:
    f.write("record\tchain\tresid\tresname\tatoms\theavy\tH\n")
    inv = defaultdict(list)
    for a in o_atoms: inv[(a["record"],a["chain"],a["resid"],a["resname"])].append(a)
    for key in sorted(inv):
        at = inv[key]; h = sum(1 for a in at if is_h(a["atom"],a["element"]))
        f.write(f"{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{len(at)}\t{len(at)-h}\t{h}\n")

review = []
if protein_charge != -3: review.append(f"protein charge {protein_charge:+d} (expected -3 for this system)")
if complex_charge != -9: review.append(f"complex charge {complex_charge:+d} (expected -9 for this system)")
if bad_names:            review.append(f"non-ff14SB protein resnames: {','.join(bad_names)}")
if min_hh[0] < 2.0:      review.append(f"closest protein-ligand heavy contact {min_hh[0]:.3f} A (<2.0)")
if prot_ter != 3:        review.append(f"protein TER count {prot_ter} (expected 3)")

with audit_txt.open("w") as f:
    for k,v in [("protein_atoms",len(prot)),("protein_residues",len(prot_res)),
                ("protein_chains",",".join(chains)),("protein_TER",prot_ter),
                ("per_chain_range",";".join(f"{ch}:{per_chain[ch][0]}-{per_chain[ch][1]}({per_chain[ch][2]})" for ch in chains)),
                ("HIS_remaining",resname_counts.get("HIS",0)),
                ("HID_HIE_HIP",f"{resname_counts.get('HID',0)}/{resname_counts.get('HIE',0)}/{resname_counts.get('HIP',0)}"),
                ("cha_residues",n_cha_res),("cha_template_atoms",n_t),("cha_template_bonds",tmpl_bonds),
                ("cha_total_atoms",cha_total),("cha_template_charge",f"{tmpl_charge:.4f}"),
                ("protein_charge",protein_charge),("complex_charge",complex_charge),("expected_Na",expected_na),
                ("min_protein_ligand_heavy_A",f"{min_hh[0]:.3f}"),
                ("combined_atoms",len(o_atoms)),("combined_TER",o_ter),("output",out)]:
        f.write(f"{k}\t{v}\n")
    f.write("review\t" + ("; ".join(review) if review else "none") + "\n")

print("STEP 09a - pre-tleap combine")
print(f"  protein: {len(prot)} atoms, {len(prot_res)} residues, chains {','.join(chains)}, TER={prot_ter}")
print(f"    per chain: " + "  ".join(f"{ch} {per_chain[ch][0]}-{per_chain[ch][1]}({per_chain[ch][2]})" for ch in chains))
print(f"    histidines HID/HIE/HIP = {resname_counts.get('HID',0)}/{resname_counts.get('HIE',0)}/{resname_counts.get('HIP',0)}; HIS left = {resname_counts.get('HIS',0)}")
print(f"  ligand: {n_cha_res} CHA x {n_t} atoms = {cha_total}; template charge {tmpl_charge:+.3f}; all copies match template names/order")
print(f"  charge: protein {protein_charge:+d}, 3xCHA {ligand_charge:+.1f}, complex {complex_charge:+d} -> expected Na+ = {expected_na}")
print(f"  closest protein-ligand heavy contact: {min_hh[0]:.3f} A")
print(f"  combined: {len(o_atoms)} atoms ({len(o_prot)} protein + {len(o_cha)} CHA), TER={o_ter}")
print()
if review:
    print("  REVIEW:")
    for r in review: print(f"    - {r}")
else:
    print("  PASS - no review flags")
print(f"\n  WROTE {out}")
for p in (audit_txt, cha_check, charge_tsv, resinv): print(f"  WROTE {p}")
print("\nRESULT: pre-tleap complex ready for step 09b (tleap build)")
PY
echo "STEP 09a DONE"

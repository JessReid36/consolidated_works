#!/usr/bin/env bash
# Step 03 - clean and validate the repaired A/B/C protein before ligand placement.
# Resolves alternate locations (keep highest occupancy), removes any hydrogens
# (element column is authoritative; a removal breakdown is printed so we can
# confirm only true H are dropped), and validates heavy-atom structure:
# no altlocs/H remain, backbone continuity intact, active-site residues present.
set -euo pipefail

root="$HOME/system_development"
prot_in="$root/02_preparation/protein_only/abc_repaired.pdb"
outdir="$root/02_preparation/protein_only"
admin="$root/00_admin"

python3 - "$prot_in" "$outdir" "$admin" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict, Counter

prot_in, outdir, admin = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
prot_out    = outdir / "abc_repaired_clean.pdb"
cleanup_rep = admin  / "step03_cleanup_report.tsv"
audit_rep   = admin  / "step03_validation_audit.txt"

def parse(line):
    return {
        "serial": int(line[6:11]), "name_field": line[12:16], "name": line[12:16].strip(),
        "altloc": line[16], "resname": line[17:20].strip(), "chain": line[21],
        "resid": int(line[22:26]), "icode": line[26],
        "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
        "occ": float(line[54:60]) if line[54:60].strip() else 1.0,
        "bfac": float(line[60:66]) if line[60:66].strip() else 0.0,
        "element": line[76:78].strip() if len(line) >= 78 else "",
    }

# Authoritative H detection: element column == H. Name-based fallback ONLY when
# the element column is empty, and reported separately so we can inspect it.
def h_by_element(a):
    return a["element"].upper() == "H"
def h_by_name_fallback(a):
    if a["element"]:              # element present -> trust it, no fallback
        return False
    n = a["name"]
    return n.startswith("H") or (len(n) > 1 and n[0].isdigit() and n[1].upper() == "H")

atoms = [parse(l) for l in prot_in.read_text().splitlines() if l.startswith("ATOM  ")]
assert atoms, "no ATOM records parsed"

n_initial = len(atoms)
n_elem_present = sum(1 for a in atoms if a["element"])
h_elem = [a for a in atoms if h_by_element(a)]
h_name = [a for a in atoms if h_by_name_fallback(a)]

print(f"initial atoms: {n_initial}")
print(f"atoms with element column populated: {n_elem_present}/{n_initial}")
print(f"hydrogens by ELEMENT column (H): {len(h_elem)}")
print(f"hydrogens by NAME fallback (element col empty): {len(h_name)}")

# Show what the NAME fallback would catch, grouped by element we'd infer -
# this is the diagnostic that reveals over-matching (heavy atoms named H*).
if h_name:
    by_resname_atom = Counter((a["resname"], a["name"]) for a in h_name)
    print("\nNAME-fallback matches (potential over-match - inspect):")
    for (rn, an), c in sorted(by_resname_atom.items()):
        print(f"  {rn} {an}: {c}")

# Decision: remove only element-column H. Do NOT trust name fallback for deletion
# unless element column is entirely absent (then we fall back, but report loudly).
if n_elem_present == 0:
    print("\nWARNING: element column empty throughout - using name fallback for H removal")
    is_h = h_by_name_fallback
else:
    is_h = h_by_element

# Resolve altlocs per atom identity (highest occupancy; tie -> altloc A, then blank).
def agroup(a): return (a["chain"], a["resid"], a["icode"], a["resname"], a["name"])
def altrank(a):
    pref = {"A": 0, " ": 1, "B": 2}.get(a["altloc"], 3)
    return (-a["occ"], pref, a["serial"])

groups = defaultdict(list)
for a in atoms:
    groups[agroup(a)].append(a)

selected, altloc_decisions = [], []
for key, recs in groups.items():
    if len(recs) == 1:
        selected.append(recs[0]); continue
    recs = sorted(recs, key=altrank)
    selected.append(recs[0])
    altloc_decisions.append((key, recs[0]["altloc"].strip() or "blank",
                             ",".join((r["altloc"].strip() or "blank") for r in recs[1:])))

heavy = [a for a in selected if not is_h(a)]
for a in heavy:
    a["altloc"] = " "

def infer_elem(a):
    if a["element"]: return a["element"]
    return "".join(c for c in a["name"] if c.isalpha())[:1].upper() or "X"

def sort_key(a):
    return ({"A":0,"B":1,"C":2}.get(a["chain"],9), a["resid"], a["icode"], a["serial"])

heavy.sort(key=sort_key)
with prot_out.open("w") as f:
    serial, chain = 1, None
    for a in heavy:
        if chain is not None and a["chain"] != chain:
            f.write("TER\n")
        chain = a["chain"]
        f.write(f"ATOM  {serial:5d} {a['name_field']}{a['altloc']}{a['resname']:>3s} "
                f"{a['chain']:1s}{a['resid']:4d}{a['icode']:1s}   "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                f"{a['occ']:6.2f}{a['bfac']:6.2f}          {infer_elem(a):>2s}\n")
        serial += 1
    f.write("TER\nEND\n")

# Re-audit the written file independently.
fa = [parse(l) for l in prot_out.read_text().splitlines() if l.startswith("ATOM  ")]
final_h   = [a for a in fa if h_by_element(a)]
final_alt = [a for a in fa if a["altloc"].strip()]
final_oxt = [a for a in fa if a["name"] == "OXT"]

standard = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
            "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}
nonstd = sorted(set(a["resname"] for a in fa if a["resname"] not in standard))

by_res = defaultdict(list)
for a in fa: by_res[(a["chain"], a["resid"])].append(a["name"])
dups = [(k, [n for n, c in Counter(v).items() if c > 1]) for k, v in by_res.items()
        if any(c > 1 for c in Counter(v).values())]

xyz = defaultdict(dict)
for a in fa: xyz[(a["chain"], a["resid"])][a["name"]] = (a["x"], a["y"], a["z"])
def d(p, q): return math.sqrt(sum((p[i]-q[i])**2 for i in range(3)))
bad_links, checked = [], 0
for ch in "ABC":
    for r in range(1, 127):
        C = xyz.get((ch, r), {}).get("C"); N = xyz.get((ch, r+1), {}).get("N")
        if C and N:
            checked += 1
            dd = d(C, N)
            if dd < 1.15 or dd > 1.70: bad_links.append((ch, r, r+1, dd))

required = [7,57,59,60,63,73,74,75,78,90,108,115]
seen = defaultdict(set)
for a in fa: seen[a["chain"]].add(a["resid"])
active_missing = [(ch, r) for ch in "ABC" for r in required if r not in seen[ch]]

counts = Counter(a["chain"] for a in fa)
res_per_chain = {ch: len(set(a["resid"] for a in fa if a["chain"] == ch)) for ch in "ABC"}

with cleanup_rep.open("w") as f:
    f.write("category\tchain\tresid\tresname\tatom\tkept_altloc\tdropped\n")
    for (ch, r, ic, rn, an), kept, dropped in altloc_decisions:
        f.write(f"altloc\t{ch}\t{r}\t{rn}\t{an}\t{kept}\t{dropped}\n")

with audit_rep.open("w") as f:
    f.write(f"input\t{prot_in}\noutput\t{prot_out}\n")
    f.write(f"initial_atoms\t{n_initial}\n")
    f.write(f"element_column_populated\t{n_elem_present}\n")
    f.write(f"hydrogens_removed_by_element\t{len(h_elem)}\n")
    f.write(f"altloc_groups_resolved\t{len(altloc_decisions)}\n")
    f.write(f"final_heavy_atoms\t{len(fa)}\n")
    f.write(f"final_hydrogens\t{len(final_h)}\n")
    f.write(f"final_altlocs\t{len(final_alt)}\n")
    f.write(f"final_OXT\t{len(final_oxt)}\n")
    f.write(f"nonstandard_resnames\t{','.join(nonstd) if nonstd else 'none'}\n")
    f.write(f"duplicate_atom_residues\t{len(dups)}\n")
    f.write(f"peptide_links_checked\t{checked}\n")
    f.write(f"bad_peptide_links\t{len(bad_links)}\n")
    f.write(f"active_site_missing\t{len(active_missing)}\n")

print(f"\naltloc groups resolved: {len(altloc_decisions)}")
print(f"final heavy atoms: {len(fa)}")
print(f"final hydrogens: {len(final_h)}  final altlocs: {len(final_alt)}  OXT: {len(final_oxt)}")
print(f"nonstandard resnames: {nonstd if nonstd else 'none'}")
print("\nchain summary:")
for ch in "ABC":
    print(f"  {ch}: {counts[ch]} atoms, {res_per_chain[ch]} residues")
print(f"\npeptide links checked: {checked}, bad: {len(bad_links)}")
for ch, r1, r2, dd in bad_links:
    print(f"  {ch} {r1}->{r2}: d={dd:.3f}")

# Hard gates.
if final_h:        sys.exit("FAIL hydrogens remain")
if final_alt:      sys.exit("FAIL altlocs remain")
if nonstd:         sys.exit(f"FAIL nonstandard residues: {nonstd}")
if dups:           sys.exit(f"FAIL duplicate atoms: {dups[:5]}")
if bad_links:      sys.exit("FAIL peptide continuity broken")
if active_missing: sys.exit(f"FAIL active-site residues missing: {active_missing}")
print("\nall validation gates passed")
print(f"WROTE {prot_out}")
PY
echo "STEP 03 DONE"

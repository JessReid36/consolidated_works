#!/usr/bin/env bash
# Step 06 - freeze the accepted pre-protonation complex.
#
# The CP2K-derived chorismate mol2 files use duplicate generic atom names
# (repeated C, H, O), which AMBER/tleap cannot handle. This step renames ligand
# atoms uniquely (C1, H1, O1, ...) while preserving coordinates, atom types,
# charges, and bond connectivity (net charge -2), rebuilds the combined complex
# with the safe names, audits for duplicates/charge, and freezes the accepted
# files with checksums.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
protein="$root/02_preparation/protein_only/abc_repaired_clean.pdb"
ligdir="$root/02_preparation/ligand_placement"
outdir="$root/02_preparation/accepted_preprotonation"
admin="$root/00_admin"
mkdir -p "$outdir" "$admin"

python3 - "$protein" "$ligdir" "$outdir" "$admin" <<'PY'
import sys
from pathlib import Path
from collections import Counter, defaultdict

protein, ligdir, outdir, admin = (Path(sys.argv[1]), Path(sys.argv[2]),
                                  Path(sys.argv[3]), Path(sys.argv[4]))

combined_out = outdir / "abc_cha_preprotonation.pdb"
audit_tsv    = admin  / "step06_ligand_mol2_audit.tsv"
pdb_audit    = admin  / "step06_combined_pdb_audit.txt"

lig_specs = [
    dict(label="cha_a", src=ligdir/"cha_a_placed.mol2", out=outdir/"cha_a.mol2",
         rename=outdir/"cha_a_renaming.tsv", chain="A", resid=201),
    dict(label="cha_b", src=ligdir/"cha_b_placed.mol2", out=outdir/"cha_b.mol2",
         rename=outdir/"cha_b_renaming.tsv", chain="B", resid=201),
    dict(label="cha_c", src=ligdir/"cha_c_placed.mol2", out=outdir/"cha_c.mol2",
         rename=outdir/"cha_c_renaming.tsv", chain="C", resid=201),
]

def element_of(atom_type, atom_name):
    e = atom_type.split(".")[0]
    e = "".join(c for c in e if c.isalpha())
    if e:
        return e.upper() if len(e) == 1 else e[:2].capitalize()
    e = "".join(c for c in atom_name if c.isalpha())
    return e[:1].upper() if e else "X"

def read_mol2(path):
    name, atoms, bonds, sec = "CHA", [], [], None
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>MOLECULE"):
            sec = "MOL"; 
            if i+1 < len(lines): name = lines[i+1].strip()
            continue
        if line.startswith("@<TRIPOS>ATOM"): sec = "ATOM"; continue
        if line.startswith("@<TRIPOS>BOND"): sec = "BOND"; continue
        if line.startswith("@<TRIPOS>"): sec = None; continue
        if sec == "ATOM" and len(line.split()) >= 9:
            p = line.split()
            atoms.append(dict(id=int(p[0]), old=p[1], x=float(p[2]), y=float(p[3]),
                              z=float(p[4]), type=p[5], charge=float(p[8])))
        elif sec == "BOND" and len(line.split()) >= 4:
            p = line.split()
            bonds.append(dict(id=int(p[0]), a1=int(p[1]), a2=int(p[2]), order=p[3]))
    return name, atoms, bonds

def rename_unique(atoms):
    counts = defaultdict(int); out = []
    for a in atoms:
        el = element_of(a["type"], a["old"]); counts[el] += 1
        b = dict(a); b["element"] = el; b["name"] = f"{el}{counts[el]}"
        if len(b["name"]) > 4: sys.exit(f"FAIL name too long: {b['name']}")
        out.append(b)
    dups = [n for n, c in Counter(a["name"] for a in out).items() if c > 1]
    if dups: sys.exit(f"FAIL duplicate names after rename: {dups}")
    return out

def write_mol2(path, atoms, bonds, name="CHA"):
    with path.open("w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{name}\n")
        f.write(f"{len(atoms):5d} {len(bonds):5d}     1     0     0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for a in atoms:
            f.write(f"{a['id']:7d} {a['name']:<8s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} "
                    f"{a['type']:<8s} {1:4d} {name:<8s} {a['charge']:10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for b in bonds:
            f.write(f"{b['id']:6d} {b['a1']:5d} {b['a2']:5d} {b['order']}\n")

audit_rows, all_ligands = [], []
for s in lig_specs:
    if not s["src"].exists(): sys.exit(f"FAIL missing placed ligand: {s['src']}")
    name, atoms, bonds = read_mol2(s["src"])
    if len(atoms) != 24: sys.exit(f"FAIL {s['label']} expected 24 atoms, got {len(atoms)}")
    if len(bonds) != 24: sys.exit(f"FAIL {s['label']} expected 24 bonds, got {len(bonds)}")
    renamed = rename_unique(atoms)
    write_mol2(s["out"], renamed, bonds)
    with s["rename"].open("w") as f:
        f.write("id\told_name\tnew_name\telement\ttype\tcharge\n")
        for a in renamed:
            f.write(f"{a['id']}\t{a['old']}\t{a['name']}\t{a['element']}\t{a['type']}\t{a['charge']:.6f}\n")
    net = sum(a["charge"] for a in renamed)
    heavy = [a for a in renamed if a["element"] != "H"]
    audit_rows.append(dict(label=s["label"], atoms=len(renamed), heavy=len(heavy),
                           bonds=len(bonds), net=net, rounded=round(net),
                           first=",".join(a["name"] for a in renamed[:10])))
    all_ligands.append((s, renamed))

with audit_tsv.open("w") as f:
    f.write("label\tatoms\theavy\tbonds\tnet_charge\trounded\tfirst_10_names\n")
    for r in audit_rows:
        f.write(f"{r['label']}\t{r['atoms']}\t{r['heavy']}\t{r['bonds']}\t"
                f"{r['net']:.6f}\t{r['rounded']}\t{r['first']}\n")

# rebuild combined PDB with renamed ligand atoms
prot_lines, natoms = [], 0
for line in protein.read_text().splitlines():
    if line.startswith("ATOM  "): prot_lines.append(line); natoms += 1
    elif line.startswith("TER"): prot_lines.append(line)

serial = natoms + 1
with combined_out.open("w") as f:
    for line in prot_lines: f.write(line + "\n")
    for s, atoms in all_ligands:
        for a in atoms:
            f.write(f"HETATM{serial:5d} {a['name'][:4]:<4s} CHA "
                    f"{s['chain']:1s}{s['resid']:4d}    "
                    f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                    f"{1.00:6.2f}{0.00:6.2f}          {a['element']:>2s}\n")
            serial += 1
        f.write("TER\n")
    f.write("END\n")

# audit combined PDB for ligand duplicate names
lig_names = defaultdict(list); nhet = 0
for line in combined_out.read_text().splitlines():
    if line.startswith("HETATM"):
        nhet += 1
        lig_names[(line[21], int(line[22:26]))].append(line[12:16].strip())
bad = [(k, [n for n, c in Counter(v).items() if c > 1]) for k, v in lig_names.items()
       if any(c > 1 for c in Counter(v).values())]
pdb_audit.write_text(
    f"combined_pdb\t{combined_out}\nprotein_atoms\t{natoms}\nhetatm_atoms\t{nhet}\n"
    f"ligand_residues\t{len(lig_names)}\nligand_duplicate_name_residues\t{len(bad)}\n")

print("ligand mol2 audit:")
for r in audit_rows:
    print(f"  {r['label']}: atoms={r['atoms']} heavy={r['heavy']} bonds={r['bonds']} "
          f"net_charge={r['net']:.6f} rounded={r['rounded']} first={r['first']}")
print(f"\ncombined PDB: protein_atoms={natoms} hetatm={nhet} "
      f"ligand_residues={len(lig_names)} duplicate_name_residues={len(bad)}")

for r in audit_rows:
    if r["rounded"] != -2: sys.exit(f"FAIL {r['label']} charge {r['rounded']} != -2")
if bad: sys.exit(f"FAIL duplicate ligand names remain: {bad}")
print("\nall freeze gates passed")
print(f"WROTE {combined_out}")
PY

# acceptance decision (carries corrected Arg63 provenance findings)
cat > "$admin/step06_acceptance_decision.txt" << 'EOF'
Step 06 - accepted pre-protonation complex

Accepted files (frozen):
  02_preparation/accepted_preprotonation/abc_cha_preprotonation.pdb
  02_preparation/accepted_preprotonation/cha_a.mol2
  02_preparation/accepted_preprotonation/cha_b.mol2
  02_preparation/accepted_preprotonation/cha_c.mol2

Protein scaffold:
  Paper-frame A/B/C trimer (2CHT core + 1DBF terminal repair), residues 1-127
  per chain, heavy-atom only, no altlocs, no duplicate atom names, backbone
  continuity intact, all active-site residues present.

Ligand:
  Chorismate placed by superposition onto the crystallographic TSA pose
  (TSA-fit RMSD 0.08-0.15 A). Net charge -2. CP2K duplicate atom names (C,H,O)
  renamed uniquely (C1,H1,O1,...) for AMBER/tleap safety; coordinates, atom
  types, charges, and bonds preserved.

Active-site architecture (corrected, proven):
  The active site is inter-subunit; each site is completed by residues from the
  adjacent chain (Arg63', Lys60', Thr74', Cys75', Phe57', Ala59'). Contact
  measurement across all chains recovered this with zero mismatches.

  The site-C Arg63' rotamer (guanidinium swung away, ~6-8 A vs ~3 A at A/B) is
  INHERITED CRYSTALLOGRAPHIC heterogeneity, not a workflow artifact:
    - raw 2CHT shows Arg63' conformational heterogeneity across its 12 TSA sites
      (guanidinium 2.9-7.3 A; 7/12 swung); the cha_c source site was already
      swung in the crystal.
    - repaired A/B/C Arg63 side chains are identical to raw 2CHT (RMSD 0.000 A).
  Expected to relax during MD equilibration.
EOF
cat "$admin/step06_acceptance_decision.txt"

echo
echo "=== step 06: checksums ==="
sha256sum \
  "$outdir/abc_cha_preprotonation.pdb" \
  "$outdir/cha_a.mol2" "$outdir/cha_b.mol2" "$outdir/cha_c.mol2" \
  > "$admin/sha256_step06_accepted_preprotonation.txt"
cat "$admin/sha256_step06_accepted_preprotonation.txt"

echo
echo "=== step 06: output files ==="
for f in \
  "$outdir/abc_cha_preprotonation.pdb" \
  "$outdir/cha_a.mol2" "$outdir/cha_b.mol2" "$outdir/cha_c.mol2" \
  "$outdir/cha_a_renaming.tsv" "$outdir/cha_b_renaming.tsv" "$outdir/cha_c_renaming.tsv" \
  "$admin/step06_acceptance_decision.txt" \
  "$admin/step06_ligand_mol2_audit.tsv" \
  "$admin/step06_combined_pdb_audit.txt" \
  "$admin/sha256_step06_accepted_preprotonation.txt"
do
  [[ -s "$f" ]] || { echo "FAIL missing/empty: $f"; exit 1; }
  echo "PASS $f"
done
echo "STEP 06 DONE"

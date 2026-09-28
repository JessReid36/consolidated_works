#!/usr/bin/env bash
# Step 04 - place chorismate into the repaired A/B/C active sites.
#
# Method: each chorismate ligand is positioned by least-squares (Kabsch)
# superposition onto the crystallographic transition-state-analogue (TSA) pose
# of its target A/B/C active site. CP2K ligand coordinates are registered to the
# J/K/L copy of 2CHT (step 01b); all 12 chains are crystallographically equivalent
# copies, so superposing from J/K/L onto A/B/C introduces no bias.
#
# Active site is INTER-SUBUNIT (Chook et al. 1993/1994): each site sits at the
# interface of two adjacent monomers. The active-site residue set is the published
# Agbaglo/DeYonker QM-cluster set (their SI, scheme S3): same-subunit Arg7, Arg90,
# Glu78, Tyr108, Leu115; adjacent-subunit (') Arg63', Lys60', Val73', Thr74',
# Cys75', Phe57', Ala59'. Arg116 is NOT in that set - it is second shell (measured
# here at ~5.5-6.8 A vs 2.6-3.4 A for true contacts) and is reported separately for
# provenance only. Catalytic-contact distances are measured against ALL chains,
# reporting the closest copy and its chain, so the same-chain vs cross-chain (')
# contribution is explicit and checkable against the literature.
#
# Distances are CLOSEST HEAVY-ATOM approaches on unrelaxed placeholder coords -
# they show neighbourhood/architecture, not precise interaction geometry. Not gated.
#
# NOTE: combined PDB carries no bond/charge records; the placed .mol2 files remain
# authoritative for connectivity and formal charge (chorismate net charge = -2).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
raw2cht="$root/01_inputs/structures/2cht_raw.pdb"
protein="$root/02_preparation/protein_only/abc_repaired_clean.pdb"
ligdir="$root/01_inputs/ligands"
outdir="$root/02_preparation/ligand_placement"
admin="$root/00_admin"
mkdir -p "$outdir" "$admin"

python3 - "$raw2cht" "$protein" "$ligdir" "$outdir" "$admin" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict
import numpy as np

raw2cht, protein, ligdir, outdir, admin = (Path(sys.argv[1]), Path(sys.argv[2]),
    Path(sys.argv[3]), Path(sys.argv[4]), Path(sys.argv[5]))

transform_rep = admin  / "step04_placement_transform_report.tsv"
contact_rep   = admin  / "step04_placement_contact_report.txt"
audit_rep     = admin  / "step04_placement_audit.txt"
catalytic_rep = admin  / "step04_catalytic_contact_report.tsv"
combined_pdb  = outdir / "abc_with_chorismate_unprotonated.pdb"

# Active-site residues = the published Agbaglo/DeYonker QM-cluster set (SI, S3),
# with EXPECTED chain origin relative to the substrate's own subunit. "same" = own
# chain; "cross" = adjacent chain ('). Code measures against ALL chains and reports
# the observed closest chain, so these labels are tested.
active_site = [
    (7,   "ARG", "same"),  (90,  "ARG", "same"),  (78,  "GLU", "same"),
    (108, "TYR", "same"),  (115, "LEU", "same"),
    (63,  "ARG", "cross"), (60,  "LYS", "cross"), (73,  "VAL", "cross"),
    (74,  "THR", "cross"), (75,  "CYS", "cross"), (57,  "PHE", "cross"),
    (59,  "ALA", "cross"),
]

# Second-shell residues: measured and reported for provenance, but NOT part of the
# active-site/QM-cluster set. Arg116 sits ~5.5-6.8 A from the substrate here and is
# absent from the Agbaglo set, so it is classified second shell rather than a contact.
second_shell = [
    (116, "ARG", "same"),
]

placements = [
    dict(lig="liga.mol2", out="cha_a_placed.mol2", src=("K",210), tgt=("A",203),
         lig_chain="A", lig_resid=201, resname="CHA", label="cha_a"),
    dict(lig="ligb.mol2", out="cha_b_placed.mol2", src=("J",212), tgt=("B",201),
         lig_chain="B", lig_resid=201, resname="CHA", label="cha_b"),
    dict(lig="ligc.mol2", out="cha_c_placed.mol2", src=("L",211), tgt=("C",202),
         lig_chain="C", lig_resid=201, resname="CHA", label="cha_c"),
]

def read_tsa(path):
    g = defaultdict(dict)
    for line in path.read_text().splitlines():
        if line.startswith("HETATM") and line[17:20].strip() == "TSA":
            key = (line[21], int(line[22:26])); atom = line[12:16].strip()
            g[key][atom] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return g

def read_mol2_atoms(path):
    atoms, sec = [], None
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"): sec = line.strip(); continue
        if sec == "@<TRIPOS>ATOM" and len(line.split()) >= 6:
            p = line.split()
            atoms.append(dict(id=int(p[0]), name=p[1], x=float(p[2]), y=float(p[3]),
                              z=float(p[4]), type=p[5],
                              charge=float(p[8]) if len(p) > 8 else 0.0))
    return atoms

def read_mol2_bonds(path):
    bonds, sec = [], None
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"): sec = line.strip(); continue
        if sec == "@<TRIPOS>BOND" and line.strip():
            bonds.append(line)
    return bonds

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    C = (P - Pc).T @ (Q - Qc)
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1, 1, d]) @ Wt
    rmsd = math.sqrt((((P - Pc) @ U + Qc - Q) ** 2).sum() / len(P))
    return Pc, Qc, U, rmsd

def elem_of(a):
    e = a["type"].split(".")[0]
    e = "".join(c for c in e if c.isalpha())
    return (e[:2].capitalize() if len(e) > 1 else e.upper()) if e else "X"

def is_h(name, elem): return elem.upper() == "H" or name.strip().startswith("H")
def dist(a, b): return math.sqrt(float(((a - b) ** 2).sum()))

def read_protein(path):
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM  "):
            out.append(dict(atom=line[12:16].strip(), resname=line[17:20].strip(),
                chain=line[21], resid=int(line[22:26]),
                xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                element=line[76:78].strip() if len(line) >= 78 else ""))
    return out

tsa = read_tsa(raw2cht)
prot = read_protein(protein)
assert prot, "no protein atoms"

placed_all, trows, clines, alines, catrows = [], [], [], [], []

for pl in placements:
    if pl["src"] not in tsa: sys.exit(f"FAIL missing source TSA {pl['src']}")
    if pl["tgt"] not in tsa: sys.exit(f"FAIL missing target TSA {pl['tgt']}")
    src, tgt = tsa[pl["src"]], tsa[pl["tgt"]]
    common = sorted(set(src) & set(tgt))
    if len(common) < 6: sys.exit(f"FAIL too few common TSA atoms {pl['src']}->{pl['tgt']}: {len(common)}")

    P = np.array([src[n] for n in common]); Q = np.array([tgt[n] for n in common])
    Pc, Qc, U, rmsd = kabsch(P, Q)

    atoms = read_mol2_atoms(ligdir / pl["lig"])
    if not atoms: sys.exit(f"FAIL no atoms in {pl['lig']}")
    placed = []
    for a in atoms:
        q = (np.array([a["x"], a["y"], a["z"]]) - Pc) @ U + Qc
        b = dict(a); b["x"], b["y"], b["z"] = map(float, q); placed.append(b)

    bonds = read_mol2_bonds(ligdir / pl["lig"])
    with (outdir / pl["out"]).open("w") as f:
        f.write("@<TRIPOS>MOLECULE\n" + pl["resname"] + "\n")
        f.write(f"{len(placed):5d} {len(bonds):5d}     1     0     0\nSMALL\nUSER_CHARGES\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for a in placed:
            f.write(f"{a['id']:7d} {a['name']:<8s} {a['x']:10.4f} {a['y']:10.4f} {a['z']:10.4f} "
                    f"{a['type']:<8s} {1:4d} {pl['resname']:<8s} {a['charge']:10.6f}\n")
        f.write("@<TRIPOS>BOND\n")
        for bl in bonds: f.write(bl + "\n")

    centre = np.array([[a["x"], a["y"], a["z"]] for a in placed]).mean(0)
    tgt_centre = np.array(list(tgt.values())).mean(0)
    offset = dist(centre, tgt_centre)

    heavy = [a for a in placed if not is_h(a["name"], elem_of(a))]
    heavy_xyz = [np.array([a["x"], a["y"], a["z"]]) for a in heavy]

    close, min_d, min_pair = {}, None, None
    for L in heavy_xyz:
        for pa in prot:
            if is_h(pa["atom"], pa["element"]): continue
            dd = dist(L, pa["xyz"])
            if min_d is None or dd < min_d:
                min_d = dd; min_pair = (pa["chain"], pa["resid"], pa["resname"], pa["atom"], dd)
            if dd <= 4.0:
                k = (pa["chain"], pa["resid"], pa["resname"])
                if k not in close or dd < close[k]: close[k] = dd

    # catalytic contacts: for each active-site (and second-shell) residue, find the
    # closest copy ACROSS ALL CHAINS and report which chain + same/cross label. The
    # shell tag keeps second-shell Arg116 in the record without treating it as a
    # defining active-site contact.
    for shell, residue_set in (("active_site", active_site), ("second_shell", second_shell)):
        for resid, exp_resname, exp_origin in residue_set:
            cands = [pa for pa in prot if pa["resid"] == resid
                     and not is_h(pa["atom"], pa["element"])]
            if not cands:
                catrows.append((pl["label"], pl["lig_chain"], shell, resid, exp_resname,
                                exp_origin, "NA", "NA", "NA", "absent")); continue
            best = None
            for L in heavy_xyz:
                for pa in cands:
                    dd = dist(L, pa["xyz"])
                    if best is None or dd < best[0]:
                        best = (dd, pa["chain"], pa["resname"], pa["atom"])
            dd, pchain, prname, patom = best
            observed = "same" if pchain == pl["lig_chain"] else "cross"
            match = "ok" if observed == exp_origin else "MISMATCH"
            catrows.append((pl["label"], pl["lig_chain"], shell, resid, prname, exp_origin,
                            pchain, observed, f"{dd:.3f}", f"{patom}:{match}"))

    trows.append((pl["label"], pl["lig"], f"{pl['src'][0]}{pl['src'][1]}",
                  f"{pl['tgt'][0]}{pl['tgt'][1]}", len(common), rmsd, offset))
    clines.append(f"contacts within 4.0 A of {pl['label']}")
    clines.append("-" * 50)
    for (ch, r, rn), dd in sorted(close.items()):
        clines.append(f"{ch} {r:4d} {rn:3s}  min_dist={dd:6.3f}")
    clines.append("")
    alines += [f"{pl['label']}_atoms\t{len(placed)}",
               f"{pl['label']}_heavy\t{len(heavy)}",
               f"{pl['label']}_tsa_fit_rmsd_A\t{rmsd:.4f}",
               f"{pl['label']}_centre_offset_A\t{offset:.3f}",
               f"{pl['label']}_min_prot_lig_heavy_A\t{min_d:.3f}"]
    if min_pair:
        pc, pr, prn, pan, dd = min_pair
        alines.append(f"{pl['label']}_min_pair\tprot:{pc}:{pr}:{prn}:{pan}:{dd:.3f}")
    if min_d is not None and min_d < 0.80:
        sys.exit(f"FAIL severe overlap for {pl['label']}: {min_d:.3f} A")
    placed_all.append((pl, placed))

with transform_rep.open("w") as f:
    f.write("label\tinput\tsource_tsa\ttarget_tsa\tcommon_atoms\tfit_rmsd_A\tcentre_offset_A\n")
    for lab, lig, s, t, n, r, o in trows:
        f.write(f"{lab}\t{lig}\t{s}\t{t}\t{n}\t{r:.4f}\t{o:.3f}\n")
contact_rep.write_text("\n".join(clines) + "\n")
audit_rep.write_text("step04_placement_audit\n" + "\n".join(alines) + "\n")

with catalytic_rep.open("w") as f:
    f.write("# Active-site set = published Agbaglo QM-cluster residues (S3); 'cross' residues are contributed by the adjacent chain (').\n")
    f.write("# shell=active_site are the QM-cluster residues; shell=second_shell (Arg116) is measured for provenance, not a defining contact.\n")
    f.write("ligand\tlig_chain\tshell\tresid\tresname\texpected_origin\tobserved_chain\tobserved_origin\tclosest_heavy_A\tpartner_atom_match\n")
    for row in catrows:
        f.write("\t".join(str(x) for x in row) + "\n")

with combined_pdb.open("w") as f:
    for line in protein.read_text().splitlines():
        if line.startswith("ATOM  ") or line.startswith("TER"):
            f.write(line + "\n")
    serial = len(prot) + 1
    for pl, placed in placed_all:
        for a in placed:
            f.write(f"HETATM{serial:5d} {a['name'][:4]:<4s} {pl['resname']:>3s} "
                    f"{pl['lig_chain']:1s}{pl['lig_resid']:4d}    "
                    f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                    f"{1.00:6.2f}{0.00:6.2f}          {elem_of(a):>2s}\n")
            serial += 1
        f.write("TER\n")
    f.write("END\n")

print("placement report:")
for lab, lig, s, t, n, r, o in trows:
    print(f"  {lab}: {lig} {s}->{t} | fit RMSD {r:.4f} A | centre offset {o:.3f} A")

print("\ncatalytic-residue closest heavy-atom approach across all chains:")
print(f"  {'site':4s} {'shell':12s} {'res':7s} {'exp':5s} {'obs_chain':9s} {'obs':5s} {'dist_A':7s} {'check'}")
for pl in placements:
    for row in catrows:
        if row[0] != pl["label"]: continue
        _, lch, shell, resid, rname, exp, obch, obs, dd, match = row
        flag = "" if match.endswith("ok") else "  <-- unexpected"
        print(f"  {lch:4s} {shell:12s} {rname}{resid:<4d} {exp:5s} {obch:9s} {obs:5s} {dd:>7s} {match}{flag}")

print(f"\nWROTE {catalytic_rep}")
print(f"WROTE {combined_pdb}")
PY
echo "STEP 04 DONE"

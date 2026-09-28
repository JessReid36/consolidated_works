#!/usr/bin/env bash
# Step 02 - repair missing terminal residues of the 2CHT A/B/C trimer using 1DBF.
#
# Rationale: 2CHT chains A/B/C lack ordered terminal residues. Following the
# reference MD-preparation protocol, missing terminal residues are transferred
# from 1DBF, which is globally superposed onto 2CHT A/B/C first. Grafted termini
# that do not meet the backbone at a physical C-N distance are closed by minimal
# rigid translation of the terminal block; this is geometric pre-conditioning
# only - bond angles/dihedrals at the graft are relaxed later by minimisation
# and MD equilibration, and the affected joins must be re-checked post-minimisation.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

root="$HOME/system_development"
struct="$root/01_inputs/structures"
outdir="$root/02_preparation/protein_only"
admin="$root/00_admin"
mkdir -p "$outdir" "$admin"

python3 - "$struct/2cht_raw.pdb" "$struct/1dbf_raw.pdb" "$outdir" "$admin" <<'PY'
import sys, math
from pathlib import Path
from collections import defaultdict
import numpy as np

pdb2, pdb1, outdir, admin = (Path(sys.argv[1]), Path(sys.argv[2]),
                             Path(sys.argv[3]), Path(sys.argv[4]))
target_chains = ["A", "B", "C"]

raw_out       = outdir / "2cht_abc_raw_protein.pdb"
aligned_out   = outdir / "1dbf_abc_aligned_to_2cht.pdb"
unadj_out     = outdir / "abc_repaired_unadjusted.pdb"
final_out     = outdir / "abc_repaired.pdb"
align_rep     = admin  / "step02_alignment_report.tsv"
cont_rep      = admin  / "step02_peptide_continuity_report.tsv"
adj_rep       = admin  / "step02_join_adjustment_report.tsv"
src_rep       = admin  / "step02_residue_source_report.tsv"
decision_rep  = admin  / "step02_repair_decision.txt"

# Terminal residues to take from aligned 1DBF (paper-frame patch definition).
patch_resids = {
    "A": set([1] + list(range(116, 128))),
    "B": set([1] + list(range(115, 128))),
    "C": set([1, 2] + list(range(115, 128))),
}
# Graft joins permitted to receive rigid-translation closure. Any OTHER bad
# join is treated as a hard failure (guards against silent bad geometry).
allowed_joins = {("A",115,116), ("B",1,2), ("B",114,115), ("C",2,3), ("C",114,115)}
TARGET_CN = 1.329  # canonical peptide C-N bond length (A)

def parse(path):
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM  "):
            continue
        atoms.append({
            "name": line[12:16], "altloc": line[16], "resname": line[17:20],
            "chain": line[21], "resid": int(line[22:26]), "icode": line[26],
            "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
            "occ": line[54:60] if len(line) >= 60 else "  1.00",
            "bfac": line[60:66] if len(line) >= 66 else "  0.00",
            "element": line[76:78] if len(line) >= 78 else "",
        })
    return atoms

def rkey(a):  return (a["chain"], a["resid"], a["icode"], a["resname"].strip())
def akey(a):  return (a["chain"], a["resid"], a["name"].strip())
def rsort(k): return (target_chains.index(k[0]), k[1], k[2], k[3])

def fmt(a, serial):
    return (f"ATOM  {serial:5d} {a['name']}{a['altloc']}{a['resname']:>3s} "
            f"{a['chain']:1s}{a['resid']:4d}{a['icode']:1s}   "
            f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            f"{a['occ']:>6s}{a['bfac']:>6s}          {a['element']:>2s}")

def write(path, resdict):
    serial, chain = 1, None
    with path.open("w") as f:
        for k in sorted(resdict, key=rsort):
            if chain is not None and k[0] != chain:
                f.write("TER\n")
            chain = k[0]
            for a in resdict[k]:
                f.write(fmt(a, serial) + "\n"); serial += 1
        f.write("TER\nEND\n")

def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    C = (P - Pc).T @ (Q - Qc)
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    U = V @ np.diag([1, 1, d]) @ Wt
    rmsd = math.sqrt((((P - Pc) @ U + Qc - Q) ** 2).sum() / len(P))
    return Pc, Qc, U, rmsd

def xform(a, Pc, Qc, U):
    p = np.array([a["x"], a["y"], a["z"]])
    q = (p - Pc) @ U + Qc
    b = dict(a); b["x"], b["y"], b["z"] = map(float, q); return b

def dist(a, b): return math.sqrt(float(((a - b) ** 2).sum()))

a2 = [a for a in parse(pdb2) if a["chain"] in target_chains]
a1 = [a for a in parse(pdb1) if a["chain"] in target_chains]
assert a2 and a1, "missing A/B/C atoms in input"

existing = defaultdict(list)
for a in a2:
    b = dict(a); b["_src"] = "2cht"; existing[rkey(b)].append(b)
write(raw_out, existing)

m2 = {akey(a): a for a in a2}
m1 = {akey(a): a for a in a1}
src, tgt = [], []
for ch in target_chains:
    for r in range(1, 128):
        for nm in ("N", "CA", "C", "O"):
            k = (ch, r, nm)
            if k in m1 and k in m2 and m1[k]["resname"].strip() == m2[k]["resname"].strip():
                src.append([m1[k]["x"], m1[k]["y"], m1[k]["z"]])
                tgt.append([m2[k]["x"], m2[k]["y"], m2[k]["z"]])
assert len(src) >= 300, f"too few common backbone atoms: {len(src)}"
Pc, Qc, U, rmsd = kabsch(np.array(src), np.array(tgt))
print(f"1DBF->2CHT A/B/C backbone superposition: {len(src)} atoms, RMSD {rmsd:.4f} A")
align_rep.write_text("alignment\tcommon_backbone_atoms\tRMSD_A\n"
                     f"1dbf_abc_to_2cht_abc\t{len(src)}\t{rmsd:.4f}\n")
assert rmsd <= 2.0, f"alignment RMSD too high: {rmsd:.4f}"

aligned = defaultdict(list)
for a in a1:
    b = xform(a, Pc, Qc, U); b["_src"] = "1dbf_aligned"; aligned[rkey(b)].append(b)
write(aligned_out, aligned)

# Build repaired structure: 2CHT core, terminal patch residues from aligned 1DBF.
repaired = defaultdict(list)
for k, v in existing.items():
    repaired[k] = [dict(a) for a in v]
for k in list(repaired):
    ch, r, ic, rn = k
    if ch in patch_resids and r in patch_resids[ch]:
        del repaired[k]
for k, v in aligned.items():
    ch, r, ic, rn = k
    if ch in patch_resids and r in patch_resids[ch]:
        repaired[k] = [dict(a) for a in v]
write(unadj_out, repaired)

def coords(resdict):
    d = defaultdict(dict)
    for k, v in resdict.items():
        for a in v:
            d[(k[0], k[1])][a["name"].strip()] = np.array([a["x"], a["y"], a["z"]])
    return d

def continuity(resdict):
    c = coords(resdict); bad = []; checked = 0
    by_ch = defaultdict(list)
    for (ch, r) in c: by_ch[ch].append(r)
    for ch in sorted(by_ch):
        rs = sorted(set(by_ch[ch]))
        for r1, r2 in zip(rs[:-1], rs[1:]):
            if r2 != r1 + 1:
                bad.append((ch, r1, r2, "noncontiguous", None)); continue
            C, N = c[(ch, r1)].get("C"), c[(ch, r2)].get("N")
            if C is None or N is None:
                bad.append((ch, r1, r2, "missing_C_or_N", None)); continue
            d = dist(C, N); checked += 1
            if d < 1.15 or d > 1.70:
                bad.append((ch, r1, r2, "abnormal_CN", d))
    return checked, bad

def shift_block(resdict, ch, resids, s):
    for k, v in resdict.items():
        if k[0] == ch and k[1] in resids:
            for a in v:
                a["x"] += float(s[0]); a["y"] += float(s[1]); a["z"] += float(s[2])

checked0, bad0 = continuity(repaired)
print(f"\nunadjusted repair: {checked0} links checked, {len(bad0)} bad")
for ch, r1, r2, why, d in bad0:
    print(f"  {ch} {r1}->{r2}: {why}" + (f" d={d:.3f}" if d else ""))
    if (ch, r1, r2) not in allowed_joins:
        sys.exit(f"FAIL unexpected bad join {ch} {r1}->{r2} - not in allowed graft set")

# Close allowed terminal joins by minimal rigid translation of the graft block.
adj_rows = []
for cycle in range(1, 6):
    checked, bad = continuity(repaired)
    if not bad:
        break
    ch, r1, r2, why, d = bad[0]
    c = coords(repaired)
    C, N = c[(ch, r1)]["C"], c[(ch, r2)]["N"]
    if r1 <= 2:                                   # N-terminal graft: move residues 1..r1
        direction = C - N; desired = N + direction / np.linalg.norm(direction) * TARGET_CN
        s = desired - C; block = set(range(1, r1 + 1)); label = f"{ch}:1-{r1}"
    else:                                         # C-terminal graft: move r2..127
        direction = N - C; desired = C + direction / np.linalg.norm(direction) * TARGET_CN
        s = desired - N; block = set(range(r2, 128)); label = f"{ch}:{r2}-127"
    shift_block(repaired, ch, block, s)
    adj_rows.append((cycle, ch, r1, r2, d, float(np.linalg.norm(s)), label))

checkedF, badF = continuity(repaired)
adj_rep.write_text("cycle\tchain\tr1\tr2\tinitial_CN_A\tshift_mag_A\tmoved_block\n" +
    "".join(f"{c}\t{ch}\t{r1}\t{r2}\t{d:.3f}\t{sm:.3f}\t{lab}\n"
            for c, ch, r1, r2, d, sm, lab in adj_rows))
cont_rep.write_text("stage\tlinks_checked\tbad_links\n"
                    f"unadjusted\t{checked0}\t{len(bad0)}\n"
                    f"adjusted\t{checkedF}\t{len(badF)}\n")

if badF:
    for ch, r1, r2, why, d in badF:
        print(f"  STILL BAD {ch} {r1}->{r2}: {why}" + (f" d={d:.3f}" if d else ""))
    sys.exit("FAIL adjusted repair still has bad peptide links")

write(final_out, repaired)
with src_rep.open("w") as f:
    f.write("chain\tresid\tresname\tsource\tatoms\n")
    for k in sorted(repaired, key=rsort):
        srcs = ",".join(sorted(set(a["_src"] for a in repaired[k])))
        f.write(f"{k[0]}\t{k[1]}\t{k[3]}\t{srcs}\t{len(repaired[k])}\n")

decision_rep.write_text(
    "Step 02 - A/B/C terminal repair using 1DBF\n\n"
    f"Global backbone superposition RMSD: {rmsd:.4f} A ({len(src)} atoms)\n"
    f"Unadjusted bad peptide links: {len(bad0)}\n"
    f"Adjusted bad peptide links: {len(badF)}\n"
    f"Peptide links checked: {checkedF}\n\n"
    "Terminal graft joins were closed by minimal rigid translation to a canonical\n"
    "C-N distance. This is geometric pre-conditioning; graft bond angles/dihedrals\n"
    "are relaxed by downstream minimisation and MD and must be re-checked then.\n")

print("\nterminal join closures:")
for c, ch, r1, r2, d, sm, lab in adj_rows:
    print(f"  {ch} {r1}->{r2}: {d:.3f} A -> {TARGET_CN} A; moved {lab} by {sm:.3f} A")

# Structural summary + active-site presence.
required = [7,57,59,60,63,73,74,75,78,90,108,115]
seen = defaultdict(set); counts = defaultdict(int)
for k in repaired:
    for a in repaired[k]:
        seen[k[0]].add((k[1], k[3])); counts[k[0]] += 1
print("\nchain summary:")
for ch in target_chains:
    rs = sorted(r for r, rn in seen[ch])
    print(f"  {ch}: {counts[ch]} atoms, {len(rs)} residues, {min(rs)}-{max(rs)}")
missing = [(ch, r) for ch in target_chains for r in required
           if not any(rr == r for rr, rn in seen[ch])]
if missing:
    sys.exit(f"FAIL active-site residues missing: {missing}")
print("active-site residues (7,57,59,60,63,73,74,75,78,90,108,115): all present in A/B/C")
print(f"\nWROTE {final_out}")
PY
echo "STEP 02 DONE"

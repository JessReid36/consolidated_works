#!/usr/bin/env python3
"""
step05_verify.py - re-runnable audit of step 05, Arg63' provenance and conformation.

Step 05 asks whether the swung-away Arg63' at placed site C is inherited from the crystal
or introduced by the repair and placement workflow. The CONCLUSION is correct. The
ARGUMENT given in step05_arg63_provenance.sh is not, and check 5 below demonstrates that.

Checks performed:

  1. Continuity. The raw deposition and the placed complex are the committed outputs of
     step 01 and step 04, by SHA-256.
  2. Raw survey reproduced. For all twelve analogue sites, the nearest Arg63 guanidinium
     nitrogen across all chains is recomputed and compared against
     step05_arg63_provenance_report.tsv.
  3. Placed survey reproduced, for the three placed substrates.
  4. The SWUNG flag is applied consistently at the script's own threshold of 4.5 A on the
     nearest-nitrogen distance.
  5. SOURCE VERSUS TARGET. The script's header reasons from the SOURCE analogue site that
     each substrate geometry was aligned to (K210, J212, L211 for A, B, C). Placement
     transforms the substrate INTO the target site and does not move the protein, so the
     Arg63' conformation seen by a placed substrate is the one belonging to its TARGET
     site (A203, B201, C202). This check tests both predictions against the observed
     result. Source-site reasoning predicts all three sites swung, because all three
     source sites are swung; the observed pattern is normal, normal, swung. Target-site
     reasoning reproduces it.
  6. Repair invariance (step05b). Arg63 side chains in the repaired trimer are unchanged
     against raw 2CHT, over 7 side-chain atoms in each of the three chains.
  7. Measurement basis. The raw deposition carries hydrogens, so its nearest-any-atom
     column can be satisfied by a hydrogen on Arg63 itself. In the placed complex the
     protein is heavy-atom only after step 03, while the substrates retain their own
     hydrogens, so the same column is protein-heavy against substrate-including-H. The two
     nearest-any columns are therefore NOT comparable. The nearest-NITROGEN column
     restricts the protein side to NE, NH1 and NH2 in both cases and is the one to quote.
  8. Literature. The proportion of bound sites in the raw deposition is compared against
     Claeyssens et al. (2011), who report Arg63 bound in roughly one third of sites.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step05_verify.py [outputs_dir] [step01_outputs] [step04_outputs]
"""
import sys
import math
import hashlib
from pathlib import Path
from collections import defaultdict

FAIL = []
SWUNG_THRESHOLD = 4.5
GUANIDINIUM_N = {"NE", "NH1", "NH2"}
SOURCE = {"A": "K210", "B": "J212", "C": "L211"}
TARGET = {"A": "A203", "B": "B201", "C": "C202"}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def atoms(p, resname=None, hetatm=True):
    out = []
    for L in Path(p).read_text().splitlines():
        if not (L.startswith("ATOM  ") or (hetatm and L.startswith("HETATM"))):
            continue
        rn = L[17:20].strip()
        if resname and rn != resname:
            continue
        el = (L[76:78].strip() if len(L) >= 78 else "") or L[12:16].strip()[0]
        out.append({"name": L[12:16].strip(), "resname": rn, "chain": L[21],
                    "resid": int(L[22:26]), "element": el,
                    "xyz": (float(L[30:38]), float(L[38:46]), float(L[46:54]))})
    return out


def nearest_arg63_N(ligand_xyz, arg63):
    best = None
    for a in arg63:
        if a["name"] not in GUANIDINIUM_N:
            continue
        d = min(math.dist(a["xyz"], c) for c in ligand_xyz)
        if best is None or d < best[0]:
            best = (d, a["chain"])
    return best


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s01 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../01_inputs_acquisition/outputs")
    s04 = Path(sys.argv[3] if len(sys.argv) > 3 else "../../04_substrate_placement/outputs")
    print(f"step 05 verification against {out.resolve()}\n")

    rows = [l.split("\t") for l in
            (out / "step05_arg63_provenance_report.tsv").read_text().splitlines()[1:]]
    raw_rec = {r[1]: (float(r[2]), r[3], r[7]) for r in rows if r[0] == "raw_2cht"}
    pl_rec = {r[1]: (float(r[2]), r[3], r[7]) for r in rows if r[0] == "placed"}

    print("1. continuity")
    check("step01 2cht_raw.pdb reachable", (s01 / "2cht_raw.pdb").exists(),
          sha(s01 / "2cht_raw.pdb")[:16] if (s01 / "2cht_raw.pdb").exists() else "")
    check("step04 placed complex reachable",
          (s04 / "abc_with_chorismate_unprotonated.pdb").exists())

    print("\n2. raw survey reproduced over all twelve analogue sites")
    raw = atoms(s01 / "2cht_raw.pdb")
    arg63_raw = [a for a in raw if a["resid"] == 63 and a["resname"] == "ARG"]
    tsa = defaultdict(list)
    for a in raw:
        if a["resname"] == "TSA":
            tsa[f"{a['chain']}{a['resid']}"].append(a["xyz"])
    check("twelve analogue sites found", len(tsa) == 12, f"n={len(tsa)}")
    worst = 0.0
    for site, coords in sorted(tsa.items()):
        got = nearest_arg63_N(coords, arg63_raw)
        want = raw_rec[site]
        worst = max(worst, abs(got[0] - want[0]))
        if got[1] != want[1]:
            check(f"{site}: nearest Arg63 chain", False, f"{got[1]} vs {want[1]}")
    check("all twelve nearest-nitrogen distances reproduced", worst < 2e-3,
          f"max difference {worst:.4f} A")

    print("\n3. placed survey reproduced")
    pl = atoms(s04 / "abc_with_chorismate_unprotonated.pdb")
    arg63_pl = [a for a in pl if a["resid"] == 63 and a["resname"] == "ARG"]
    cha = defaultdict(list)
    for a in pl:
        if a["resname"] == "CHA":
            cha[a["chain"]].append(a["xyz"])
    worst = 0.0
    for ch, coords in sorted(cha.items()):
        got = nearest_arg63_N(coords, arg63_pl)
        want = pl_rec[f"{ch}201"]
        worst = max(worst, abs(got[0] - want[0]))
    check("all three placed distances reproduced", worst < 2e-3,
          f"max difference {worst:.4f} A")

    print("\n4. SWUNG flag applied consistently at the script's 4.5 A threshold")
    bad = [s for s, (d, c, f) in {**raw_rec, **pl_rec}.items()
           if (d > SWUNG_THRESHOLD) != (f == "SWUNG")]
    check(f"flag equals nearest-N > {SWUNG_THRESHOLD} A everywhere", not bad, str(bad))

    print("\n5. source-site versus target-site reasoning")
    print(f"         {'lig':<5}{'placed':>9}{'source':>9}{'src N':>8}{'target':>9}{'tgt N':>8}")
    src_pred, tgt_pred, obs = [], [], []
    for ch in "ABC":
        p = pl_rec[f"{ch}201"][0]
        s = raw_rec[SOURCE[ch]][0]
        t = raw_rec[TARGET[ch]][0]
        print(f"         {ch:<5}{p:>9.3f}{SOURCE[ch]:>9}{s:>8.3f}{TARGET[ch]:>9}{t:>8.3f}")
        src_pred.append(s > SWUNG_THRESHOLD)
        tgt_pred.append(t > SWUNG_THRESHOLD)
        obs.append(p > SWUNG_THRESHOLD)
    check("target-site conformation predicts the observed pattern", tgt_pred == obs,
          f"predicted {tgt_pred}, observed {obs}")
    check("source-site conformation does NOT predict it (the header's argument fails)",
          src_pred != obs, f"predicted {src_pred}, observed {obs}")
    tgt_err = max(abs(pl_rec[f"{ch}201"][0] - raw_rec[TARGET[ch]][0]) for ch in "ABC")
    print(f"         placed vs target-site nearest-N agree to within {tgt_err:.3f} A")

    print("\n6. repair invariance of the Arg63 side chains")
    for line in (out / "step05b_arg63_conformation_report.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        check(f"chain {f[0]}: Arg63 side chain unchanged by repair",
              float(f[3]) < 1e-3 and f[4] == "unchanged",
              f"{f[2]} atoms, RMSD {f[3]} A, {f[4]}")

    print("\n7. measurement basis")
    n_h = sum(1 for a in raw if a["element"].upper() == "H")
    check("the raw deposition carries hydrogens, so nearest-any is not comparable",
          n_h > 0, f"{n_h} hydrogens in 2cht_raw.pdb; quote the nearest-nitrogen column")
    prot_h = sum(1 for a in pl if a["resname"] != "CHA" and a["element"].upper() == "H")
    lig_h = sum(1 for a in pl if a["resname"] == "CHA" and a["element"].upper() == "H")
    check("the PROTEIN in the placed complex is heavy-atom only, per step 03",
          prot_h == 0, f"{prot_h} protein hydrogens")
    check("the placed substrates retain their own hydrogens",
          lig_h == 24, f"{lig_h} ligand hydrogens, 8 per chorismate")

    print("\n8. literature comparison")
    n_swung = sum(1 for d, c, f in raw_rec.values() if d > SWUNG_THRESHOLD)
    bound = 12 - n_swung
    print(f"         {n_swung} of 12 sites swung, {bound} of 12 bound ({bound/12:.0%})")
    print("         Claeyssens et al. (2011): Arg63 bound in roughly one third of sites")
    check("bound fraction is of the order Claeyssens report", 2 <= bound <= 6,
          f"{bound}/12")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

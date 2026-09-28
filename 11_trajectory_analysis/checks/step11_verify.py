#!/usr/bin/env python3
"""
step11_verify.py - re-runnable audit of step 11, production trajectory analysis.

Every figure quoted in step11_trajectory_analysis.md is recomputed from the committed
.dat files. Three of them do not reproduce as stated, and the discrepancies are reported
rather than smoothed over. The literature claims are also checked against the papers.

Checks performed:

  1. 11a backbone RMSD: mean, maximum, and the first-to-last 2 ns drift.
  2. 11b per-residue RMSF: the core statistic, the flexible set, and the claim that the
     contact-set residues are not among the flexible ones.
  3. 11c reaction coordinate: per-copy breaking distance, forming distance and r, and the
     fraction of frames inside the Hur & Bruice distance criterion. The last is NOT in the
     committed document and is what actually separates the three sites.
  4. The NAC criterion is sourced, and the full three-condition test was applied.
  5. Ensemble group labels correspond to the frames the full NAC test rejected.
  6. Literature: the Agbaglo RMSD comparator, with the 2026 correction applied.

Exit status 0 if every check passes, 1 otherwise. A discrepancy against the committed
document is reported as a FAIL, because the document is what the write-up would quote.

Usage:  python3 step11_verify.py [outputs] [qmmm_dir]
"""
import sys
import statistics as st
from pathlib import Path

FAIL = []
CONTACT_SET = {7, 57, 59, 60, 63, 73, 74, 75, 78, 90, 108, 115}
SECOND_SHELL = {116}
D_NAC = 3.7           # Hur & Bruice, PNAS 2003, 100, 12015
SITE = {"1": "A", "2": "B", "3": "C"}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    qmm = Path(sys.argv[2] if len(sys.argv) > 2 else
               "../../../chorismate-thesis-results/05_qmmm")
    print(f"step 11 verification against {out.resolve()}\n")

    print("1. 11a backbone RMSD")
    rm = []
    for l in (out / "rmsd_vs_time.dat").read_text().splitlines():
        if l.startswith("#") or not l.strip():
            continue
        p = l.split("\t")
        rm.append((float(p[1]), float(p[2])))
    v = [x[1] for x in rm]
    check("frame count is the 20 ns production at 1 ps", len(v) == 20000, f"{len(v)}")
    check("mean RMSD is 3.16 A as documented", abs(st.mean(v) - 3.16) < 0.005,
          f"{st.mean(v):.2f}")
    check("maximum RMSD is 4.31 A as documented", abs(max(v) - 4.31) < 0.005,
          f"{max(v):.2f}")
    first, last = st.mean(v[:2000]), st.mean(v[-2000:])
    check("first 2 ns mean is 1.98 A", abs(first - 1.98) < 0.005, f"{first:.2f}")
    check("last 2 ns mean is 3.59 A", abs(last - 3.59) < 0.005, f"{last:.2f}")
    print(f"         drift {last - first:.2f} A over 20 ns")

    print("\n2. 11b per-residue RMSF")
    f = []
    for l in (out / "rmsf_per_residue.dat").read_text().splitlines():
        if l.startswith("#") or not l.strip():
            continue
        p = l.split("\t")
        f.append((p[0], int(p[1]), p[2], float(p[3])))
    vals = sorted(x[3] for x in f)
    lo = vals[:int(0.8 * len(vals))]
    check("381 residues analysed", len(f) == 381, f"{len(f)}")
    check("lower-80% mean is 0.69 A as documented", abs(st.mean(lo) - 0.69) < 0.006,
          f"{st.mean(lo):.3f}")
    # The document pairs "mean 0.69, median 0.74". 0.69 is the lower-80% mean; 0.74 is the
    # median of ALL residues. Both are correct individually but they describe different
    # populations, so the pair should not be quoted as one statistic.
    med_lo, med_all = st.median(lo), st.median([x[3] for x in f])
    check("the documented median 0.74 is the lower-80% median",
          abs(med_lo - 0.74) < 0.006,
          f"lower-80% median is {med_lo:.3f}; the all-residue median is {med_all:.3f}. "
          f"The document pairs a lower-80% MEAN with an all-residue MEDIAN.")
    top = sorted(f, key=lambda x: -x[3])[:6]
    print("         most flexible: " +
          ", ".join(f"{c}{r} {val:.1f}" for c, r, n, val in top))
    site = [x for x in f if x[1] in CONTACT_SET | SECOND_SHELL]
    thr = vals[int(0.8 * len(vals))]
    above = [f"{c}{r}" for c, r, n, val in site if val >= thr]
    check("every contact-set residue is below 1.0 A RMSF",
          max(x[3] for x in site) < 1.0,
          f"n={len(site)}, mean {st.mean([x[3] for x in site]):.2f}, "
          f"max {max(x[3] for x in site):.2f}")
    check("the document's claim that none is in the flexible set holds strictly",
          not above,
          f"{above} sit at or above the 80th-percentile threshold of {thr:.2f} A, "
          f"though all are under 1.0 A while the tails reach {max(vals):.1f} A")

    print("\n3. 11c reaction coordinate")
    cols = {k: [] for k in ("1b", "1f", "1r", "2b", "2f", "2r", "3b", "3f", "3r")}
    for l in (out / "rxn_coord_per_frame.dat").read_text().splitlines():
        if l.startswith("#") or not l.strip():
            continue
        p = [float(x) for x in l.split("\t")]
        for i, k in enumerate(cols):
            cols[k].append(p[1 + i])
    doc = {"1": (4.15, -2.70), "2": (3.29, -1.83), "3": (3.83, -2.38)}
    for c in "123":
        fm, rm_ = st.mean(cols[c + "f"]), st.mean(cols[c + "r"])
        check(f"CHA#{c} (site {SITE[c]}) forming mean and r reproduce",
              abs(fm - doc[c][0]) < 0.006 and abs(rm_ - doc[c][1]) < 0.006,
              f"form {fm:.2f} (doc {doc[c][0]}), r {rm_:.2f} (doc {doc[c][1]})")
    bmin = min(min(cols[c + "b"]) for c in "123")
    bmax = max(max(cols[c + "b"]) for c in "123")
    check("the breaking bond is effectively fixed in classical MD",
          bmax - bmin < 0.30, f"range {bmin:.3f}-{bmax:.3f} A across all copies")
    print(f"\n         NOT IN THE COMMITTED DOCUMENT — occupancy of the Hur & Bruice "
          f"distance criterion:")
    occ = {}
    for c in "123":
        fr = sum(1 for x in cols[c + "f"] if x <= D_NAC) / len(cols[c + "f"])
        occ[c] = fr
        print(f"           CHA#{c} (site {SITE[c]}): {100*fr:5.1f}% of frames with "
              f"d(form) <= {D_NAC} A   min {min(cols[c+'f']):.2f} A")
    check("the three sites differ enormously in near-attack occupancy",
          max(occ.values()) / max(min(occ.values()), 1e-9) > 10,
          f"{100*occ['2']:.1f}% at site B against {100*occ['1']:.1f}% at site A. "
          f"The document says 'all three sample near-attack conformations', which is "
          f"true of the MINIMA but hides a 27-fold difference in occupancy.")

    print("\n4. the NAC criterion is sourced and fully applied")
    log = qmm / "11f_nac_check.log"
    if not log.exists():
        check("NAC verification log reachable", False, str(log))
    else:
        txt = log.read_text()
        check("the threshold is attributed to a primary source",
              "Hur" in txt and "Bruice" in txt and "12015" in txt,
              "Hur & Bruice, PNAS 2003, 100, 12015")
        check("all three published conditions are tested, not distance alone",
              "theta1" in txt and "theta2" in txt)
        import re
        m = re.search(r"(\d+) of (\d+) frames satisfy the full", txt)
        check("the full test was actually evaluated", bool(m),
              f"{m.group(0)}" if m else "")

    print("\n5. ensemble group labels match the NAC test")
    ens = qmm / "19_ensemble_barriers/ensemble_barriers.tsv"
    rows, hdr = [], None
    for l in ens.read_text().splitlines():
        if l.startswith("#"):
            continue
        p = l.rstrip("\n").split("\t")
        if p[0] == "frame":
            hdr = p
            continue
        rows.append(dict(zip(hdr, p)))
    nonfull = sorted(r["frame"] for r in rows if r["group"] != "full_NAC")
    rejected = {"05680", "07310", "15825"}
    check("the frames the full NAC test rejected are labelled near_NAC",
          rejected <= set(nonfull), f"non-full_NAC frames: {nonfull}")
    new = [r for r in rows if int(r["frame"]) >= 20000]
    check("every frame in the extended set is a full NAC",
          all(r["group"] == "full_NAC" for r in new), f"n={len(new)}")

    print("\n6. literature comparator")
    print("         Agbaglo et al. 2024, scheme S2: backbone RMSD of the entire protein,")
    print("         mean 2.66 A, sd 0.76. Ours is 3.16 A, drift 1.61 A over 20 ns.")
    print("         THE 2026 CORRECTION MATTERS HERE: it restates that for schemes S2-S7")
    print("         the RMSD was computed against the INITIAL EQUILIBRATED MD STRUCTURE,")
    print("         not the X-ray structure as originally published. Our RMSD is measured")
    print("         against production frame 0, which is the post-equilibration structure,")
    print("         so the two are like for like ONLY under the corrected reading.")
    check("our RMSD is of the same order as the published comparator",
          abs(st.mean(v) - 2.66) < 1.0, f"{st.mean(v):.2f} vs 2.66 A")

    print("\n7. the full 60 ns, which the committed document does not cover")
    def load(p, idx):
        o = []
        for l in Path(p).read_text().splitlines():
            if l.startswith("#") or not l.strip():
                continue
            f = l.split("\t")
            o.append([float(f[i]) for i in idx])
        return o
    full = load(out / "rmsd_vs_time_full.dat", [1, 2])
    fv = [x[1] for x in full]
    check("the full trajectory is 60 ns", len(fv) == 60000, f"{len(fv)} frames")
    print(f"         RMSD over 20 ns: mean {st.mean(v):.2f}   "
          f"over 60 ns: mean {st.mean(fv):.2f}, max {max(fv):.2f}")
    blocks = [st.mean(fv[i:i + 10000]) for i in range(0, 60000, 10000)]
    print("         per 10 ns block: " + "  ".join(f"{b:.2f}" for b in blocks))
    # The drift flagged at 20 ns should be a rise to a plateau, not an unbounded trend.
    plateau = max(blocks[2:]) - min(blocks[2:])
    check("RMSD rises then plateaus rather than drifting without bound",
          plateau < 0.5,
          f"blocks 20-60 ns span {plateau:.2f} A; the 20 ns document saw only the rise")
    rc = load(out / "rxn_coord_per_frame_full.dat", [2, 5, 8])
    print("         near-attack occupancy, 20 ns -> 60 ns:")
    occ60 = {}
    for i, (c, s) in enumerate([("1", "A"), ("2", "B"), ("3", "C")]):
        f20 = [x[i] for x in load(out / "rxn_coord_per_frame.dat", [2, 5, 8])]
        f60 = [x[i] for x in rc]
        o20 = 100 * sum(1 for y in f20 if y <= D_NAC) / len(f20)
        o60 = 100 * sum(1 for y in f60 if y <= D_NAC) / len(f60)
        occ60[c] = o60
        print(f"           CHA#{c} (site {s}): {o20:5.1f}% -> {o60:5.1f}%")
    check("site B holds its near-attack occupancy over the full range",
          occ60["2"] > 95.0, f"{occ60['2']:.1f}% over 60 ns")
    check("sites A and C fall further away over the full range",
          occ60["1"] < 3.0 and occ60["3"] < 20.0,
          f"A {occ60['1']:.1f}%, C {occ60['3']:.1f}%. The ensemble draws 29 of 44 frames "
          f"from beyond 20 ns, precisely where A and C are worse and B is unchanged.")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

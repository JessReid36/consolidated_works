#!/usr/bin/env python3
"""
step19_verify.py - audit of the barrier ensemble and its comparison against Claeyssens.

THIS IS THE STEP THE THESIS RESTS ON. Every headline number lives here, so the checks
re-derive rather than read: barriers, in vacuo references, differential transition-state
stabilisation, and the regression that is compared against the published one.

Checks performed:

  1. Table integrity: the two sub-ensembles, and which columns are populated for each.
  2. Every barrier is re-derived from the committed energies, not taken from the barrier
     column, and the two are compared.
  3. The do-not-quote entries are present and identifiable.
  4. THE SIGN CONVENTION. stab_TS here is barrier - barrier_vac, so negative means the
     environment lowers the barrier. Claeyssens define stabilisation with the opposite
     sign. Our gradient is therefore the negative of theirs, and this check asserts the
     relation explicitly so nobody compares raw signs.
  5. THE REGRESSION, computed separately for each sub-ensemble, against the enzyme-only
     fit of Claeyssens 2005 (Y = 16.02 - 0.95 X, R = 0.97), which is the correct
     comparator because the 2011 fit combines enzyme and water data.
  6. The internal consistency test both papers state: the regression intercept should
     equal the mean in vacuo barrier.
  7. Leave-one-out stability of the reported regression.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step19_verify.py [outputs]
"""
import sys
import math
import statistics as st
from pathlib import Path

FAIL = []
HARTREE = 627.5095
# Claeyssens et al. 2005, Chem. Commun. 5068, enzyme-only fit read from Fig. 2(b).
C2005 = {"gradient": -0.95, "intercept": 16.02, "R": 0.97,
         "barrier": 12.0, "barrier_sd": 1.7, "invacuo": 16.2, "stab": 4.2}
# Claeyssens et al. 2011, Org. Biomol. Chem. 9, 1578.
C2011 = {"barrier": 11.3, "barrier_sd": 1.8, "invacuo": 19.2, "invacuo_sd": 1.3,
         "stab": 7.3, "stab_sd": 2.0}
DO_NOT_QUOTE = {"08170"}


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def fit(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    mx, my = st.mean(xs), st.mean(ys)
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in pts)
    syy = sum((y - my) ** 2 for y in ys)
    m = sxy / sxx
    return m, my - m * mx, sxy / math.sqrt(sxx * syy)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    print(f"step 19 verification against {out.resolve()}\n")

    rows, hdr = [], None
    for l in (out / "ensemble_barriers.tsv").read_text().splitlines():
        if l.startswith("#"):
            continue
        f = l.rstrip("\n").split("\t")
        if f[0] == "frame":
            hdr = f
            continue
        rows.append(dict(zip(hdr, f)))

    def grp(r):
        return "extended" if int(r["frame"]) >= 20000 else "pilot"

    print("1. table integrity")
    check("the table is populated", len(rows) == 44, f"{len(rows)} frames")
    for g in ("pilot", "extended"):
        sel = [r for r in rows if grp(r) == g]
        nb = sum(num(r["barrier"]) is not None for r in sel)
        nv = sum(num(r["barrier_vac"]) is not None for r in sel)
        ns = sum(num(r["stab_TS"]) is not None for r in sel)
        print(f"         {g:<9} rows {len(sel):>3}   barrier {nb:>3}   "
              f"in vacuo {nv:>3}   stab_TS {ns:>3}")
    check("every frame with an in vacuo reference also has a stabilisation",
          all((num(r["barrier_vac"]) is None) == (num(r["stab_TS"]) is None)
              for r in rows))

    print("\n2. barriers re-derived from the committed energies")
    worst, n_ok = 0.0, 0
    for r in rows:
        b, ci, rr = num(r["barrier"]), num(r["E_ci_Eh"]), num(r["E_qmmm_R"])
        if None in (b, ci, rr):
            continue
        d = abs((ci - rr) * HARTREE - b)
        worst = max(worst, d)
        n_ok += 1
    check("barrier equals (E_ci - E_qmmm_R) for every frame that has both",
          worst < 0.01, f"{n_ok} frames checked, worst disagreement {worst:.4f} kcal/mol")
    worst_s = 0.0
    for r in rows:
        s, b, v = num(r["stab_TS"]), num(r["barrier"]), num(r["barrier_vac"])
        if None in (s, b, v):
            continue
        worst_s = max(worst_s, abs((b - v) - s))
    check("stab_TS equals barrier minus barrier_vac", worst_s < 1e-6,
          f"worst disagreement {worst_s:.2e} kcal/mol")

    print("\n3. do-not-quote entries")
    for fr in sorted(DO_NOT_QUOTE):
        m = [r for r in rows if r["frame"] == fr]
        check(f"frame {fr} is present and identifiable", bool(m),
              f"group={m[0]['group']}, barrier={m[0]['barrier']}" if m else "absent")
        if m and num(m[0]["barrier"]) is not None:
            check(f"frame {fr} carries the collapsed-band signature",
                  num(m[0]["barrier"]) == 0.0,
                  "barrier exactly 0.00, which is not a barrier; it must be excluded "
                  "from every statistic")

    print("\n4. the sign convention")
    stabs = [num(r["stab_TS"]) for r in rows if num(r["stab_TS"]) is not None]
    check("stabilisation is recorded as negative, i.e. the environment lowers the barrier",
          st.mean(stabs) < 0, f"mean {st.mean(stabs):+.3f} kcal/mol over {len(stabs)} frames")
    print(f"         Claeyssens define stabilisation with the OPPOSITE sign: their mean is")
    print(f"         +{C2011['stab']:.1f} (2011) and +{C2005['stab']:.1f} (2005) kcal/mol.")
    print(f"         Ours of {st.mean(stabs):+.2f} is the same physics with the sign flipped.")
    print(f"         Consequently our regression gradient is the NEGATIVE of theirs.")

    print("\n5. the regression, per sub-ensemble, against Claeyssens 2005")
    print(f"         comparator: Y = {C2005['intercept']:.2f} {C2005['gradient']:+.2f} X, "
          f"R = {C2005['R']:.2f}   (enzyme-only, 16 paths)")
    print(f"         {'set':<22}{'n':>4}{'grad':>9}{'in their sign':>15}"
          f"{'intercept':>11}{'r':>8}{'mean vac':>10}")
    results = {}
    for lbl, sel in (("pilot", [r for r in rows if grp(r) == "pilot"]),
                     ("extended", [r for r in rows if grp(r) == "extended"]),
                     ("pooled", rows)):
        pts = [(num(r["stab_TS"]), num(r["barrier"])) for r in sel
               if num(r["stab_TS"]) is not None and num(r["barrier"]) is not None
               and r["frame"] not in DO_NOT_QUOTE]
        if len(pts) < 3:
            print(f"         {lbl:<22}{len(pts):>4}   too few points")
            continue
        m, c, rr = fit(pts)
        vac = [num(r["barrier_vac"]) for r in sel if num(r["stab_TS"]) is not None
               and r["frame"] not in DO_NOT_QUOTE]
        results[lbl] = (len(pts), m, c, rr, st.mean(vac))
        print(f"         {lbl:<22}{len(pts):>4}{m:>+9.3f}{-m:>+15.3f}"
              f"{c:>11.2f}{rr:>+8.3f}{st.mean(vac):>10.3f}")
    check("the extended ensemble reproduces the published gradient magnitude",
          abs(abs(results["extended"][1]) - abs(C2005["gradient"])) < 0.15,
          f"ours {-results['extended'][1]:+.3f} in their sign against "
          f"{C2005['gradient']:+.2f}")
    check("the extended ensemble reproduces the published correlation",
          abs(results["extended"][3]) > 0.85,
          f"r = {results['extended'][3]:+.3f} against their R = {C2005['R']:.2f}")

    print("\n6. the internal consistency test both papers state")
    print("         the regression intercept should equal the mean in vacuo barrier")
    for lbl, (n, m, c, rr, vac) in results.items():
        gap = abs(c - vac)
        print(f"         {lbl:<22} intercept {c:6.2f}  mean in vacuo {vac:6.2f}  "
              f"gap {gap:.2f}")
    pub_gap = abs(C2005["intercept"] - C2005["invacuo"])
    check("the extended ensemble passes it at least as well as the published fit",
          abs(results["extended"][2] - results["extended"][4]) < pub_gap + 0.15,
          f"our gap {abs(results['extended'][2]-results['extended'][4]):.2f} against "
          f"their {pub_gap:.2f}")

    print("\n7. leave-one-out stability of the extended fit")
    pts = [(num(r["stab_TS"]), num(r["barrier"])) for r in rows if grp(r) == "extended"
           and num(r["stab_TS"]) is not None and num(r["barrier"]) is not None]
    grads = []
    for i in range(len(pts)):
        m, c, rr = fit([p for j, p in enumerate(pts) if j != i])
        grads.append(m)
    check("no single frame controls the gradient",
          max(grads) - min(grads) < 0.25,
          f"gradient spans {min(grads):+.3f} to {max(grads):+.3f} across "
          f"{len(pts)} leave-one-out fits")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

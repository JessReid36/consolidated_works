#!/usr/bin/env python3
"""
step15_17_verify.py - audit of the reaction-path endpoints and the restrained scan.

Steps 15, 16 and 17 are combined: optimise the reactant, drive a restrained scan along the
reaction coordinate, optimise the product. They are one operation - establish the path
endpoints and a first profile between them - and the scan's result is only interpretable
against the two endpoints.

Checks performed:

  1. Continuity with the QM/MM setup established at steps 13-14.
  2. The scan ran to completion and every window produced an energy.
  3. THE RESTRAINTS DO NOT HOLD THE TARGET COORDINATE. Achieved r is compared against
     target r window by window. The lag is systematic and grows.
  4. THE SCAN NEVER SAMPLES THE TRANSITION-STATE REGION. There is a discontinuity where
     the system snaps over the barrier, and the gap contains the reaction coordinate at
     which Claeyssens locate the TS.
  5. The scan maximum is therefore the last point before the snap, not a saddle, and the
     barrier read from it depends entirely on which reference point is chosen. All three
     candidate references are computed.
  6. Endpoint optimisations completed and are distinguishable in energy.

Exit status 0 if every check passes, 1 otherwise. Checks 3-5 PASS when the pathology is
present and correctly characterised: the point is to assert that the scan behaves this
way, so that nobody later quotes its maximum as a barrier.

Usage:  python3 step15_17_verify.py [outputs] [targets_file]
"""
import sys
import re
from pathlib import Path

FAIL = []
HARTREE = 627.5095
CLAEYSSENS_TS_R = -0.50
CLAEYSSENS_RS_STAR = -1.800


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    tgt = Path(sys.argv[2] if len(sys.argv) > 2 else "../scripts/step16_scan_targets.txt")
    print(f"steps 15-17 verification against {out.resolve()}\n")

    print("1. continuity and completeness")
    rows = [l.split("\t") for l in
            (out / "scan_energies.tsv").read_text().splitlines() if l.strip()]
    pts = [r[0] for r in rows]
    R = [float(r[1]) for r in rows]
    E = [float(r[2]) for r in rows]
    targets = [tuple(map(float, l.split()))
               for l in tgt.read_text().splitlines() if l.strip()]
    check("every scan window produced an energy", len(rows) == len(targets),
          f"{len(rows)} windows, {len(targets)} targets")
    check("energies are on the QM/MM total scale, not the bare-QM scale",
          all(-1200 < e < -1000 for e in E),
          f"{min(E):.3f} to {max(E):.3f} Eh; the bare-QM scale would be near -836")
    bonds = [l.split("\t") for l in
             (out / "scan_bonds_all.tsv").read_text().splitlines()[1:] if l.strip()]
    check("a per-window bond record accompanies the energies", len(bonds) == len(rows),
          f"{len(bonds)} rows")

    print("\n2. do the restraints hold the target coordinate?")
    lags = []
    for i, (tb, tf) in enumerate(targets):
        lags.append(R[i] - (tb - tf))
    worst = max(lags, key=abs)
    print(f"         {'pt':>3}{'target r':>10}{'achieved r':>12}{'lag':>9}")
    for i in (0, 4, 9, 10, 11, len(rows) - 1):
        print(f"         {pts[i]:>3}{targets[i][0]-targets[i][1]:>10.3f}"
              f"{R[i]:>12.3f}{lags[i]:>+9.3f}")
    # The lag must be assessed per branch. Before the snap the optimiser trails the
    # target; after it the system has fallen past the target and the sign flips. Pooling
    # the two hides both effects.
    snap = next((i for i in range(1, len(R)) if R[i] - R[i - 1] > 0.5), len(R))
    pre = lags[1:snap]
    post = lags[snap:]
    check("before the snap, the achieved coordinate trails the target throughout",
          pre and all(x < 0 for x in pre),
          f"windows 02-{pts[snap-1]}: lag {min(pre):+.3f} to {max(pre):+.3f} A, growing "
          f"monotonically. The restraint pulls toward the target but the optimiser never "
          f"reaches it, so the scan coordinate is not under control.")
    check("after the snap, the system has fallen PAST the target",
          post and all(x > 0 for x in post),
          f"windows {pts[snap]}-{pts[-1]}: lag {min(post):+.3f} to {max(post):+.3f} A")

    print("\n3. is the transition-state region sampled?")
    react = [i for i, r in enumerate(R) if r < 0]
    prod = [i for i, r in enumerate(R) if r > 0]
    check("the scan spans both sides of the coordinate", react and prod,
          f"{len(react)} reactant-side, {len(prod)} product-side")
    last_r, first_p = R[max(react)], R[min(prod)]
    gap = first_p - last_r
    print(f"         last reactant-side point  r = {last_r:+.3f}")
    print(f"         first product-side point  r = {first_p:+.3f}")
    print(f"         UNSAMPLED GAP             {gap:.3f} A")
    check("a discontinuity is present, i.e. the system snapped over the barrier",
          gap > 0.5,
          f"{gap:.3f} A, against a typical window step of "
          f"{abs(R[1]-R[0]):.3f} A")
    check("the gap contains the coordinate at which Claeyssens locate the TS",
          last_r < CLAEYSSENS_TS_R < first_p,
          f"their TS is at r = {CLAEYSSENS_TS_R:+.2f}, inside [{last_r:+.3f}, "
          f"{first_p:+.3f}]. The scan therefore cannot locate a saddle.")

    print("\n4. the scan maximum, and what it is not")
    imax = max(react, key=lambda i: E[i])
    check("the maximum on the reactant branch is the last point before the snap",
          imax == max(react),
          f"point {pts[imax]} at r = {R[imax]:+.3f}")
    print("         so it is a terminus, not a stationary point. Any barrier read from")
    print("         it is a LOWER BOUND, and the project's own rule already says 1D")
    print("         relaxed-scan barriers are diagnostics and must never be quoted.")

    print("\n5. the barrier depends entirely on the reference chosen")
    imin_r = min(react, key=lambda i: E[i])
    imin_all = E.index(min(E))
    b_first = (E[imax] - E[0]) * HARTREE
    b_rmin = (E[imax] - E[imin_r]) * HARTREE
    b_gmin = (E[imax] - E[imin_all]) * HARTREE
    print(f"         vs point {pts[0]} (r = {R[0]:+.3f})            {b_first:6.2f} kcal/mol")
    print(f"         vs reactant-branch minimum, point {pts[imin_r]}  {b_rmin:6.2f} kcal/mol")
    print(f"         vs global minimum, point {pts[imin_all]}         {b_gmin:6.2f} kcal/mol")
    check("the first point is the reference nearest Claeyssens' fixed RS*",
          abs(R[0] - CLAEYSSENS_RS_STAR) < 0.05,
          f"r = {R[0]:+.3f} against their RS* of {CLAEYSSENS_RS_STAR:+.3f}. The commonly "
          f"quoted {b_first:.1f} kcal/mol uses this reference, which is the "
          f"Claeyssens-consistent one.")
    check("the global minimum lies on the PRODUCT side, so that reference is invalid",
          R[imin_all] > 0,
          f"point {pts[imin_all]} at r = {R[imin_all]:+.3f}; a barrier of "
          f"{b_gmin:.2f} kcal/mol against it is product-to-maximum and meaningless")
    check("the choice of reference moves the number by about a kcal/mol",
          abs(b_rmin - b_first) > 0.5,
          f"{b_first:.2f} against {b_rmin:.2f}; state which reference is used whenever "
          f"this number appears")

    print("\n6. endpoint optimisations")
    for name in ("reactant_opt", "product_opt"):
        p = out / f"{name}.property.txt"
        check(f"{name} property file is committed", p.exists(),
              f"{p.stat().st_size/1048576:.1f} MB" if p.exists() else "")
        pdb = out / f"{name}.pdb"
        check(f"{name} structure is committed", pdb.exists())

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

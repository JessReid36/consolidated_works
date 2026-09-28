#!/usr/bin/env python3
"""
theta2_autocorr.py - correlation time of the near-attack angle theta2, the only
observable shown here to predict the barrier.

WHY THIS AND NOT THE DISTANCE
    step11e_autocorr.py measures the correlation time of the forming C1-C6
    distance and the Arg90-O13 contact, describing them as "the two quantities
    that define catalytic competence", and finds 3.2 and 2.5 ps. On that basis
    frames spaced 855 to 1680 ps apart are effectively independent.

    Across the thirteen frames for which an in vacuo barrier has been computed,
    neither of those quantities predicts it:

        theta2 deviation                 r = +0.657
        forming C1-C6 at selection       r = -0.202
        theta1 deviation                 r = +0.076
        Arg90-O13 at selection           r = -0.409

    Grossfield et al. (2018) define the correlation time as a property of an
    observable, "In time-series data of a random quantity x(t)", and advise
    comparing "auto-correlation times for individual observables". The one that
    governs sampling is the slowest among those the result depends on, so the
    independence claim needs the correlation time of theta2, not of a distance
    that does not track the barrier.

WHAT IT DOES
    Reads both trajectory files as one continuous run, computes theta2 for the
    substrate at each frame using nac_angles() imported from
    step11f_nac_check.py so the angle definition is the published one rather
    than a reimplementation, and reports:

      - the integrated correlation time by Sokal's automatic windowing, the same
        estimator step11e uses, so the two numbers are comparable;
      - the same for theta1 and the forming distance, as a control: if those
        reproduce step11e's values the difference is the observable and not the
        method;
      - the effective sample size, and the spacing required for independence;
      - Chodera's equilibration scan on theta2, choosing the discard point that
        maximises the effective sample size rather than reading it off block
        means by eye, which is how the present 11 ns cut was chosen.

    A stride is available since theta2 need not be evaluated every picosecond to
    resolve a correlation time of tens or hundreds of picoseconds; the default
    reads every frame so the fast end is not missed.

Usage:
    python3 theta2_autocorr.py [--stride 1] [--max-lag 8000]
Run under the batch system: it reads 60000 frames and uses numpy.
"""
import argparse
import importlib.util
import os
import struct
import sys

import numpy as np

HOME = os.path.expanduser("~")
NAC = f"{HOME}/system_development/phase1_system_dev/step11f_nac_check.py"


def load_nac():
    """Import the committed near-attack module, so the angle definition and the
    element assignment are the published ones and not a second implementation."""
    if not os.path.exists(NAC):
        sys.exit(f"FAIL: {NAC} not found; copy it to the cluster first")
    spec = importlib.util.spec_from_file_location("nac", NAC)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for fn in ("nac_angles",):
        if not hasattr(m, fn):
            sys.exit(f"FAIL: {NAC} has no {fn}()")
    return m


def nc_open(path):
    """Record layout of a netCDF classic file. Same reader as the committed
    analysis scripts."""
    fsz = os.path.getsize(path)
    data = open(path, "rb").read(1 << 16)
    assert data[:3] == b"CDF"
    off64 = data[3] == 2
    pos = [4]

    def u32():
        v = struct.unpack_from(">I", data, pos[0])[0]; pos[0] += 4; return v

    def off():
        v = struct.unpack_from(">Q" if off64 else ">I", data, pos[0])[0]
        pos[0] += 8 if off64 else 4
        return v

    def nm():
        n = u32(); s = data[pos[0]:pos[0]+n].decode("ascii", "replace")
        pos[0] += n + ((4 - (n % 4)) % 4); return s

    numrecs = u32()
    tag = u32(); nd = u32(); dims = []
    if tag == 0x0A:
        for _ in range(nd):
            dims.append((nm(), u32()))

    def skip_att():
        t = u32(); na = u32()
        if t == 0x0C:
            for _ in range(na):
                nm(); tp = u32(); k = u32()
                nb = k * {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}.get(tp, 1)
                pos[0] += nb + ((4 - (nb % 4)) % 4)

    skip_att(); rec = []; tag = u32(); nv = u32()
    if tag == 0x0B:
        for _ in range(nv):
            vn = nm(); ndv = u32(); dimids = [u32() for _ in range(ndv)]
            skip_att(); u32(); vs = u32(); bg = off()
            if ndv > 0 and dims[dimids[0]][1] == 0:
                rec.append({"name": vn, "vsize": vs, "begin": bg})
    dimd = {n: l for n, l in dims}
    recsize = sum(v["vsize"] for v in rec)
    frames = numrecs if numrecs != 0xFFFFFFFF else \
        (fsz - min(v["begin"] for v in rec)) // recsize
    return {"frames": frames, "recsize": recsize, "rec": rec,
            "atom": dimd.get("atom")}


def integrated_tau(x, max_lag):
    """Sokal's automatic windowing, c = 5. The same estimator step11e uses, so
    the numbers are directly comparable."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    var = np.dot(x, x) / n
    if var <= 0:
        return float("nan"), float("nan")
    size = 1
    while size < 2 * n:
        size *= 2
    f = np.fft.rfft(x, size)
    acf = np.fft.irfft(f * np.conjugate(f), size)[:n].real
    acf /= acf[0]
    tau = 1.0
    window = min(max_lag, n - 1)
    for k in range(1, window + 1):
        tau += 2.0 * acf[k]
        if k >= 5 * tau:
            break
    ess = n / (2.0 * tau) if tau > 0 else float("nan")
    return tau, ess


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prmtop",
                    default=f"{HOME}/system_development/03_amber/tleap_build/"
                            f"complex_solvated.prmtop")
    ap.add_argument("--traj",
                    default=f"{HOME}/system_development/04_amber_md/"
                            f"10c_production/prod.nc")
    ap.add_argument("--traj2",
                    default=f"{HOME}/system_development/04_amber_md/"
                            f"10d_production_extend/prod_ext.nc")
    ap.add_argument("--qm-lo", type=int, default=6207)
    ap.add_argument("--qm-hi", type=int, default=6231)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-lag", type=int, default=8000)
    ap.add_argument("--out",
                    default=f"{HOME}/system_development/04_amber_md/"
                            f"11_analysis/theta_per_frame.dat")
    args = ap.parse_args()

    nac = load_nac()
    print(f"angle definition imported from {os.path.basename(NAC)}")

    # element symbols for the substrate, from the topology
    els = []
    cur = None
    for ln in open(args.prmtop, errors="replace"):
        if ln.startswith("%FLAG"):
            cur = ln.split()[1]; continue
        if ln.startswith("%"):
            continue
        if cur == "ATOM_NAME":
            for i in range(0, len(ln.rstrip()), 4):
                tok = ln[i:i+4].strip()
                if tok:
                    els.append(tok)
    els = els[args.qm_lo:args.qm_hi]
    els = [(e[0] if e[0] in "CHNO" else "C") for e in els]
    print(f"substrate: {len(els)} atoms, {''.join(els)}")

    H1 = nc_open(args.traj)
    H2 = nc_open(args.traj2)
    c1 = next(v for v in H1["rec"] if v["name"] == "coordinates")
    c2 = next(v for v in H2["rec"] if v["name"] == "coordinates")
    N1, N2 = H1["frames"], H2["frames"]
    natom = H1["atom"]
    NT = N1 + N2
    print(f"trajectories: {N1} + {N2} = {NT} frames, one per picosecond")

    def sub(i):
        if i < N1:
            path, hdr, cv, j = args.traj, H1, c1, i
        else:
            path, hdr, cv, j = args.traj2, H2, c2, i - N1
        o = cv["begin"] + j * hdr["recsize"]
        a = np.memmap(path, dtype=">f4", mode="r", offset=o, shape=(natom, 3))
        return np.array(a[args.qm_lo:args.qm_hi], float)

    idx = list(range(0, NT, args.stride))
    d_s, t1_s, t2_s = [], [], []
    bad = 0
    for k, i in enumerate(idx):
        r = nac.nac_angles(els, sub(i))
        if r is None:
            bad += 1
            continue
        d, t1, t2 = r
        d_s.append(d); t1_s.append(t1); t2_s.append(t2)
        if k and k % 10000 == 0:
            print(f"  {k}/{len(idx)}")
    if bad:
        print(f"  {bad} frames gave no angle and were skipped")
    print(f"  {len(t2_s)} frames measured\n")

    with open(args.out, "w") as fh:
        fh.write("# frame\td_form\ttheta1\ttheta2\n")
        for i, d, a, b in zip(idx, d_s, t1_s, t2_s):
            fh.write(f"{i}\t{d:.3f}\t{a:.2f}\t{b:.2f}\n")
    print(f"wrote {args.out}\n")

    dt = args.stride
    print("=" * 78)
    print("Correlation times, Sokal windowing, the same estimator as step11e")
    print("=" * 78)
    print(f"{'observable':<26}{'tau (frames)':>14}{'tau (ps)':>12}"
          f"{'eff. N':>10}{'r with barrier':>16}")
    print("-" * 78)
    rvals = {"theta2": "+0.657", "theta1": "+0.076", "d_form": "-0.202"}
    taus = {}
    for lbl, series in (("theta2 (predicts barrier)", t2_s),
                        ("theta1", t1_s),
                        ("forming C1-C6", d_s)):
        tau, ess = integrated_tau(series, args.max_lag)
        key = "theta2" if "theta2" in lbl else ("theta1" if "theta1" in lbl else "d_form")
        taus[lbl] = tau * dt
        print(f"{lbl:<26}{tau:>14.1f}{tau*dt:>12.1f}{ess:>10.0f}{rvals[key]:>16}")
    print()
    print("  The r column is the correlation of each quantity with the in vacuo")
    print("  barrier across the thirteen frames for which one has been computed.")
    print("  Only theta2 predicts it, so its correlation time is the one that")
    print("  governs whether frames are independent for this purpose.")

    t2_tau = taus["theta2 (predicts barrier)"]
    print()
    print("=" * 78)
    print("What spacing does theta2 require?")
    print("=" * 78)
    print(f"  tau(theta2) = {t2_tau:.1f} ps")
    print(f"  independence at 5 tau: {5*t2_tau:.0f} ps between frames")
    print(f"  available equilibrated trajectory, from 11000 ps: {NT-11000} ps")
    print(f"  independent configurations at that spacing: "
          f"{(NT-11000)/(5*t2_tau):.0f}")
    print()
    print("  step11e reported 3.2 ps for the forming distance and concluded the")
    print("  spacing was 269 times the correlation time. Compare that with the")
    print("  figure above, which is for the observable that tracks the barrier.")

    # ------------------------------------- Chodera's equilibration detection
    print()
    print("=" * 78)
    print("Where does equilibration end, by Chodera's criterion, for theta2?")
    print("=" * 78)
    print("  Scan the discard point and take the value leaving the most")
    print("  effectively uncorrelated samples. The present 11000 ps cut was read")
    print("  off backbone-deviation block means by eye.")
    print()
    print(f"  {'t0 / ps':>9}{'tau / ps':>11}{'N_eff':>10}")
    rows = []
    for t0 in range(0, min(NT // 2, 30000), 2000):
        seg = t2_s[t0 // dt:]
        if len(seg) < 2000:
            break
        tau, ess = integrated_tau(seg, args.max_lag)
        if tau != tau or tau <= 0:
            continue
        rows.append((t0, tau * dt, ess))
    if rows:
        best = max(rows, key=lambda r: r[2])
        for t0, tau, ess in rows:
            mark = "   <- maximum" if (t0, tau, ess) == best else ""
            print(f"  {t0:>9}{tau:>11.1f}{ess:>10.0f}{mark}")
        print()
        print(f"  Chodera's choice: discard {best[0]} ps. The present cut is 11000.")
        if abs(best[0] - 11000) > 4000:
            print("  These differ materially. The thirty-frame selection was made")
            print("  from the 11000 ps cut and should be reconsidered.")
        else:
            print("  These agree closely enough that the existing selection stands.")
        if best[0] > 0.4 * NT:
            print()
            print("  WARNING: the chosen discard covers more than two fifths of the")
            print("  record. Grossfield et al. note this can happen when a run is")
            print("  too short to sample transitions between metastable states, in")
            print("  which case the method 'can simply result in restricting the")
            print("  production region to the last sampled metastable basin'.")


if __name__ == "__main__":
    main()

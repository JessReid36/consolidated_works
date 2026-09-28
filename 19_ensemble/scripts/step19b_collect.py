#!/usr/bin/env python3
"""
step19b_collect.py - collect the ensemble QM/MM barriers and report the
distribution, so the Phase-1 barrier can be quoted as a mean with a spread
rather than as a single path.

Reads the converged NEB climbing-image barrier from each frame directory under
05_qmmm/19_ensemble/, joins it to the step-12a selection manifest, and reports:

  * mean, standard deviation, and range of the barrier across frames
  * the standard error on the mean, which is what a comparison with a literature
    average or with experiment should be quoted against
  * where the fully characterised frame (820) sits within the distribution
  * correlation of the barrier with the Arg90-O13 contact distance and with the
    near-attack C6-C1 distance recorded at selection time

The last of these is a real test rather than a decoration. The literature
attributes transition-state stabilisation in BsCM largely to Arg90, so if the
barrier varies across frames mainly through electrostatic contact, a negative
correlation between barrier and Arg90 proximity is the expected signature.
Finding it would explain the frame-to-frame spread mechanistically; not finding
it would say the spread comes from somewhere else and is worth understanding
before the mean is quoted.

BARRIER SOURCE
--------------
The converged NEB climbing image, read from the ORCA output. This is the
validated proxy: on frame 820 it gave 15.94 kcal/mol against the fully
characterised reduced-region value of 16.00, agreeing to 0.06. Only converged
bands are counted; an unconverged NEB is reported and excluded rather than
silently averaged in.

Usage
-----
  python3 step19b_collect.py
  python3 step19b_collect.py --ensdir /path/to/19_ensemble --manifest /path/to/selection_manifest.tsv
"""

import argparse
import os
import re
import sys
from pathlib import Path

HARTREE2KCAL = 627.5094740631
REFERENCE_FRAME = 820
REFERENCE_FULL = 16.00     # fully characterised reduced-region barrier, frame 820
REFERENCE_NEB = 15.94      # its NEB-CI proxy, for the validation statement
EXPERIMENT_DH = 12.7       # experimental activation enthalpy, kcal/mol


def read_manifest(path):
    rows = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        rows[int(f[1])] = {
            "ps": float(f[2]),
            "form": float(f[5]),
            "r": float(f[6]),
            "arg90": float(f[7]),
        }
    return rows


def neb_barrier(out_path):
    """
    Converged climbing-image barrier, in kcal/mol.

    Returns (value, 'ok') or (None, reason). The energies are taken from the
    frozen/active energy table that ORCA prints for the converged band; the
    trajectory comment lines are NOT used, because Knarr writes a
    spring-augmented internal energy there that is not the physical barrier.
    """
    if not out_path.exists():
        return None, "no output"
    txt = out_path.read_text(errors="replace")

    if "ORCA TERMINATED NORMALLY" not in txt:
        # An unfinished job has no termination line either, so distinguish the
        # two rather than calling a running job a failure.
        if re.search(r"OPTIMIZATION CYCLE|LBFGS", txt):
            return None, "still running"
        return None, "did not terminate normally"
    if not re.search(r"THE NEB OPTIMIZATION HAS CONVERGED", txt, re.I):
        return None, "band not converged"

    # last block of image energies: lines ending in '@' carry the HEI column
    # ORCA 6.0.1 marks the climbing image with "<= CI" at the end of its row in
    # the PATH SUMMARY table, whose columns are:
    #   Image  Dist.(Ang.)  E(Eh)  dE(kcal/mol)  max(|Fp|)  RMS(Fp)
    # The barrier is the dE of the CI row. Older ORCA marked that row with "@",
    # which is what this pattern used to look for and which 6.0.1 never writes.
    hits = re.findall(r"^\s*\d+\s+[-\d.]+\s+[-\d.]+\s+([-\d.]+)\s+\S+\s+\S+\s*<=\s*CI\s*$",
                      txt, re.M)
    if not hits:
        # fall back to the maximum dE in the PATH SUMMARY table, which is the
        # same number whenever the climbing image sits at the peak
        blk = re.search(r"PATH SUMMARY.*?(?=Straight line distance)", txt, re.S)
        if blk:
            rows = re.findall(r"^\s*\d+\s+[-\d.]+\s+[-\d.]+\s+([-\d.]+)\s",
                              blk.group(0), re.M)
            if rows:
                hits = [max(rows, key=float)]
    if hits:
        try:
            return float(hits[-1]), "ok"
        except ValueError:
            pass

    # fall back to the explicit climbing-image statement
    m = re.findall(r"climbing image.*?([\d.]+)\s*kcal", txt, re.I | re.S)
    if m:
        return float(m[-1]), "ok"

    return None, "converged but no barrier found in output"


def read_harvest(path):
    """Optimised reactant geometry per frame, from step19d_harvest.py.

    These are the distances in the structure each barrier was actually computed
    from, as distinct from the manifest values, which describe the molecular
    dynamics snapshot the optimisation started from.
    """
    out = {}
    p = Path(path)
    if not p.exists():
        return out
    lines = [l.split("\t") for l in p.read_text().splitlines() if l.strip()]
    hdr = [h.strip() for h in lines[0]]
    try:
        i_f = hdr.index("frame")
        i_b = hdr.index("break_O3_C4")
        i_m = hdr.index("form_C1_C6")
    except ValueError:
        return out
    for r in lines[1:]:
        try:
            out[int(r[i_f])] = {"opt_break": float(r[i_b]), "opt_form": float(r[i_m])}
        except (ValueError, IndexError):
            continue
    return out


def stats(vals):
    n = len(vals)
    mean = sum(vals) / n
    if n > 1:
        sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5
        sem = sd / n ** 0.5
    else:
        sd = sem = float("nan")
    return mean, sd, sem


def pearson(x, y):
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / (sxx * syy) ** 0.5


def main():
    ap = argparse.ArgumentParser()
    home = os.path.expanduser("~")
    ap.add_argument("--ensdir", default=f"{home}/system_development/05_qmmm/19_ensemble")
    ap.add_argument("--manifest",
                    default=f"{home}/system_development/05_qmmm/12_frame_selection/"
                            f"selection_manifest.tsv")
    ap.add_argument("--harvest",
                    default=f"{home}/system_development/05_qmmm/19_ensemble_harvest/"
                            f"reactant_summary.tsv",
                    help="optimised reactant geometries, from step19d_harvest.py")
    ap.add_argument("--out", default="ensemble_barriers.tsv")
    args = ap.parse_args()

    man = read_manifest(args.manifest)
    harv = read_harvest(args.harvest)
    ensdir = Path(args.ensdir)

    results, skipped = [], []

    # frame 820 is carried in from its own directory, by its NEB value, so that
    # every entry in the distribution is the same kind of measurement
    results.append((REFERENCE_FRAME, REFERENCE_NEB, man.get(REFERENCE_FRAME, {})))

    for d in sorted(ensdir.glob("frame_*")):
        try:
            fr = int(d.name.split("_")[1])
        except (IndexError, ValueError):
            continue
        if fr == REFERENCE_FRAME:
            continue
        val, why = neb_barrier(d / "neb.out")
        if val is None:
            skipped.append((fr, why))
        else:
            results.append((fr, val, man.get(fr, {})))

    results.sort()

    print("  as selected: distances measured on the molecular dynamics snapshot")
    print("  as optimised: distances in the reactant the barrier was computed from")
    print()
    print(f"{'frame':>7}{'barrier':>9}{'Arg90 sel':>11}{'form sel':>10}"
          f"{'break opt':>11}{'form opt':>10}")
    print("-" * 58)
    for fr, val, m in results:
        h = harv.get(fr, {})
        tag = "  <- fully characterised" if fr == REFERENCE_FRAME else ""
        print(f"{fr:>7}{val:>9.2f}"
              f"{m.get('arg90', float('nan')):>11.3f}"
              f"{m.get('form', float('nan')):>10.3f}"
              f"{h.get('opt_break', float('nan')):>11.3f}"
              f"{h.get('opt_form', float('nan')):>10.3f}{tag}")

    if skipped:
        print("\nnot included:")
        for fr, why in skipped:
            print(f"  frame {fr}: {why}")

    vals = [v for _, v, _ in results]
    MIN_FOR_STATS = 6
    if len(vals) < MIN_FOR_STATS:
        print(f"\n{len(vals)} barrier(s) available. Summary statistics are withheld")
        print(f"below {MIN_FOR_STATS} frames: a mean, standard deviation and standard")
        print("error computed from a handful of frames invite being read as the")
        print("ensemble result, and the frames that finish first are not a random")
        print("sample of the ensemble.")
        if len(vals) >= 2:
            print(f"\n  observed so far: {min(vals):.2f} to {max(vals):.2f} kcal/mol, "
                  f"spread {max(vals)-min(vals):.2f}")
        sys.exit(0)

    mean, sd, sem = stats(vals)
    print(f"\nENSEMBLE  (n = {len(vals)})")
    print(f"  mean          {mean:.2f} kcal/mol")
    print(f"  sd            {sd:.2f}")
    print(f"  sem           {sem:.2f}   <- quote comparisons against this")
    print(f"  range         {min(vals):.2f} to {max(vals):.2f} "
          f"(spread {max(vals)-min(vals):.2f})")

    ref = [v for f, v, _ in results if f == REFERENCE_FRAME][0]
    below = sum(1 for v in vals if v < ref)
    print(f"\n  frame {REFERENCE_FRAME}: {ref:.2f} by NEB-CI, "
          f"{REFERENCE_FULL:.2f} fully characterised")
    print(f"    sits above {below} of {len(vals)-1} other frames; "
          f"{(ref-mean)/sd:+.2f} sd from the ensemble mean" if sd == sd else "")
    print(f"    proxy validation on this frame: "
          f"{abs(REFERENCE_FULL-REFERENCE_NEB):.2f} kcal/mol")

    print(f"\n  vs experimental dH‡ {EXPERIMENT_DH}: "
          f"ensemble mean is {mean-EXPERIMENT_DH:+.2f} kcal/mol")
    print(f"    (a potential-energy barrier is not dH‡; zero-point and thermal")
    print(f"     corrections are needed before this is a like-for-like comparison)")

    MIN_FOR_CORR = 8
    have = [(fr, v, m) for fr, v, m in results if m]
    if len(have) < MIN_FOR_CORR:
        print(f"\n  Correlations are withheld below {MIN_FOR_CORR} frames.")
    else:
        b = [v for _, v, _ in have]
        print("\n  CORRELATIONS across the ensemble")
        print("  Measured on the OPTIMISED reactant, the structure each barrier was")
        print("  computed from. This is the mechanistic question.")
        for key, label, expect in [
            ("opt_form", "forming C1-C6", "positive: further from attack, higher barrier"),
            ("opt_break", "breaking O3-C4", "no expectation; this bond is nearly invariant"),
        ]:
            pairs = [(harv[fr][key], v) for fr, v, _ in have
                     if fr in harv and key in harv[fr]]
            if len(pairs) >= MIN_FOR_CORR:
                r = pearson([a for a, _ in pairs], [c for _, c in pairs])
                print(f"    barrier vs {label:<18} r = {r:+.3f}  (n = {len(pairs)})"
                      f"   {expect}")
            else:
                print(f"    barrier vs {label:<18} only {len(pairs)} frames harvested")
        print("\n  Measured AS SELECTED, on the molecular dynamics snapshot. This asks")
        print("  whether the sampled conformation predicts the barrier, which bears on")
        print("  whether the frame-selection criteria were informative.")
        for key, label, expect in [
            ("arg90", "Arg90-O13 contact", "positive: looser contact, higher barrier"),
            ("form", "near-attack C6-C1", "positive: further from attack, higher barrier"),
        ]:
            x = [m[key] for _, _, m in have]
            r = pearson(x, b)
            print(f"    barrier vs {label:<18} r = {r:+.3f}  (n = {len(x)})   {expect}")
        print("\n    A clear positive correlation with the Arg90 contact would say the")
        print("    frame-to-frame spread is electrostatic in origin, consistent with the")
        print("    literature attribution of TS stabilisation to that residue. Absence of")
        print("    one says the spread comes from elsewhere and is worth understanding")
        print("    before the mean is quoted.")

    with open(args.out, "w") as fh:
        fh.write("frame\tprod_ps\tbarrier_kcal\tArg90_O13\tform_C6_C1\tr\n")
        for fr, val, m in results:
            fh.write(f"{fr}\t{m.get('ps','')}\t{val:.3f}\t{m.get('arg90','')}"
                     f"\t{m.get('form','')}\t{m.get('r','')}\n")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()

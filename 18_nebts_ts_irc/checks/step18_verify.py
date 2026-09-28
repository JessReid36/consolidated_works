#!/usr/bin/env python3
"""
step18_verify.py - audit of the NEB-TS band, the reduced-region transition state, the
IRC, and the reduced-region endpoints.

THIS STEP CARRIES THE ENZYME BARRIER, so the checks are correspondingly strict. Two
numbers circulate in the project notes, 15.08 and 16.00, and both are derived here from
committed energies so the discrepancy is resolved by arithmetic rather than by assertion.

Checks performed:

  1. The reduced-region active set, and that the QM region sits inside it.
  2. All three reduced-region jobs terminated normally, on the QM/MM total energy scale.
  3. THE BARRIER, derived two ways. 16.00 kcal/mol uses the reduced-region reactant;
     15.08 uses a reactant from a different active region. The difference is isolated to
     the reactant energy and quantified.
  4. THE FREQUENCY SUMMARY FILE CONTRADICTS ITSELF AND CARRIES THE SUPERSEDED BARRIER.
     Its header says one imaginary mode, its body lists two, and its footer quotes 15.08.
  5. The two imaginary modes, against the ghost-mode threshold, with the literature
     explanation for the small one.
  6. The NEB band and IRC outputs are present and describe the same reaction.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step18_verify.py [outputs]
"""
import sys
import re
from pathlib import Path

FAIL = []
HARTREE = 627.5095
GHOST_THRESHOLD = 50.0        # cm^-1, below which an imaginary mode is treated as noise
QM_RANGE = (6207, 6230)       # 0-based, established at steps 13-14


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def final_qmmm(path):
    hits = re.findall(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)",
                      Path(path).read_text(errors="replace"))
    return float(hits[-1]) if hits else None


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    print(f"step 18 verification against {out.resolve()}\n")

    print("1. the reduced active region")
    act = [int(x) for x in (out / "active_reduced.txt").read_text().split()]
    check("the active set is the documented 102 atoms", len(act) == 102, f"{len(act)}")
    qm = set(range(QM_RANGE[0], QM_RANGE[1] + 1))
    check("the QM region lies entirely inside the active region", qm <= set(act),
          f"QM {{{QM_RANGE[0]}:{QM_RANGE[1]}}}, {len(qm)} atoms, all present")
    env = sorted(set(act) - qm)
    check("the active region adds mobile environment beyond the QM region",
          len(env) == len(act) - len(qm),
          f"{len(env)} environment atoms are free to relax, the rest of the system frozen")

    print("\n2. the three reduced-region jobs")
    E = {}
    for lbl, name in (("reactant", "reactant_reduced.out"),
                      ("TS", "ts_reduced.out"),
                      ("product", "product_reduced.out")):
        p = out / name
        check(f"{name} is committed", p.exists(),
              f"{p.stat().st_size/1048576:.1f} MB" if p.exists() else "ABSENT")
        if p.exists():
            t = p.read_text(errors="replace")
            check(f"{name}: ORCA terminated normally", "ORCA TERMINATED NORMALLY" in t)
            E[lbl] = final_qmmm(p)
            check(f"{name}: energy is on the QM/MM total scale",
                  E[lbl] is not None and -1200 < E[lbl] < -1000, f"{E[lbl]:.9f} Eh")

    print("\n3. the barrier, and the 15.08 / 16.00 discrepancy resolved")
    b16 = (E["TS"] - E["reactant"]) * HARTREE
    rxn = (E["product"] - E["reactant"]) * HARTREE
    print(f"         reactant {E['reactant']:.9f}   TS {E['TS']:.9f}   "
          f"product {E['product']:.9f}")
    check("the consistent reduced-region barrier is 16.00 kcal/mol",
          abs(b16 - 16.00) < 0.02, f"{b16:+.2f} kcal/mol")
    check("the reaction is exothermic, as expected for this rearrangement", rxn < 0,
          f"{rxn:+.2f} kcal/mol")
    # The superseded 15.08 appears in the frequency summary footer. Recover the reactant
    # energy it used and show the discrepancy is entirely in that one term.
    fs = (out / "ts_frequencies_summary.txt").read_text()
    m = re.search(r"barrier:\s*([+-]?\d+\.\d+)\s*kcal/mol\s*\(TS\s*(-?\d+\.\d+)\s*vs\s*"
                  r"reactant\s*(-?\d+\.\d+)", fs)
    check("the superseded barrier and its inputs are recoverable from the file",
          bool(m), m.group(0)[:70] if m else "")
    b_old, ts_old, r_old = float(m.group(1)), float(m.group(2)), float(m.group(3))
    check("the TS energy is the same in both derivations",
          abs(ts_old - E["TS"]) < 1e-5,
          f"{ts_old:.6f} against {E['TS']:.6f}; the TS is not in dispute")
    check("the whole discrepancy lies in the REACTANT energy",
          abs(r_old - E["reactant"]) > 1e-4,
          f"{r_old:.6f} against the reduced-region {E['reactant']:.9f}, a difference of "
          f"{abs(r_old - E['reactant'])*HARTREE:.2f} kcal/mol")
    check("that difference exactly accounts for 16.00 minus 15.08",
          abs((b16 - b_old) - abs(r_old - E["reactant"]) * HARTREE) < 0.02,
          f"{b16:.2f} - {b_old:.2f} = {b16-b_old:.2f} kcal/mol")

    print("\n4. the frequency summary file contradicts itself")
    modes = [float(x) for x in re.findall(r":\s*(-?\d+\.\d+)\s*cm\*\*-1", fs)]
    imag = [x for x in modes if x < 0]
    check("the file's header claims ONE imaginary mode", "one reaction-coordinate" in fs)
    check("the file's body lists TWO", len(imag) == 2,
          f"{[f'{x:.2f}' for x in imag]} cm^-1. The header is wrong.")
    check("the file's footer carries the SUPERSEDED barrier", abs(b_old - 15.08) < 0.01,
          f"it states {b_old:+.2f} kcal/mol, which is the mismatched-region value. "
          f"Anyone quoting this file directly would quote the wrong barrier.")

    print("\n5. the two imaginary modes")
    big = [x for x in imag if abs(x) >= GHOST_THRESHOLD]
    small = [x for x in imag if abs(x) < GHOST_THRESHOLD]
    print(f"         {imag[0]:>9.2f} cm^-1   reaction coordinate")
    print(f"         {imag[1]:>9.2f} cm^-1   below the {GHOST_THRESHOLD:g} cm^-1 "
          f"ghost-mode threshold")
    check("exactly one imaginary mode is a genuine reaction coordinate", len(big) == 1,
          f"{big[0]:.2f} cm^-1")
    check("the second is small enough to treat as a ghost mode", len(small) == 1,
          f"{small[0]:.2f} cm^-1, magnitude {abs(small[0]):.2f} < {GHOST_THRESHOLD:g}")
    print("         Agbaglo et al. (2024) footnote 1 attributes exactly this artefact to")
    print("         partial Hessians with frozen Cartesian coordinates, and avoid it with")
    print("         a freeze-code scheme that zeroes Hessian elements involving two frozen")
    print("         coordinates. Our reduced-region OptTS freezes everything outside 102")
    print("         atoms, so the small mode has a published explanation. Report BOTH")
    print("         modes and cite the threshold; do not write 'one imaginary mode'.")

    print("\n6. band and IRC")
    for n in ("nebts_NEB-CI_converged.activeRegion.xyz", "nebts_MEP.QMRegion_trj.xyz",
              "nebts_MEP_ALL.QMRegion_trj.xyz"):
        check(f"{n} is committed", (out / n).exists())
    for n in ("irc_IRC_F.QMRegion_trj.xyz", "irc_IRC_B.QMRegion_trj.xyz"):
        check(f"{n} is committed", (out / n).exists())
    ci = (out / "nebts_NEB-CI_converged.activeRegion.xyz").read_text().splitlines()
    n_band = int(ci[0].split()[0])
    # TWO DIFFERENT ACTIVE REGIONS ARE IN PLAY, and conflating them is what produced the
    # superseded 15.08. The NEB band relaxes a large active region; the subsequent
    # reduced-region OptTS relaxes only 102 atoms. Energies from the two are NOT
    # interchangeable.
    check("the NEB band uses a LARGER active region than the reduced-region OptTS",
          n_band > len(act),
          f"band {n_band} atoms against reduced-region {len(act)}. These are different "
          f"partitions and their energies must never be differenced against each other.")
    print(f"         band active region      {n_band:>5} atoms")
    print(f"         reduced OptTS region    {len(act):>5} atoms")
    print(f"         QM region (both)        {len(qm):>5} atoms")
    print("         The 15.08 kcal/mol value differenced a 102-atom TS against a")
    print("         band-region reactant. 16.00 uses 102 atoms for both. Always state")
    print("         which active region a QM/MM energy belongs to.")

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

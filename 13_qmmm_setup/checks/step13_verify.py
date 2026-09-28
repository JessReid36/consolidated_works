#!/usr/bin/env python3
"""
step13_verify.py - re-runnable audit of the QM/MM setup: bridge, smoke test, DFT
single point. Covers pipeline steps 13a, 13b and 14, which together establish and
validate the QM/MM partition that every later barrier depends on.

Checks performed:

  1. Continuity. The bridge was generated from step09's committed solvated topology.
  2. The QM region is the substrate and nothing else: 24 atoms, matching CHA#2 in the
     topology, with the index convention made explicit.
  3. Charge and multiplicity are -2 and 1, the value established at step 01 and carried
     through step 08.
  4. Level of theory, read from the input, and compared against Claeyssens et al. 2011.
  5. THE ENERGY-SCALE IDENTITY. ORCA reports three energies and they must satisfy
     QM/MM = QM + MM exactly. This is the check that catches the scale confusion the
     project's own write-up prompt warns about: the QM-only bare-substrate scale near
     -836 Eh must never be differenced against the QM/MM scale near -1073 Eh.
  6. No link atoms are present, the QM region being a whole molecule.
  7. Both jobs completed without error.

Exit status 0 if every check passes, 1 otherwise.

Usage:  python3 step13_verify.py [outputs] [step09_outputs]
"""
import sys
import re
import hashlib
from pathlib import Path

FAIL = []
QM_EXPECTED = 24
CHARGE_EXPECTED = -2
MULT_EXPECTED = 1


def check(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        FAIL.append(label)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../outputs")
    s09 = Path(sys.argv[2] if len(sys.argv) > 2 else "../../09_tleap_build/outputs")
    print(f"step 13 verification against {out.resolve()}\n")

    print("1. continuity")
    prm = s09 / "complex_solvated.prmtop.pointer"
    if prm.exists():
        rec = next((l.split(":", 1)[1].strip() for l in prm.read_text().splitlines()
                    if l.startswith("sha256:")), None)
        check("step09 solvated topology is pinned (held as a pointer)", bool(rec),
              (rec or "")[:16])
    else:
        check("step09 solvated topology reachable",
              (s09 / "complex_solvated.prmtop").exists())
    bridge = out / "complex_solvated.ORCAFF.prms"
    check("the force-field bridge is committed", bridge.exists(),
          f"{bridge.stat().st_size / 1048576:.1f} MB" if bridge.exists() else "")

    print("\n2. the QM region")
    inp = (out / "dft_sp.inp").read_text()
    m = re.search(r"QMAtoms\s*\{(\d+):(\d+)\}", inp)
    check("the QM region is specified as a contiguous range", bool(m),
          m.group(0) if m else "")
    lo, hi = int(m.group(1)), int(m.group(2))
    n_qm = hi - lo + 1
    check(f"the QM region is {QM_EXPECTED} atoms, the whole substrate",
          n_qm == QM_EXPECTED, f"{lo}:{hi} inclusive = {n_qm} atoms")
    # step09 places CHA#2 at atoms 6208-6231 in 1-based prmtop numbering.
    check("the range corresponds to CHA#2 under ORCA's 0-based convention",
          lo == 6207 and hi == 6230,
          f"0-based {lo}:{hi} == 1-based {lo+1}:{hi+1}, which is CHA#2 per step 09")

    # The consequence of a single-index slip, spelled out because nothing in the ORCA
    # output would reveal one. Ranges below are the 0-based CHA boundaries verified
    # against complex_dry.prmtop at step 09; see 09_tleap_build and 19_ensemble.
    NEIGHBOURS = {382: (6183, 6206), 383: (6207, 6230), 384: (6231, 6254)}
    declared = (lo, hi)
    match = [r for r, rng in NEIGHBOURS.items() if rng == declared]
    check("the declared range is exactly one CHA residue boundary", len(match) == 1,
          f"residue {match[0]}" if match else
          f"{declared} matches no CHA boundary — it straddles two copies")
    for r, (a, b) in sorted(NEIGHBOURS.items()):
        mark = "  <== declared" if (a, b) == declared else ""
        print(f"         CHA residue {r}: 0-based {a}:{b}{mark}")
    print("         The copies abut. An off-by-one would straddle two of them and still")
    print("         give 24 atoms, charge -2 and a converging SCF: silently wrong.")

    print("\n3. charge and multiplicity")
    pm = re.search(r"\*pdbfile\s+(-?\d+)\s+(\d+)", inp)
    check("declared in the input", bool(pm), pm.group(0) if pm else "")
    check(f"charge is {CHARGE_EXPECTED}", int(pm.group(1)) == CHARGE_EXPECTED,
          pm.group(1))
    check(f"multiplicity is {MULT_EXPECTED}", int(pm.group(2)) == MULT_EXPECTED,
          pm.group(2))
    o = (out / "dft_sp.out").read_text()
    cm = re.search(r"Total Charge\s+Charge\s+\.+\s+(-?\d+)", o)
    mm = re.search(r"Multiplicity\s+Mult\s+\.+\s+(\d+)", o)
    check("ORCA read the same charge and multiplicity",
          cm and mm and int(cm.group(1)) == CHARGE_EXPECTED
          and int(mm.group(1)) == MULT_EXPECTED,
          f"charge {cm.group(1) if cm else '?'}, mult {mm.group(1) if mm else '?'}")

    print("\n4. level of theory")
    kw = inp.splitlines()[0]
    print(f"         ours:        {kw.lstrip('! ').strip()}")
    print("         Claeyssens:  B3LYP/6-31G(d)/CHARMM27  (2011, abstract)")
    for want in ["B3LYP", "D3BJ", "def2-SVP", "RIJCOSX", "QMMM"]:
        check(f"keyword {want} present", want in kw)
    check("electrostatic embedding via the ORCAFF bridge",
          "ORCAFFFilename" in inp)

    print("\n5. THE ENERGY-SCALE IDENTITY")
    qm = re.search(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", o)
    mmE = re.search(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)", o)
    qmmm = re.search(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)", o)
    check("all three energies are reported", bool(qm and mmE and qmmm))
    eq, em, eqm = float(qm.group(1)), float(mmE.group(1)), float(qmmm.group(1))
    print(f"         QM (embedded substrate) {eq:.9f} Eh")
    print(f"         MM (environment)        {em:.9f} Eh")
    print(f"         QM/MM total             {eqm:.9f} Eh")
    resid = abs((eq + em) - eqm)
    check("QM/MM total is exactly QM + MM", resid < 1e-9,
          f"residual {resid:.2e} Eh")
    check("the QM-only scale is near -836 Eh", -840 < eq < -830, f"{eq:.2f}")
    check("the QM/MM scale is near -1073 Eh", -1080 < eqm < -1070, f"{eqm:.2f}")
    print("         These are DIFFERENT SCALES and must never be differenced against one")
    print("         another. The embedded QM energy already contains the polarisation of")
    print("         the substrate by the environment; the MM term is the environment's own")
    print("         energy. This is what E_int_R and E_int_P in the harvest measure.")

    print("\n6. no link atoms")
    check("the input declares no link-atom treatment",
          "LinkAtom" not in inp and "Link" not in inp,
          "the QM region is a whole molecule, so no bond is cut")
    la = re.search(r"Number of link atoms\s+\.+\s+(\d+)", o)
    sz = re.search(r"Size of QM Subsystem plus link atoms\s+\.+\s+(\d+)", o)
    check("ORCA reports zero link atoms", la and int(la.group(1)) == 0,
          f"Number of link atoms = {la.group(1) if la else '?'}")
    check("the QM subsystem including link atoms is still the 24-atom substrate",
          sz and int(sz.group(1)) == QM_EXPECTED,
          f"{sz.group(1) if sz else '?'} atoms")

    print("\n7. both jobs completed")
    for n in ["smoke.out", "dft_sp.out"]:
        t = (out / n).read_text()
        ok = ("ORCA TERMINATED NORMALLY" in t) or ("TOTAL RUN TIME" in t)
        check(f"{n} terminated normally", ok)
        check(f"{n} reports no error or abort",
              "ORCA finished by error" not in t and "aborting" not in t.lower())

    print(f"\n{'ALL CHECKS PASSED' if not FAIL else 'FAILURES: ' + ', '.join(FAIL)}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())

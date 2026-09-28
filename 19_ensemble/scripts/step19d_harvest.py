#!/usr/bin/env python3
"""
step19d_harvest.py - extract the committable record from converged ensemble
reactant optimisations, before any scratch is cleared.

WHY THIS EXISTS
---------------
Each ensemble reactant optimisation leaves roughly 2.5 GB in its frame
directory, almost all of it ORCA scratch. The part that matters for the thesis
is small: the optimised QM-region geometry, the final energy, the reacting-bond
distances, the region definitions and the convergence trace. This script writes
those to a directory that can be committed, so that clearing scratch afterwards
loses nothing of record.

WHERE THE REGION DEFINITIONS COME FROM
--------------------------------------
ORCA writes `<base>-qm_atoms.tmp` and `<base>-active_atoms.tmp` while running
and deletes them on normal termination, so they are absent from every converged
frame. The durable record is `reactant_opt.inp`, which declares both regions and
is retained. This script parses `QMAtoms` from that file and copies the whole
input alongside the extracted geometry.

The QMAtoms block indexes atoms from zero, so index n selects the (n+1)th ATOM
record of the pdb. For this system `QMAtoms {6207:6230}` selects pdb serials
6208-6231, which is CHA residue 383, matching the selection made at step 12b.
The residues actually selected are written into the xyz comment line so that the
mapping can be checked rather than assumed.

THIS SCRIPT DELETES NOTHING.

Usage
-----
  python3 step19d_harvest.py                 # all converged frames
  python3 step19d_harvest.py 5680 7310       # named frames
"""

import os
import re
import sys
import math
import shutil

HOME = os.path.expanduser("~")
ROOT = f"{HOME}/system_development/05_qmmm/19_ensemble"
OUT = f"{HOME}/system_development/05_qmmm/19_ensemble_harvest"

# zero-based indices within the extracted QM region, from reacting_atoms.tsv
I_O3, I_C4, I_C1, I_C6 = 7, 8, 0, 12


def converged(d):
    f = os.path.join(d, "reactant_opt.out")
    if not os.path.exists(f):
        return False
    with open(f, errors="replace") as fh:
        return "The minimization has converged" in fh.read()


def qm_indices(inp_path):
    """Parse the QMAtoms block. Handles ranges (a:b), lists, and mixtures."""
    txt = open(inp_path, errors="replace").read()
    m = re.search(r"QMAtoms\s*\{([^}]*)\}", txt)
    if not m:
        return None, None
    spec = " ".join(m.group(1).split())
    idx = []
    for tok in spec.replace(",", " ").split():
        if ":" in tok:
            a, b = tok.split(":")
            idx.extend(range(int(a), int(b) + 1))
        else:
            idx.append(int(tok))
    return sorted(set(idx)), spec


def extract(d, o):
    inp = os.path.join(d, "reactant_opt.inp")
    pdb = os.path.join(d, "reactant_opt.pdb")
    if not (os.path.exists(inp) and os.path.exists(pdb)):
        return None, "missing reactant_opt.inp or reactant_opt.pdb"

    idx, spec = qm_indices(inp)
    if idx is None:
        return None, "no QMAtoms block in reactant_opt.inp"

    atoms = [l for l in open(pdb, errors="replace")
             if l.startswith(("ATOM", "HETATM"))]
    if max(idx) >= len(atoms):
        return None, f"QMAtoms index {max(idx)} exceeds {len(atoms)} pdb atoms"

    sel = [atoms[i] for i in idx]
    res = sorted({(l[17:20].strip(), l[22:26].strip()) for l in sel})

    xyz = os.path.join(o, "reactant_qm.xyz")
    coords = []
    with open(xyz, "w") as fh:
        fh.write(f"{len(sel)}\n")
        fh.write(f"optimised QM region from reactant_opt.pdb; "
                 f"QMAtoms {{{spec}}}; residues {res}\n")
        for l in sel:
            el = (l[76:78].strip() or l[12:16].strip()[0]).capitalize()
            x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
            coords.append((x, y, z))
            fh.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")

    def dist(i, j):
        a, b = coords[i], coords[j]
        return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))

    if len(coords) <= max(I_O3, I_C4, I_C1, I_C6):
        return None, f"only {len(coords)} atoms extracted"
    return (dist(I_O3, I_C4), dist(I_C1, I_C6), res), None


def main():
    args = sys.argv[1:]
    frames = []
    if args:
        frames = [f"frame_{int(a):05d}" for a in args]
    else:
        for name in sorted(os.listdir(ROOT)):
            if name.startswith("frame_") and converged(os.path.join(ROOT, name)):
                frames.append(name)

    if not frames:
        sys.exit("no converged frames found")
    os.makedirs(OUT, exist_ok=True)
    print(f"frames: {' '.join(f[6:] for f in frames)}\n")

    rows = []
    for name in frames:
        d = os.path.join(ROOT, name)
        o = os.path.join(OUT, name)
        if not os.path.isdir(d):
            print(f"  skip {name}: no directory"); continue
        if not converged(d):
            print(f"  skip {name}: not converged"); continue
        os.makedirs(o, exist_ok=True)

        shutil.copy(os.path.join(d, "reactant_opt.inp"),
                    os.path.join(o, "reactant_opt.inp"))
        csv = os.path.join(d, "reactant_opt-minimize-ener.csv")
        if os.path.exists(csv):
            shutil.copy(csv, os.path.join(o, "minimize-ener.csv"))

        got, err = extract(d, o)
        if err:
            print(f"  {name[6:]}: EXTRACTION FAILED - {err}")
            brk = frm = None
            res = []
        else:
            brk, frm, res = got

        out = open(os.path.join(d, "reactant_opt.out"), errors="replace").read()
        # ORCA writes three energies for a QM/MM run:
        #   FINAL SINGLE POINT ENERGY (MM)      the classical subsystem
        #   FINAL SINGLE POINT ENERGY           the QM subsystem alone
        #   FINAL SINGLE POINT ENERGY (QM/MM)   the total, and the only one from
        #                                       which a barrier may be formed
        # All three are recorded. The QM-only value is kept because comparing it
        # with the bare-substrate energy measures how far the protein field
        # polarises the substrate.
        def last(pattern):
            hits = re.findall(pattern, out)
            return hits[-1] if hits else "NA"

        E_total = last(r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)")
        E_qm = last(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
        E_mm = last(r"FINAL SINGLE POINT ENERGY \(MM\)\s+(-?\d+\.\d+)")

        # consistency: the total should be the sum of the two parts
        note = ""
        try:
            if abs((float(E_qm) + float(E_mm)) - float(E_total)) > 1e-6:
                note = "  SUM MISMATCH"
        except ValueError:
            note = "  incomplete energies"
        E = E_total
        steps = sum(1 for _ in open(csv)) if os.path.exists(csv) else "NA"
        wt = re.search(r"TOTAL RUN TIME:\s*(.+)", out)
        wt = wt.group(1).strip() if wt else "NA"

        rows.append((name[6:], E_total, E_qm, E_mm, steps, wt, brk, frm))
        if brk is not None:
            print(f"  {name[6:]}  E(QM/MM)={E_total}  E(QM)={E_qm}  steps={steps}  "
                  f"break={brk:.3f}  form={frm:.3f}  res={res}{note}")

    summ = os.path.join(OUT, "reactant_summary.tsv")
    with open(summ, "w") as fh:
        fh.write("frame\tE_QMMM_Eh\tE_QM_Eh\tE_MM_Eh\tsteps\twalltime"
                 "\tbreak_O3_C4\tform_C1_C6\tr\n")
        for f, Et, Eq, Em, s, w, b, m in rows:
            bs = f"{b:.3f}" if b is not None else "NA"
            ms = f"{m:.3f}" if m is not None else "NA"
            rs = f"{b-m:.3f}" if (b is not None and m is not None) else "NA"
            fh.write(f"{f}\t{Et}\t{Eq}\t{Em}\t{s}\t{w}\t{bs}\t{ms}\t{rs}\n")

    print(f"\nwrote {summ}")
    ok = sum(1 for r in rows if r[6] is not None)
    print(f"{ok} of {len(rows)} frames extracted successfully")
    print("\nNOTE: E_QMMM is the total and is the only energy from which a barrier "
          "may be\nformed. Absolute values are not comparable between frames, since "
          "each carries a\ndifferent solvent configuration in the MM term; only "
          "differences within a frame\nare meaningful.")
    print("\nReview the summary, then commit the harvest before clearing any scratch.")


if __name__ == "__main__":
    main()

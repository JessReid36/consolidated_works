#!/usr/bin/env python3
"""
step19e_harvest_barriers.py - assemble the committable record of the barrier
ensemble.

WHAT THIS COLLECTS, AND WHY EACH IS NEEDED
    Nothing may be quoted in the write-up that is not traceable to a committed
    file. Four things are therefore extracted for every frame that produced a
    converged climbing image:

      1. The barrier, together with the path summary table it was read from, so
         that the number can be checked against its source without rerunning.
      2. The in vacuo single-point energies at the reactant, transition state and
         product geometries, and the extracted substrate geometries themselves,
         so that the separation of intrinsic chemistry from environment can be
         reproduced.
      3. The transition-state geometry from the converged climbing image.
      4. A single table joining all of the above, which is what the write-up
         quotes from.

    The QM/MM energies are taken from the optimisation outputs rather than
    recomputed. The barrier is read from the row of the path summary carrying the
    climbing-image marker, which in ORCA 6.0.1 is "<= CI".

WHAT IS NOT COLLECTED
    Wavefunction and density binaries, and scratch files. They are large and are
    regenerable from the inputs, which are collected.

Usage:
    python3 step19e_harvest_barriers.py                 # default paths
    python3 step19e_harvest_barriers.py --out <dir>
"""
import argparse
import math
import os
import re
import shutil
import sys

H = 627.5094740631

# The ten frames satisfying the full three-condition near-attack criterion.
# Frames failing theta2, and any frame without a converged climbing image, are
# harvested too but flagged, so that the decision to report or exclude them is
# made in the write-up and not silently here.
FULL_NAC = ["02450", "04085", "09025", "09900", "10775", "11630", "12485", "14155", "17505", "19185", "20000", "21634", "23268", "24883", "26495", "33320", "34991", "36665", "38344", "40485", "41786", "42436", "43087", "43738", "44388", "45688", "46990", "47641", "48292", "49592", "50892", "52192", "53493", "54144", "54795", "55446", "56098", "57397", "58698", "59999"]
NEAR_NAC = ["05680", "07310", "15825"]
OTHER = ["08170"]


def path_summary(txt):
    """The PATH SUMMARY block, as written, and the climbing-image barrier."""
    blk = re.search(r"-+\s*\n\s*PATH SUMMARY.*?(?=Straight line distance)",
                    txt, re.S)
    if not blk:
        return None, None
    body = blk.group(0)
    ci = re.search(r"^\s*(\d+)\s+[-\d.]+\s+(-?\d+\.\d+)\s+([-\d.]+)\s+\S+\s+\S+"
                   r"\s*<=\s*CI\s*$", body, re.M)
    if ci:
        return body, (int(ci.group(1)), float(ci.group(2)), float(ci.group(3)))
    return body, None


def last_energy(path, qmmm=False):
    if not os.path.exists(path):
        return None
    pat = (r"FINAL SINGLE POINT ENERGY \(QM/MM\)\s+(-?\d+\.\d+)" if qmmm
           else r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)\s*$")
    hits = re.findall(pat, open(path, errors="replace").read(), re.M)
    return float(hits[-1]) if hits else None


def vac_energy(path):
    if not os.path.exists(path):
        return None
    t = open(path, errors="replace").read()
    if "TERMINATED NORMALLY" not in t:
        return None
    hits = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", t)
    return float(hits[-1]) if hits else None


def stats(v):
    v = [x for x in v if x is not None and x == x]
    if len(v) < 2:
        return (float("nan"),) * 2
    m = sum(v) / len(v)
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def main():
    home = os.path.expanduser("~")
    ap = argparse.ArgumentParser()
    ap.add_argument("--ensdir", default=f"{home}/system_development/05_qmmm/19_ensemble")
    ap.add_argument("--vacdir", default=f"{home}/system_development/05_qmmm/20_invacuo")
    ap.add_argument("--out", default=f"{home}/system_development/05_qmmm/19_ensemble_barriers")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    rows = []

    for group, frames in (("full_NAC", FULL_NAC),
                          ("near_NAC", NEAR_NAC),
                          ("other", OTHER)):
        for f in frames:
            d = f"{args.ensdir}/frame_{f}"
            fd = f"{args.out}/frame_{f}"
            os.makedirs(fd, exist_ok=True)
            rec = {"frame": f, "group": group}

            # ---- 1. the barrier and the table it came from
            nebout = f"{d}/neb.out"
            if os.path.exists(nebout):
                txt = open(nebout, errors="replace").read()
                body, ci = path_summary(txt)
                if body:
                    with open(f"{fd}/path_summary.txt", "w") as fh:
                        fh.write(f"# from {nebout}\n")
                        fh.write("# the barrier is the dE(kcal/mol) of the row "
                                 "marked '<= CI'\n\n")
                        fh.write(body)
                if ci:
                    rec["ci_image"], rec["E_ci_Eh"], rec["barrier"] = ci
                else:
                    rec["barrier"] = None
                conv = "THE NEB OPTIMIZATION HAS CONVERGED" in txt
                rec["neb_converged"] = "yes" if conv else "no"
                # the band's first segment, a diagnostic of the guess
                m = re.search(r"D\(\s*0-\s*1\)\s*=\s*([\d.]+)", txt)
                rec["D01"] = float(m.group(1)) if m else None
                for src in ("neb.inp",):
                    if os.path.exists(f"{d}/{src}"):
                        shutil.copy(f"{d}/{src}", f"{fd}/{src}")

            # ---- 2 and 3. geometries and in vacuo energies
            for tag, name in (("R", "reactant"), ("TS", "transition_state"),
                              ("P", "product")):
                g = f"{args.vacdir}/{f}_{tag}.xyz"
                if os.path.exists(g):
                    shutil.copy(g, f"{fd}/{name}_qm.xyz")
                o = f"{args.vacdir}/sp_{f}_{tag}.out"
                rec[f"E_vac_{tag}"] = vac_energy(o)
                i = f"{args.vacdir}/sp_{f}_{tag}.inp"
                if os.path.exists(i):
                    shutil.copy(i, f"{fd}/invacuo_{name}.inp")

            # ---- the QM/MM energies, from the optimisations
            rec["E_qmmm_R"] = last_energy(f"{d}/reactant_opt.out", qmmm=True)
            rec["E_qmmm_P"] = last_energy(f"{d}/product_opt.out", qmmm=True)
            rec["E_qm_R"] = last_energy(f"{d}/reactant_opt.out")
            rec["E_qm_P"] = last_energy(f"{d}/product_opt.out")

            # ---- derived
            vR, vT, vP = rec.get("E_vac_R"), rec.get("E_vac_TS"), rec.get("E_vac_P")
            rec["barrier_vac"] = (vT - vR) * H if (vR and vT) else None
            # ---- stab_TS: the differential TS stabilisation of Claeyssens et
            # al. (Chem. Commun. 2005, 5068; Org. Biomol. Chem. 2011, 9, 1578).
            # Their E_INTERACTION(TS) - E_INTERACTION(reactant) reduces exactly
            # to the QM/MM barrier minus the in vacuo barrier at the same
            # geometries: the MM-only terms cancel in that difference. That is
            # why it can be formed here without an embedded QM energy at the
            # climbing image, which the path summary does not report.
            # Negative means the environment lowers the barrier.
            b, bv = rec.get("barrier"), rec.get("barrier_vac")
            rec["stab_TS"] = (float(b) - float(bv)) if (b is not None and bv is not None) else None

            # ---- E_int_*: the EMBEDDED QM energy of a structure minus its
            # gas-phase energy at the same geometry. Interaction plus
            # polarisation of the QM region by the MM charges. Large and
            # negative for this dianion in a cationic site.
            # NOT a stabilisation. NEVER difference these against any quantity
            # built from a QM/MM TOTAL energy (e.g. E_ci_Eh): the total carries
            # the MM energy of ~55680 atoms and the result is meaningless.
            # The former stab_TS_total column did exactly that and was wrong by
            # about -240 Eh; it has been removed.
            qR, qP = rec.get("E_qm_R"), rec.get("E_qm_P")
            rec["E_int_R"] = (qR - vR) * H if (qR and vR) else None
            rec["E_int_P"] = (qP - vP) * H if (qP and vP) else None
            rows.append(rec)

    # ------------------------------------------------------------ the table
    cols = ["frame", "group", "neb_converged", "ci_image", "D01", "barrier",
            "barrier_vac", "stab_TS", "E_ci_Eh", "E_qmmm_R", "E_qmmm_P",
            "E_qm_R", "E_qm_P", "E_vac_R", "E_vac_TS", "E_vac_P",
            "E_int_R", "E_int_P"]
    out = f"{args.out}/ensemble_barriers.tsv"
    with open(out, "w") as fh:
        fh.write("# barrier ensemble, harvested by step19e_harvest_barriers.py\n")
        fh.write("# barrier: kcal/mol, from the '<= CI' row of the path summary\n")
        fh.write("# barrier_vac: kcal/mol, in vacuo at the same geometries\n")
        fh.write("# stab_TS: kcal/mol, barrier - barrier_vac. The differential TS\n")
        fh.write("#   stabilisation of Claeyssens et al. Negative = the environment\n")
        fh.write("#   lowers the barrier. This is the column to quote.\n")
        fh.write("# E_int_R, E_int_P: kcal/mol, EMBEDDED QM energy minus gas-phase\n")
        fh.write("#   energy at the same geometry. Interaction plus polarisation,\n")
        fh.write("#   NOT a stabilisation. Never difference against E_ci_Eh or any\n")
        fh.write("#   other QM/MM total.\n")
        fh.write("# energies: Eh.\n")
        fh.write("# group: full_NAC frames are those reported in the write-up\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                ("" if r.get(c) is None else
                 (f"{r[c]:.6f}" if isinstance(r.get(c), float) else str(r[c])))
                for c in cols) + "\n")

    # ------------------------------------------------------------- summary
    print(f"wrote {out}")
    print(f"{'frame':<8}{'group':<10}{'NEB':>5}{'barrier':>9}{'vac':>9}{'D(0-1)':>9}")
    print("-" * 50)
    for r in rows:
        print(f"{r['frame']:<8}{r['group']:<10}{r.get('neb_converged','-'):>5}"
              f"{(r['barrier'] if r.get('barrier') is not None else float('nan')):>9.2f}"
              f"{(r['barrier_vac'] if r.get('barrier_vac') is not None else float('nan')):>9.2f}"
              f"{(r['D01'] if r.get('D01') is not None else float('nan')):>9.3f}")

    full = [r for r in rows if r["group"] == "full_NAC" and r.get("barrier")]
    if len(full) >= 6:
        mb, sb = stats([r["barrier"] for r in full])
        mv, sv = stats([r["barrier_vac"] for r in full])
        print(f"\nfull near-attack set, n = {len(full)}")
        print(f"  QM/MM barrier     {mb:7.2f} +/- {sb:5.2f} kcal/mol")
        print(f"  in vacuo barrier  {mv:7.2f} +/- {sv:5.2f} kcal/mol")
        print(f"  sigma in kJ/mol   {sb*4.184:7.2f}")
        print("\n  Ryde (2017) recommends the standard deviation as the")
        print("  convergence criterion and surveys 24 QM/MM studies spanning")
        print("  0.6 to 97 kJ/mol, 73% of which exceed 10 kJ/mol.")
    print("\nCommit the output directory, then the write-up may quote from it.")


if __name__ == "__main__":
    main()

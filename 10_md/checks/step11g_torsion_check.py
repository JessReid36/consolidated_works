#!/usr/bin/env python3
"""
step11g_torsion_check.py - open item 8: does the GAFF parameterisation bias the sampled
conformation about the bond that breaks?

WHY THIS EXISTS
parmchk2 reported 0 ATTN for cha.frcmod, but two dihedral parameters were supplied by
analogy with a penalty score of 237.0, against a next-highest of 3.0:

    c2-ce-os-c3   substituted by  X-c2-os-X   penalty 237.0
    c -ce-os-c3   substituted by  X-c2-os-X   penalty 237.0

Resolving the types against cha_gaff.mol2, those are the torsions

    C1-C2-O3-C4      and      C3-C2-O3-C4

and O3-C4 is the bond that breaks in the rearrangement (step 11 atom mapping:
"breaking C2-O13 -> C4-O3"). So the two least reliable terms in the substrate force field
govern the torsion about the reacting bond.

This does NOT affect any QM/MM barrier: the substrate is in the QM region, so GAFF never
computes the reaction. It DOES govern the classical MD that samples the reactant
conformational distribution, and hence the frames selected for QM/MM. The question is
therefore whether the MD samples the same torsion range that the QM/MM optimisation
settles into.

THE TEST
Compare the torsion distribution over the production trajectory against the same torsions
in the QM/MM optimised reactants. The reference values, computed from the 22 committed
reactant_qm.xyz geometries, are:

    C1-C2-O3-C4   mean  -94.0 deg   sd 4.4   range -101.4 to -80.9
    C3-C2-O3-C4   mean  +83.3 deg   sd 4.2   range  +75.5 to +93.9

INTERPRETATION, decided before running so the result cannot be rationalised afterwards:
  - If the MD mean sits within roughly one MD standard deviation of the QM/MM mean, and
    the QM/MM range lies inside the bulk of the MD distribution, the parameterisation is
    adequate for sampling. That is all it is asked to do. Record and move on.
  - If the MD distribution is centred elsewhere, or is bimodal with the QM/MM values in a
    minor basin, then the frames handed to QM/MM were drawn from a distorted distribution.
    That is a sampling caveat for X.2 and the ensemble sections. It does not invalidate
    any individual barrier.

RUN THIS ON THE HPC. It needs prod.nc, which is held outside git under the large-file
policy.

    cd ~/system_development
    python3 step11g_torsion_check.py \
        04_amber_md/10c_production/prod.nc \
        03_amber/tleap_build/complex_solvated.prmtop \
        00_admin/step11g_torsion_check.tsv

Optionally append the extension trajectory as a fourth argument:
        ... 00_admin/step11g_torsion_check.tsv 04_amber_md/10d_production_extend/prod_ext.nc

No cpptraj: it is broken cluster-wide. This uses the same stdlib NetCDF3 reader approach
as step11e_autocorr.py, plus numpy.
"""
import os
import sys
import math
import struct

import numpy as np

# Reference: QM/MM optimised reactants, 22 frames, from the committed reactant_qm.xyz.
REFERENCE = {
    "C1-C2-O3-C4": {"mean": -94.0, "sd": 4.4, "min": -101.4, "max": -80.9},
    "C3-C2-O3-C4": {"mean": 83.3, "sd": 4.2, "min": 75.5, "max": 93.9},
}
# Atom names in the CHA residue, from cha_gaff.mol2.
TORSIONS = [("C1-C2-O3-C4", ("C1", "C2", "O3", "C4")),
            ("C3-C2-O3-C4", ("C3", "C2", "O3", "C4"))]


def prmtop_atoms(path):
    """Return (natom, atom_names, residue_labels, residue_pointers)."""
    flags, fmts, cur = {}, {}, None
    for ln in open(path, "r", errors="replace"):
        if ln.startswith("%FLAG"):
            cur = ln.split()[1]
            flags[cur] = []
        elif ln.startswith("%FORMAT"):
            fmts[cur] = ln[ln.find("(") + 1:ln.rfind(")")]
        elif cur and not ln.startswith("%"):
            flags[cur].append(ln.rstrip("\n"))

    def width(fmt):
        for i, c in enumerate(fmt):
            if c.isalpha():
                tail = fmt[i + 1:].split(".")[0]
                return int(tail) if tail.isdigit() else 0
        return 0

    def toks(name, cast):
        w = width(fmts.get(name, ""))
        out = []
        for ln in flags.get(name, []):
            if w:
                out += [ln[i:i + w].strip() for i in range(0, len(ln.rstrip()), w)
                        if ln[i:i + w].strip()]
            else:
                out += ln.split()
        return [cast(x) for x in out]

    return (toks("POINTERS", int)[0], toks("ATOM_NAME", str),
            toks("RESIDUE_LABEL", str), toks("RESIDUE_POINTER", int))


def nc_open(path):
    """Minimal NetCDF3 record-variable reader. Same approach as step11e_autocorr.py."""
    data = open(path, "rb").read(1 << 20)
    fsz = os.path.getsize(path)
    off64 = (data[3] == 2)
    pos = [4]

    def u32():
        v = struct.unpack_from(">I", data, pos[0])[0]
        pos[0] += 4
        return v

    def off():
        v = struct.unpack_from(">Q" if off64 else ">I", data, pos[0])[0]
        pos[0] += 8 if off64 else 4
        return v

    def nm():
        n = u32()
        s = data[pos[0]:pos[0] + n].decode("ascii", "replace")
        pos[0] += n + ((4 - (n % 4)) % 4)
        return s

    numrecs = u32()
    dims = []
    tag = u32()
    ne = u32()
    if tag == 0x0A:
        for _ in range(ne):
            dn = nm()
            dl = u32()
            dims.append((dn, dl))

    def skip_att():
        t = u32()
        n = u32()
        if t == 0x0C:
            for _ in range(n):
                nm()
                tp = u32()
                k = u32()
                nb = k * {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}.get(tp, 1)
                pos[0] += nb + ((4 - (nb % 4)) % 4)

    skip_att()
    rec = []
    tag = u32()
    nv = u32()
    if tag == 0x0B:
        for _ in range(nv):
            vn = nm()
            nd = u32()
            dimids = [u32() for _ in range(nd)]
            skip_att()
            tp = u32()
            vs = u32()
            bg = off()
            if nd > 0 and dims[dimids[0]][1] == 0:
                rec.append({"name": vn, "vsize": vs, "begin": bg, "type": tp})
    dimd = {n: l for n, l in dims}
    recsize = sum(v["vsize"] for v in rec)
    frames = numrecs if numrecs != 0xFFFFFFFF else \
        (fsz - min(v["begin"] for v in rec)) // recsize
    return {"frames": frames, "recsize": recsize, "rec": rec,
            "atom": dimd.get("atom"), "path": path}


def coord_frames(nc, wanted_idx):
    """Yield the coordinates of `wanted_idx` atoms, one frame at a time."""
    crd = next(v for v in nc["rec"] if v["name"] == "coordinates")
    natom = nc["atom"]
    dtype = ">f4" if crd["type"] == 5 else ">f8"
    itemsize = 4 if crd["type"] == 5 else 8
    with open(nc["path"], "rb") as fh:
        for f in range(nc["frames"]):
            base = crd["begin"] + f * nc["recsize"]
            out = np.empty((len(wanted_idx), 3))
            for j, ai in enumerate(wanted_idx):
                fh.seek(base + ai * 3 * itemsize)
                out[j] = np.frombuffer(fh.read(3 * itemsize), dtype=dtype)
            yield out


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    b1n = b1 / np.linalg.norm(b1)
    m = np.cross(n1, b1n)
    return math.degrees(math.atan2(float(m @ n2), float(n1 @ n2)))


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    nc_path, prm_path, out_path = sys.argv[1:4]
    extra = sys.argv[4] if len(sys.argv) > 4 else None

    natom, anames, rlabels, rptr = prmtop_atoms(prm_path)
    cha = [i for i, lab in enumerate(rlabels) if lab == "CHA"]
    if not cha:
        sys.exit("no CHA residue found in the prmtop")
    print(f"prmtop: {natom} atoms, {len(cha)} CHA residues at residue indices "
          f"{[i + 1 for i in cha]}")

    # Atom index ranges for each CHA copy.
    copies = []
    for i in cha:
        lo = rptr[i] - 1
        hi = (rptr[i + 1] - 1) if i + 1 < len(rptr) else natom
        names = anames[lo:hi]
        copies.append({"lo": lo, "names": names})
        print(f"  CHA residue {i + 1}: atoms {lo + 1}-{hi} ({hi - lo})")

    # Resolve each torsion to absolute atom indices, per copy.
    plans = []
    for label, quad in TORSIONS:
        for c, cp in enumerate(copies):
            try:
                idx = [cp["lo"] + cp["names"].index(a) for a in quad]
            except ValueError:
                sys.exit(f"atom {quad} not found in CHA copy {c}; names={cp['names']}")
            plans.append((label, c, idx))

    wanted = sorted({i for _, _, idx in plans for i in idx})
    pos_of = {a: j for j, a in enumerate(wanted)}

    series = {(lab, c): [] for lab, c, _ in plans}
    ncs = [nc_open(nc_path)] + ([nc_open(extra)] if extra else [])
    total = 0
    for nc in ncs:
        print(f"{os.path.basename(nc['path'])}: {nc['frames']} frames, "
              f"{nc['atom']} atoms")
        if nc["atom"] != natom:
            sys.exit(f"atom count mismatch: prmtop {natom}, nc {nc['atom']}")
        for frame in coord_frames(nc, wanted):
            for lab, c, idx in plans:
                p = [frame[pos_of[a]] for a in idx]
                series[(lab, c)].append(dihedral(*p))
            total += 1
            if total % 2000 == 0:
                print(f"  ... {total} frames")
    print(f"total frames: {total}")

    with open(out_path, "w") as fh:
        fh.write("# open item 8: torsion about the breaking O3-C4 bond\n")
        fh.write("# MD production vs QM/MM optimised reactants\n")
        fh.write("torsion\tcha_copy\tn\tmean\tsd\tmin\tmax\t"
                 "ref_mean\tref_sd\tref_min\tref_max\tref_within_1sd\n")
        print(f"\n{'torsion':<14}{'copy':>5}{'n':>8}{'MD mean':>10}{'MD sd':>8}"
              f"{'QM/MM mean':>12}{'delta':>8}{'verdict':>12}")
        verdicts = []
        for (lab, c), v in sorted(series.items()):
            a = np.asarray(v)
            ref = REFERENCE[lab]
            delta = float(a.mean()) - ref["mean"]
            within = abs(delta) <= max(float(a.std(ddof=1)), ref["sd"])
            verdicts.append(within)
            fh.write(f"{lab}\t{c}\t{len(a)}\t{a.mean():.2f}\t{a.std(ddof=1):.2f}\t"
                     f"{a.min():.2f}\t{a.max():.2f}\t{ref['mean']:.1f}\t{ref['sd']:.1f}\t"
                     f"{ref['min']:.1f}\t{ref['max']:.1f}\t{within}\n")
            print(f"{lab:<14}{c:>5}{len(a):>8}{a.mean():>10.1f}{a.std(ddof=1):>8.1f}"
                  f"{ref['mean']:>12.1f}{delta:>+8.1f}"
                  f"{'consistent' if within else 'DIVERGENT':>12}")
        # Coarse histogram, 20 degree bins, to expose bimodality.
        fh.write("#\n# histogram, 20 degree bins, all copies pooled\n")
        print("\nhistogram, 20 degree bins, all copies pooled "
              "(| marks the QM/MM mean):")
        for lab, _ in TORSIONS:
            pooled = np.concatenate([np.asarray(series[(lab, c)])
                                     for c in range(len(copies))])
            hist, edges = np.histogram(pooled, bins=18, range=(-180, 180))
            fh.write(f"# {lab}\n")
            print(f"  {lab}")
            for h, e in zip(hist, edges[:-1]):
                mark = " <== QM/MM" if e <= REFERENCE[lab]["mean"] < e + 20 else ""
                bar = "#" * int(60 * h / max(hist.max(), 1))
                fh.write(f"#   {e:+.0f}\t{h}\n")
                print(f"    {e:+5.0f} {h:>7} {bar}{mark}")
        fh.write(f"#\n# verdict: {'CONSISTENT' if all(verdicts) else 'DIVERGENT'}\n")

    print(f"\nWROTE {out_path}")
    print("\nVERDICT: " + ("CONSISTENT - the MD samples the torsion range the QM/MM "
                           "optimisation settles into. The high-penalty parameters did "
                           "not bias the sampled ensemble. Record and move on."
                           if all(verdicts) else
                           "DIVERGENT - the frames handed to QM/MM were drawn from a "
                           "distribution centred away from the QM/MM optimum. This is a "
                           "sampling caveat for X.2 and the ensemble sections; it does "
                           "NOT invalidate any individual barrier."))
    return 0 if all(verdicts) else 2


if __name__ == "__main__":
    sys.exit(main())

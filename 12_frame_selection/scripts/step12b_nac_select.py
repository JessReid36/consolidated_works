#!/usr/bin/env python3
"""
step12b_nac_select.py - select additional reactant frames using the FULL
near-attack-conformer criterion, to extend the barrier ensemble.

WHY
---
The step-12a selection applied only the distance clause of the NAC definition
(d(C1-C6) < 3.7 A). Checking the twelve selected frames against the full
three-condition criterion of Hur & Bruice (PNAS 2003, 100, 12015) - distance
plus two angular conditions on the pi-orbital alignment - showed that nine
satisfy it and three do not:

    frame  5680:  theta2 deviation 21.4 deg
    frame  7310:  theta2 deviation 24.7 deg
    frame 15825:  theta2 deviation 23.0 deg

All three fail on theta2 alone and marginally, against a 20 degree tolerance.
For comparison, the unreactive diequatorial conformers Hur & Bruice identify in
water show theta2 deviations of 110 +/- 9 degrees. These are near-NACs, not
misassigned structures.

Rather than discard them, this script selects THREE ADDITIONAL frames that pass
all three conditions, giving twelve full NACs plus three near-NACs. Nothing is
lost, and the three near-NACs become a test: if near-attack geometry is what
governs the barrier, they should sit systematically higher.

SELECTION
---------
Candidates must satisfy, at the same site used at step 12a:
  * Arg90-O13 contact < 3.2 A          (as step 12a)
  * d(C1-C6) <= 3.7 A                  (NAC condition i)
  * |theta1 - 8.2| <= 20 deg           (NAC condition ii)
  * |theta2 - 18.6| <= 20 deg          (NAC condition iii)

From that pool, frames are chosen to maximise the minimum time separation from
every already-selected frame and from each other, so the additions extend the
sampling rather than duplicating it. The integrated autocorrelation time of the
forming distance is ~15 ps, so any separation of a few hundred ps is ample; the
maximin rule simply spreads the additions as evenly as the pool allows.

Usage
-----
  python3 step12b_nac_select.py                 # add 3 frames
  python3 step12b_nac_select.py --n 3 --stride 5
"""

import argparse
import os
import re
import struct
import sys
import numpy as np

COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66}
TS_T1, TS_T2, TOL, D_NAC = 8.2, 18.6, 20.0, 3.7
ARG_N = {"NE", "NH1", "NH2"}


def parse_prmtop(path):
    flags, fmts, cur = {}, {}, None
    for ln in open(path, errors="replace"):
        ln = ln.rstrip("\n")
        if ln.startswith("%FLAG"):
            cur = ln.split()[1]; flags[cur] = []
        elif ln.startswith("%FORMAT"):
            fmts[cur] = ln[ln.find("(") + 1:ln.rfind(")")]
        elif ln.startswith("%"):
            continue
        elif cur is not None:
            flags[cur].append(ln)

    def width(f):
        m = re.match(r"\s*\d*[aAiIeEfFgG](\d+)", f)
        return int(m.group(1)) if m else None

    def toks(name, cast):
        w = width(fmts.get(name, "")); out = []
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
    data = open(path, "rb").read(1 << 20)
    fsz = os.path.getsize(path)
    off64 = (data[3] == 2)
    pos = [4]

    def u32():
        v = struct.unpack_from(">I", data, pos[0])[0]; pos[0] += 4; return v

    def off():
        v = struct.unpack_from(">Q" if off64 else ">I", data, pos[0])[0]
        pos[0] += 8 if off64 else 4; return v

    def nm():
        n = u32(); s = data[pos[0]:pos[0] + n].decode("ascii", "replace")
        pos[0] += n + ((4 - (n % 4)) % 4); return s

    numrecs = u32(); dims = []; tag = u32(); ne = u32()
    if tag == 0x0A:
        for _ in range(ne):
            dn = nm(); dl = u32(); dims.append((dn, dl))

    def skip_att():
        t = u32(); n = u32()
        if t == 0x0C:
            for _ in range(n):
                nm(); tp = u32(); k = u32()
                nb = k * {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}.get(tp, 1)
                pos[0] += nb + ((4 - (nb % 4)) % 4)

    skip_att(); rec = []; tag = u32(); nv = u32()
    if tag == 0x0B:
        for _ in range(nv):
            vn = nm(); nd = u32(); dimids = [u32() for _ in range(nd)]
            skip_att(); tp = u32(); vs = u32(); bg = off()
            if nd > 0 and dims[dimids[0]][1] == 0:
                rec.append({"name": vn, "vsize": vs, "begin": bg})
    dimd = {n: l for n, l in dims}
    recsize = sum(v["vsize"] for v in rec)
    frames = numrecs if numrecs != 0xFFFFFFFF else \
        (fsz - min(v["begin"] for v in rec)) // recsize
    return {"frames": frames, "recsize": recsize, "rec": rec, "atom": dimd.get("atom")}


def nac_angles(els, xyz, i_c1, i_c6):
    n = len(els)
    D = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(-1))
    bond = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < (COV[els[i]] + COV[els[j]]) * 1.3:
                bond[i, j] = bond[j, i] = True
    ring = []

    def dfs(p):
        if ring:
            return
        if len(p) == 6:
            if bond[p[-1], p[0]]:
                ring.extend(p)
            return
        for k in np.where(bond[p[-1]])[0]:
            if els[k] != "C" or k in p:
                continue
            dfs(p + [k])

    for s in range(n):
        if els[s] == "C":
            dfs([s])
        if ring:
            break
    if not ring:
        return None

    c_ring = i_c6 if i_c6 in ring else i_c1
    c_side = i_c1 if c_ring == i_c6 else i_c6

    def axis(c):
        nb = sorted(np.where(bond[c])[0], key=lambda k: D[c, k])
        if len(nb) < 2:
            return None
        v = np.cross(xyz[nb[0]] - xyz[c], xyz[nb[1]] - xyz[c])
        nv = np.linalg.norm(v)
        return v / nv if nv > 0 else None

    a_r, a_s = axis(c_ring), axis(c_side)
    if a_r is None or a_s is None:
        return None
    v = xyz[c_side] - xyz[c_ring]
    d = float(np.linalg.norm(v))
    u = v / d
    t1 = float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a_r, u)))))))
    t2 = float(np.degrees(np.arccos(min(1.0, abs(float(np.dot(a_s, -u)))))))
    return d, t1, t2


def main():
    ap = argparse.ArgumentParser()
    home = os.path.expanduser("~")
    ap.add_argument("--prmtop",
                    default=f"{home}/system_development/03_amber/tleap_build/complex_solvated.prmtop")
    ap.add_argument("--traj",
                    default=f"{home}/system_development/04_amber_md/10c_production/prod.nc")
    ap.add_argument("--manifest",
                    default=f"{home}/system_development/05_qmmm/12_frame_selection/selection_manifest.tsv")
    ap.add_argument("--out",
                    default=f"{home}/system_development/05_qmmm/12_frame_selection/selection_manifest_nac_extra.tsv")
    ap.add_argument("--site", default="CHA#2")
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--min-ps", type=float, default=11000.0,
                    dest="min_ps",
                    help="exclude frames earlier than this time in the "
                         "production trajectory. The backbone RMSD rises "
                         "monotonically over the first ten nanoseconds and is "
                         "stable thereafter, so earlier frames do not sample "
                         "the equilibrated ensemble. The default carries a "
                         "one-nanosecond margin beyond that point.")
    ap.add_argument("--contact-cut", type=float, default=3.2)
    args = ap.parse_args()

    natom, names, labels, resptr = parse_prmtop(args.prmtop)
    nres = len(labels)
    starts = [p - 1 for p in resptr] + [natom]
    res_of = np.empty(natom, int)
    for ri in range(nres):
        res_of[starts[ri]:starts[ri + 1]] = ri
    argN = np.array([a for a in range(natom)
                     if labels[res_of[a]] == "ARG" and names[a] in ARG_N])

    cha = []
    for ri in range(nres):
        if labels[ri] == "CHA":
            cha.append((ri, {names[a]: a for a in range(starts[ri], starts[ri + 1])}))
    site = int(args.site.split("#")[1]) - 1
    ri, amap = cha[site]

    order = sorted(amap.items(), key=lambda kv: kv[1])
    sub_idx = [a for _, a in order]
    sub_els = [("H" if nm[0] == "H" else nm[0]) for nm, _ in order]
    loc = {a: k for k, a in enumerate(sub_idx)}
    i_c1, i_c6 = loc[amap["C1"]], loc[amap["C6"]]
    a_o3, a_c4 = amap["O3"], amap["C4"]
    cha_resid = ri + 1

    H = nc_open(args.traj)
    cvar = next(v for v in H["rec"] if v["name"] == "coordinates")
    lvar = next((v for v in H["rec"] if v["name"] == "cell_lengths"), None)
    top = max(max(sub_idx), int(argN.max())) + 1

    def frame(i):
        o = cvar["begin"] + i * H["recsize"]
        full = np.array(np.memmap(args.traj, dtype=">f4", mode="r", offset=o,
                                  shape=(natom, 3))[:top], float)
        box = np.array([1e9, 1e9, 1e9])
        if lvar is not None:
            bo = lvar["begin"] + i * H["recsize"]
            box = np.array(np.memmap(args.traj, dtype=">f8", mode="r", offset=bo,
                                     shape=(3,)), float)
        return full, box

    existing = sorted(int(l.split("\t")[1])
                      for l in open(args.manifest).read().splitlines()
                      if not l.startswith("#") and l.strip())
    print(f"already selected: {existing}\n")
    print(f"scanning every {args.stride}th frame from {int(args.min_ps)} ps "
          f"onward at {args.site} against the full NAC criterion ...")

    pool = []
    # One frame is written per picosecond, so the frame index and the time in
    # picoseconds are the same number. Starting the scan at the first
    # equilibrated frame is equivalent to filtering afterwards but avoids
    # reading coordinates that would only be discarded.
    i0 = int(args.min_ps)
    if i0 > 0:
        print(f"excluding the first {i0} frames as pre-equilibration; "
              f"{H['frames'] - i0} remain")
    for i in range(i0, H["frames"], args.stride):
        if i in existing:
            continue
        full, box = frame(i)
        dvec = full[argN] - full[amap["O3"]]
        dvec -= np.round(dvec / box) * box
        contact = float(np.sqrt((dvec * dvec).sum(1)).min())
        if contact >= args.contact_cut:
            continue
        r = nac_angles(sub_els, full[sub_idx], i_c1, i_c6)
        if not r:
            continue
        d, t1, t2 = r
        if d <= D_NAC and abs(t1 - TS_T1) <= TOL and abs(t2 - TS_T2) <= TOL:
            brk = float(np.linalg.norm(full[a_c4] - full[a_o3]))
            pool.append((i, d, t1, t2, contact, brk))

    print(f"  {len(pool)} candidate frames satisfy all three NAC conditions "
          f"and the Arg90 contact cut\n")
    if len(pool) < args.n:
        sys.exit(f"FAIL: only {len(pool)} candidates; cannot select {args.n}")

    # maximin: repeatedly take the candidate furthest in time from everything chosen
    chosen = []
    ref = list(existing)
    for _ in range(args.n):
        best, bestd = None, -1
        for cand in pool:
            if cand[0] in [c[0] for c in chosen]:
                continue
            dmin = min(abs(cand[0] - r) for r in ref + [c[0] for c in chosen])
            if dmin > bestd:
                best, bestd = cand, dmin
        chosen.append(best)
        print(f"  selected frame {best[0]:>5} ({best[0]} ps)  "
              f"min separation {bestd} ps")
    chosen.sort()

    print(f"\n{'frame':>7}{'ps':>8}{'d (A)':>9}{'theta1':>9}{'dev':>7}"
          f"{'theta2':>9}{'dev':>7}{'Arg90':>9}")
    print("-" * 65)

    # The extra frames are written in the SAME eight-column layout as
    # selection_manifest.tsv, so step12b_extract.sh and step19b_collect.py read
    # them without modification. The angles are recorded alongside, in a
    # separate file, rather than changing the shared column order.
    with open(args.out, "w") as fh:
        fh.write("# idx\tframe\tprod_ps\tsite\tCHA_resid\tform_C6_C1\tr\tArg90_O13\n")
        for k, (i, d, t1, t2, c, brk) in enumerate(chosen, start=len(existing) + 1):
            print(f"{i:>7}{i:>8}{d:>9.3f}{t1:>9.1f}{abs(t1-TS_T1):>7.1f}"
                  f"{t2:>9.1f}{abs(t2-TS_T2):>7.1f}{c:>9.3f}")
            fh.write(f"{k}\t{i}\t{float(i):.1f}\t{args.site}\t{cha_resid}"
                     f"\t{d:.3f}\t{brk - d:.3f}\t{c:.3f}\n")

    ang = args.out.replace(".tsv", "_angles.tsv")
    with open(ang, "w") as fh:
        fh.write("frame\td_form\td_break\ttheta1\tdev1\ttheta2\tdev2\tArg90_O13\tNAC\n")
        for (i, d, t1, t2, c, brk) in chosen:
            fh.write(f"{i}\t{d:.3f}\t{brk:.3f}\t{t1:.1f}\t{abs(t1-TS_T1):.1f}"
                     f"\t{t2:.1f}\t{abs(t2-TS_T2):.1f}\t{c:.3f}\tyes\n")

    print(f"\nwrote {args.out}")
    print(f"wrote {ang}")
    print("\nTo bring these into the pipeline:")
    print("  1. append the rows of the new manifest to selection_manifest.tsv")
    print("  2. re-run step12b_extract.sh (re-extracting the existing twelve is")
    print("     idempotent - identical frames give identical rst7 files)")
    print("  3. bash step19a_ensemble_prepare.sh <frame> for each new frame")


if __name__ == "__main__":
    main()

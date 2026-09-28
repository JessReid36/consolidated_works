#!/usr/bin/env python3
"""
step11f_nac_check.py - verify the selected frames against the PUBLISHED
near-attack-conformer definition, not the distance criterion alone.

WHY THIS REPLACES THE EARLIER CHECK
-----------------------------------
Hur & Bruice (PNAS 2003, 100, 12015; doi:10.1073/pnas.1534873100) define a
chorismate NAC by THREE conditions, not one (their Fig. 1):

    (i)   d(C5-C16) <= 3.7 A          - van der Waals contact of the forming bond
    (ii)  |theta1 - 8.2 deg| <= 20    - approach angle vs the ring carbon's pi orbital
    (iii) |theta2 - 18.6 deg| <= 20   - side-chain pi direction vs the ring carbon

where 8.2 and 18.6 degrees are their transition-state values and NAC allows a
+/-20 degree deviation from them.

The step-12a frame filter tests condition (i) only. That matters because Hur &
Bruice specifically show a distance-only criterion admits mostly unreactive
structures: of 50,000 snapshots satisfying d <= 3.7 A in water, they observed
"only one snapshot of NAC", the rest being diequatorial conformers with the
pi-orbitals pointing away from one another (theta1 44 +/- 5, theta2 129 +/- 9).

Their 50,000-snapshot result is for chorismate in WATER, where NAC is ~1e-4 % of
the population. In the enzyme they report E.NAC as 30 % of the Michaelis complex,
so the risk is far lower here - but it must be shown, not assumed.

ANGLE DEFINITIONS
-----------------
Hur & Bruice label the forming-bond carbons C5 (in the ring) and C16 (the
terminal methylene of the enolpyruvyl side chain). In this work's convention
those are C6 (index 12, in the ring) and C1 (index 0, the methylene). Ring
membership is determined from connectivity rather than assumed.

    theta1 = angle between the ring carbon's pi axis and the vector from the
             ring carbon to the side-chain carbon
             ("C16 approaching angle relative to C5 pi orbital")

    theta2 = angle between the side-chain carbon's pi axis and the vector from
             the side-chain carbon to the ring carbon
             ("C16 pi direction relative to C5 position")

The pi axis at each sp2 carbon is the normal to the plane it forms with its two
nearest bonded neighbours. For the terminal methylene those neighbours are its
two hydrogens, which define the same plane as the heavy-atom framework.

VALIDATION
----------
Applied to this work's transition state the implementation gives theta1 = 12.9,
theta2 = 18.2 degrees, against Hur & Bruice's published TS values of 8.2 and
18.6. The theta2 agreement to 0.4 degrees, and theta1 to within 5 degrees, is
the check that the angles being measured are theirs and not something else.

Usage
-----
  python3 step11f_nac_check.py                    # all selected frames
  python3 step11f_nac_check.py --geometry ts.xyz  # a single structure
"""

import argparse
import os
import re
import struct
import sys
import numpy as np

COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66}
TS_THETA1, TS_THETA2 = 8.2, 18.6      # Hur & Bruice 2003, TS reference values
TOL = 20.0                             # their permitted deviation
D_NAC = 3.7                            # their van der Waals contact criterion
IDX_C1, IDX_C6 = 0, 12                 # this work's forming-bond atoms


def connectivity(els, xyz):
    n = len(els)
    D = np.sqrt(((xyz[:, None] - xyz[None]) ** 2).sum(-1))
    bond = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < (COV[els[i]] + COV[els[j]]) * 1.3:
                bond[i, j] = bond[j, i] = True
    return D, bond


def six_ring(els, bond):
    """Locate the cyclohexadiene ring by depth-first search over carbons."""
    found = []

    def dfs(path):
        if found:
            return
        if len(path) == 6:
            if bond[path[-1], path[0]]:
                found.extend(path)
            return
        for k in np.where(bond[path[-1]])[0]:
            if els[k] != "C" or k in path:
                continue
            dfs(path + [k])

    for s in range(len(els)):
        if els[s] == "C":
            dfs([s])
        if found:
            break
    return found


def pi_axis(c, xyz, bond, D):
    """Normal to the sp2 plane at atom c, from its two nearest bonded neighbours."""
    nb = sorted(np.where(bond[c])[0], key=lambda k: D[c, k])
    if len(nb) < 2:
        return None
    v = np.cross(xyz[nb[0]] - xyz[c], xyz[nb[1]] - xyz[c])
    nv = np.linalg.norm(v)
    return v / nv if nv > 0 else None


def nac_angles(els, xyz):
    """Return (d, theta1, theta2) for the forming bond, per Hur & Bruice."""
    D, bond = connectivity(els, xyz)
    ring = six_ring(els, bond)
    if not ring:
        return None
    c_ring = IDX_C6 if IDX_C6 in ring else IDX_C1
    c_side = IDX_C1 if c_ring == IDX_C6 else IDX_C6

    a_ring = pi_axis(c_ring, xyz, bond, D)
    a_side = pi_axis(c_side, xyz, bond, D)
    if a_ring is None or a_side is None:
        return None

    v = xyz[c_side] - xyz[c_ring]
    d = float(np.linalg.norm(v))
    u = v / d
    t1 = np.degrees(np.arccos(min(1.0, abs(float(np.dot(a_ring, u))))))
    t2 = np.degrees(np.arccos(min(1.0, abs(float(np.dot(a_side, -u))))))
    return d, float(t1), float(t2)


def verdict(d, t1, t2):
    ok_d = d <= D_NAC
    ok_1 = abs(t1 - TS_THETA1) <= TOL
    ok_2 = abs(t2 - TS_THETA2) <= TOL
    return ok_d, ok_1, ok_2, (ok_d and ok_1 and ok_2)


# ------------------------------------------------------------------ file input
def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0].split()[0])
    els, xyz = [], []
    for ln in lines[2:2 + n]:
        f = ln.split()
        els.append(f[0])
        xyz.append([float(f[1]), float(f[2]), float(f[3])])
    return els, np.asarray(xyz, float)


# ------------------------------------------------------- trajectory input
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


def main():
    ap = argparse.ArgumentParser()
    home = os.path.expanduser("~")
    ap.add_argument("--geometry", default=None,
                    help="check a single .xyz instead of the trajectory frames")
    ap.add_argument("--prmtop",
                    default=f"{home}/system_development/03_amber/tleap_build/complex_solvated.prmtop")
    ap.add_argument("--traj",
                    default=f"{home}/system_development/04_amber_md/10c_production/prod.nc")
    ap.add_argument("--manifest",
                    default=f"{home}/system_development/05_qmmm/12_frame_selection/selection_manifest.tsv")
    ap.add_argument("--site", default="CHA#2")
    ap.add_argument("--validate", default=None,
                    help="a TS .xyz; prints the comparison with the published TS angles")
    args = ap.parse_args()

    print("NAC verification against Hur & Bruice (PNAS 2003, 100, 12015)")
    print(f"  criteria: d <= {D_NAC} A, |theta1 - {TS_THETA1}| <= {TOL}, "
          f"|theta2 - {TS_THETA2}| <= {TOL}\n")

    if args.validate:
        els, xyz = read_xyz(args.validate)
        r = nac_angles(els, xyz)
        if r:
            d, t1, t2 = r
            print("IMPLEMENTATION VALIDATION against the published TS geometry")
            print(f"  this work's TS : theta1 {t1:.1f}   theta2 {t2:.1f}")
            print(f"  Hur & Bruice   : theta1 {TS_THETA1}    theta2 {TS_THETA2}")
            print(f"  difference     : {abs(t1-TS_THETA1):.1f}       {abs(t2-TS_THETA2):.1f}\n")

    if args.geometry:
        els, xyz = read_xyz(args.geometry)
        r = nac_angles(els, xyz)
        if not r:
            sys.exit("FAIL: could not locate the ring or a pi axis")
        d, t1, t2 = r
        od, o1, o2, ok = verdict(d, t1, t2)
        print(f"{args.geometry}")
        print(f"  d      {d:6.3f} A   {'PASS' if od else 'FAIL'}")
        print(f"  theta1 {t1:6.1f}     dev {abs(t1-TS_THETA1):5.1f}   {'PASS' if o1 else 'FAIL'}")
        print(f"  theta2 {t2:6.1f}     dev {abs(t2-TS_THETA2):5.1f}   {'PASS' if o2 else 'FAIL'}")
        print(f"  => {'NAC' if ok else 'NOT a NAC by the full definition'}")
        return

    # ------------------------------------------------- trajectory frames
    natom, names, labels, resptr = parse_prmtop(args.prmtop)
    nres = len(labels)
    starts = [p - 1 for p in resptr] + [natom]
    cha = []
    for ri in range(nres):
        if labels[ri] == "CHA":
            cha.append((ri, {names[a]: a for a in range(starts[ri], starts[ri + 1])}))
    site = int(args.site.split("#")[1]) - 1
    ri, amap = cha[site]
    order = sorted(amap.items(), key=lambda kv: kv[1])
    sub_idx = [a for _, a in order]
    sub_els = [("H" if nm[0] == "H" else nm[0]) for nm, _ in order]

    H = nc_open(args.traj)
    cvar = next(v for v in H["rec"] if v["name"] == "coordinates")

    def frame(i):
        o = cvar["begin"] + i * H["recsize"]
        full = np.array(np.memmap(args.traj, dtype=">f4", mode="r", offset=o,
                                  shape=(natom, 3))[:max(sub_idx) + 1], float)
        return full[sub_idx]

    rows = [l.split("\t") for l in open(args.manifest).read().splitlines()
            if not l.startswith("#") and l.strip()]

    print(f"{'frame':>7}{'d (A)':>9}{'theta1':>9}{'dev':>7}{'theta2':>9}{'dev':>7}"
          f"{'verdict':>10}")
    print("-" * 58)
    npass = 0
    for r in rows:
        fr = int(r[1])
        xyz = frame(fr)
        res = nac_angles(sub_els, xyz)
        if not res:
            print(f"{fr:>7}   could not determine geometry")
            continue
        d, t1, t2 = res
        od, o1, o2, ok = verdict(d, t1, t2)
        npass += ok
        print(f"{fr:>7}{d:>9.3f}{t1:>9.1f}{abs(t1-TS_THETA1):>7.1f}"
              f"{t2:>9.1f}{abs(t2-TS_THETA2):>7.1f}{'NAC' if ok else 'distance only':>10}")

    print(f"\n  {npass} of {len(rows)} frames satisfy the full three-condition definition")
    if npass == len(rows):
        print("  All selected frames are NACs by the published criterion, so the")
        print("  distance-only filter used at selection did not admit unreactive")
        print("  conformers in this case. State this as a verified property rather")
        print("  than assuming the distance criterion sufficed.")
    else:
        print("  Not all frames satisfy the full definition. The distance-only filter")
        print("  has admitted structures that are not near-attack conformers, and the")
        print("  affected frames should be reported explicitly or excluded.")


if __name__ == "__main__":
    main()

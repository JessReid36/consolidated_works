#!/usr/bin/env python3
"""
step18c_extract_reduced_region.py

Build the reduced ActiveAtoms region for the transition-state optimisation
(step 18c). The region is the chorismate substrate plus the sidechains of the
first-shell catalytic residues that line the inter-subunit active site.

The active site is inter-subunit, so the correct copy of each catalytic residue
is identified by proximity to the substrate (CHA#2, residue 383) rather than by
sequence number: the trimer contains three copies of each residue and only the
copy nearest the bound substrate is catalytically relevant. Sidechain atoms are
kept and backbone atoms dropped, since only the sidechain tips (guanidinium,
carboxylate, hydroxyl) contribute to transition-state stabilisation and need to
relax during the saddle search.

Output: active_reduced.txt, a space-separated list of 0-based atom indices for
ORCA's ActiveAtoms selection (24 substrate atoms + first-shell sidechains).

cpptraj is unavailable on this cluster (missing libreadline.so.6), so the prmtop
is parsed natively.

Usage:
    python3 step18c_extract_reduced_region.py complex_solvated.prmtop reactant_opt.pdb
"""

import sys
import re
import numpy as np
from collections import defaultdict

# CHA#2 substrate is residue 383 (1-based) -> index 382 (0-based)
CHA_RES_INDEX = 382
# substrate QM atoms, 0-based, from the cha_gaff ordering
CHA_ATOMS = list(range(6207, 6231))

# first-shell catalytic residues, identified by proximity below (not hard-coded here);
# Lys60 is anticatalytic and is deliberately excluded from the movable set.
# Backbone atom names dropped when extracting sidechains:
BACKBONE = {"N", "H", "CA", "HA", "C", "O", "H1", "H2", "H3", "OXT", "HA2", "HA3"}


def read_flag(txt, flag):
    m = re.search(r"%FLAG " + flag + r"\s*\n%FORMAT\((.*?)\)\s*\n(.*?)(?=%FLAG|\Z)", txt, re.S)
    return m.group(2) if m else ""


def parse_prmtop(prmtop_path):
    txt = open(prmtop_path).read()
    nm = read_flag(txt, "ATOM_NAME").replace("\n", "")
    names = [nm[i:i + 4].strip() for i in range(0, len(nm), 4)]
    rl = read_flag(txt, "RESIDUE_LABEL").replace("\n", "")
    reslabels = [rl[i:i + 4].strip() for i in range(0, len(rl), 4)]
    resptr = [int(x) for x in read_flag(txt, "RESIDUE_POINTER").split()]
    resptr.append(len(names) + 1)  # sentinel so the last residue has an end
    return names, reslabels, resptr


def res_atom_range(resptr, i):
    """0-based atom indices belonging to residue i (0-based)."""
    return range(resptr[i] - 1, resptr[i + 1] - 1)


def read_pdb_coords(pdb_path, natom):
    recs = [l for l in open(pdb_path) if l[:6] in ("ATOM  ", "HETATM")]
    if len(recs) != natom:
        raise ValueError(f"PDB atoms {len(recs)} != prmtop {natom}")
    return np.array([[float(r[30:38]), float(r[38:46]), float(r[46:54])] for r in recs])


def main():
    prmtop_path, pdb_path = sys.argv[1], sys.argv[2]
    names, reslabels, resptr = parse_prmtop(prmtop_path)
    natom = len(names)
    xyz = read_pdb_coords(pdb_path, natom)

    cha_xyz = xyz[CHA_ATOMS]

    def min_dist_to_cha(ridx):
        ats = list(res_atom_range(resptr, ridx))
        return np.min(np.linalg.norm(xyz[ats][:, None, :] - cha_xyz[None, :, :], axis=2))

    bylabel = defaultdict(list)
    for i, lab in enumerate(reslabels):
        bylabel[lab].append(i)

    # identify the closest copy of each catalytic residue type to the substrate.
    # For this system the first shell resolves to (1-based residue numbers):
    #   ARG 134, 63, 217  (Arg90, Arg7, Arg63'); GLU 205 (Glu78); TYR 235 (Tyr108)
    # Lys60 (LYS 60) is the closest lysine but is anticatalytic -> excluded.
    chosen = []
    for lab, n_keep in [("ARG", 3), ("GLU", 1), ("TYR", 1)]:
        cand = sorted((min_dist_to_cha(i), i) for i in bylabel.get(lab, []))
        chosen += [i for _, i in cand[:n_keep]]

    sidechain = []
    for i in chosen:
        for a in res_atom_range(resptr, i):
            if names[a] not in BACKBONE:
                sidechain.append(a)

    active = sorted(set(CHA_ATOMS + sidechain))
    with open("active_reduced.txt", "w") as fh:
        fh.write(" ".join(map(str, active)))

    print(f"substrate atoms:      {len(CHA_ATOMS)}")
    print(f"catalytic sidechains: {len(sidechain)}")
    print(f"total movable:        {len(active)}")
    print("written to active_reduced.txt")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
renumber_ext_manifest.py - put frames selected from the production extension onto
the continuous production timeline.

WHY
    The production run was extended in a separate job, so the trajectory is in
    two files. Frames selected from the extension are numbered from zero within
    their own file, which collides with the numbering of the first run: frame
    13320 of the extension is a different structure from frame 13320 of prod.nc,
    and nothing in the manifest or in the extracted filename would say which is
    which.

    The two files are one trajectory, since the extension restarted from the
    first run's restart file with velocities. Adding the length of the first
    trajectory to the frame index therefore gives the frame's position in the run
    as a whole, which is both unambiguous and physically correct.

WHAT IT CHANGES
    The frame and prod_ps columns, by the length of the first trajectory. The idx
    column is renumbered to continue after the existing manifest, if one is
    given. Nothing else is touched.

    The offset is read from the first trajectory's netCDF header rather than
    assumed, so it stays correct if either file is regenerated.

NOTE ON TIME
    The reported simulation time carries a 100 ps offset from the five
    equilibration stages, which ran 10000 steps at 2 fs each and were continued
    with irest=1. The trajectory's own frames are unaffected: frame 0 of prod.nc
    is the first production frame however the clock is labelled. The prod_ps
    column here follows the frame index, not the simulation clock, and that
    convention is preserved.

Usage:
    python3 renumber_ext_manifest.py in.tsv out.tsv [--existing selection_manifest.tsv]
"""
import argparse
import os
import struct
import sys


def nc_frames(path):
    """Number of records in a classic or 64-bit-offset netCDF file.

    Read from the header rather than inferred from the file size, since a
    trajectory may be written with a different record layout.
    """
    data = open(path, "rb").read(65536)
    if data[:3] != b"CDF":
        raise ValueError("%s is not a netCDF classic file" % path)
    off64 = data[3] == 2
    pos = [4]

    def u32():
        v = struct.unpack_from(">I", data, pos[0])[0]
        pos[0] += 4
        return v

    def u64():
        v = struct.unpack_from(">Q" if off64 else ">I", data, pos[0])[0]
        pos[0] += 8 if off64 else 4
        return v

    def name():
        n = u32()
        s = data[pos[0]:pos[0] + n].decode()
        pos[0] += n + ((4 - n % 4) % 4)
        return s

    numrecs = u32()
    # dimensions
    tag = u32(); ndim = u32()
    dims = []
    if tag == 10:
        for _ in range(ndim):
            dims.append((name(), u32()))
    # global attributes, skipped
    tag = u32(); natt = u32()
    if tag == 12:
        for _ in range(natt):
            name()
            nc_type = u32(); n = u32()
            size = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}[nc_type]
            pos[0] += n * size + ((4 - (n * size) % 4) % 4)
    # variables
    tag = u32(); nvar = u32()
    recsize = 0
    begins = []
    if tag == 11:
        for _ in range(nvar):
            name()
            nd = u32()
            shape = [dims[u32()][1] for _ in range(nd)]
            tag2 = u32(); na = u32()
            if tag2 == 12:
                for _ in range(na):
                    name()
                    nc_type = u32(); n = u32()
                    size = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}[nc_type]
                    pos[0] += n * size + ((4 - (n * size) % 4) % 4)
            nc_type = u32(); vsize = u32(); begin = u64()
            if shape and shape[0] == 0:          # a record variable
                recsize += vsize
                begins.append(begin)
    if numrecs != 0xFFFFFFFF:
        return numrecs
    if not begins or recsize == 0:
        raise ValueError("%s: cannot determine the record count" % path)
    return (os.path.getsize(path) - min(begins)) // recsize


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("outfile")
    home = os.path.expanduser("~")
    ap.add_argument("--traj1",
                    default=f"{home}/system_development/04_amber_md/"
                            f"10c_production/prod.nc",
                    help="the first trajectory, whose length is the offset")
    ap.add_argument("--existing",
                    default=f"{home}/system_development/05_qmmm/"
                            f"12_frame_selection/selection_manifest.tsv",
                    help="if present, idx continues after the last row of this")
    args = ap.parse_args()

    off = nc_frames(args.traj1)
    print(f"offset from {os.path.basename(args.traj1)}: {off} frames")

    start_idx = 0
    if os.path.exists(args.existing):
        idxs = [int(l.split("\t")[0]) for l in open(args.existing)
                if l.strip() and not l.startswith("#")]
        if idxs:
            start_idx = max(idxs)
            print(f"existing manifest ends at idx {start_idx}")

    lines = open(args.infile).read().splitlines()
    out = []
    n = 0
    for l in lines:
        if not l.strip():
            continue
        if l.startswith("#"):
            out.append(l)
            continue
        c = l.split("\t")
        if len(c) < 3:
            out.append(l)
            continue
        n += 1
        c[0] = str(start_idx + n)
        old_fr = int(c[1])
        c[1] = str(old_fr + off)
        c[2] = "%.1f" % (float(c[2]) + off)
        out.append("\t".join(c))
        if n <= 3 or n == len(lines) - 1:
            print(f"  frame {old_fr:>6} -> {c[1]:>6}")

    with open(args.outfile, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print(f"\nwrote {args.outfile}: {n} frames renumbered")
    print("\nAppend these rows to the selection manifest, then run step12b_extract.sh.")
    print("The patched extractor reads frames below %d from the first" % off)
    print("trajectory and the rest from the extension.")


if __name__ == "__main__":
    main()

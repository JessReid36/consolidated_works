#!/usr/bin/env python3
"""
verify_all.py - check the whole of consolidated_works in one pass.

Run this after any extraction, any pull, or any incident. It answers three questions:

  1. Does every file match the SHA-256 recorded for it in its step's manifest?
     This is what catches a truncated download, a partial extraction, or a file that was
     silently replaced. It is the check that would have caught the reactant_reduced.out
     truncation before it was committed.
  2. Does every step's own check script still pass?
  3. Is the bookkeeping consistent - does INDEX.tsv list every step folder and no others?

Exit status 0 if everything passes, 1 otherwise.

Usage:
    python3 verify_all.py            # manifests and bookkeeping only, fast
    python3 verify_all.py --checks   # also re-run every step's check script
"""
import hashlib
import re
import subprocess
import sys
from pathlib import Path

# Manifest sections whose lines are prose or pointers, not files to hash.
NON_FILE_SECTIONS = {
    "depends", "declared", "toolchain", "note", "oversized,", "held", "literature",
    "key", "the", "NOT", "THE", "DO", "incomplete", "SIGN", "environment",
}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def resolve(step, section, recorded_path):
    """Find the local file a manifest line refers to, allowing for renames."""
    name = Path(recorded_path).name
    if section == "scripts" and name == "README.md":
        name = "local_workstation_README.md"
    if section == "stage_manifests" and name == "selection_manifest.tsv":
        name = "step12a_selection_manifest.tsv"
    candidates = [
        step / section / name,
        step / "outputs" / name,
        step / "outputs" / "stage_manifests" / name,
    ]
    # step 10 records outputs by their path under 04_amber_md/
    if recorded_path.startswith("04_amber_md/"):
        candidates.insert(0, step / "outputs" / recorded_path[len("04_amber_md/"):])
    for c in candidates:
        if c.exists():
            return c
    return candidates[0]


def main():
    root = Path(__file__).resolve().parent
    run_checks = "--checks" in sys.argv
    steps = sorted(p for p in root.iterdir()
                   if p.is_dir() and (p / "manifest.txt").exists())

    print(f"consolidated_works verification — {len(steps)} step folders\n")

    print("1. file integrity against the manifests")
    total_ok = total_bad = total_missing = 0
    for step in steps:
        ok = bad = missing = 0
        section = None
        for line in (step / "manifest.txt").read_text().splitlines():
            if line.startswith("## "):
                section = line[3:].split()[0]
                continue
            m = re.match(r"\s+(\S+)\s+(\d+)\s+([0-9a-f]{64})$", line)
            if not m or section is None or section in NON_FILE_SECTIONS:
                continue
            path, size, digest = m.group(1), int(m.group(2)), m.group(3)
            f = resolve(step, section, path)
            if not f.exists():
                missing += 1
                print(f"    MISSING  {step.name}/{section}/{Path(path).name}")
            elif f.stat().st_size != size:
                bad += 1
                print(f"    SIZE     {step.name}/{section}/{Path(path).name}  "
                      f"{f.stat().st_size} against {size} recorded  "
                      f"<-- TRUNCATED OR REPLACED")
            elif sha(f) != digest:
                bad += 1
                print(f"    HASH     {step.name}/{section}/{Path(path).name}  "
                      f"<-- CONTENT CHANGED")
            else:
                ok += 1
        total_ok += ok
        total_bad += bad
        total_missing += missing
        flag = "" if not (bad or missing) else "   <-- PROBLEM"
        print(f"    {step.name:<28} {ok:>4} ok"
              f"{'' if not bad else f', {bad} bad'}"
              f"{'' if not missing else f', {missing} missing'}{flag}")
    print(f"\n    total: {total_ok} verified, {total_bad} bad, {total_missing} missing")

    print("\n2. bookkeeping")
    index = root / "INDEX.tsv"
    listed = set()
    if index.exists():
        for line in index.read_text().splitlines()[1:]:
            if line.strip():
                listed.add(line.split("\t")[0])
    present = {s.name for s in steps}
    only_index = sorted(listed - present)
    only_disk = sorted(present - listed)
    print(f"    INDEX.tsv lists {len(listed)}, {len(present)} step folders on disk")
    if only_index:
        print(f"    listed but absent: {only_index}")
    if only_disk:
        print(f"    present but unlisted: {only_disk}")
    dupes = [s for s in listed if
             sum(1 for l in index.read_text().splitlines()[1:]
                 if l.startswith(s + "\t")) > 1]
    if dupes:
        print(f"    DUPLICATE INDEX ROWS: {sorted(set(dupes))}")

    # Open-item numbering is easy to collide when several sessions append.
    oi = root / "00_crosscutting" / "open_items.txt"
    if oi.exists():
        nums = re.findall(r"^ITEM (\d+) ", oi.read_text(), re.M)
        seen, coll = set(), []
        for n in nums:
            if n in seen:
                coll.append(n)
            seen.add(n)
        print(f"    open_items.txt: {len(nums)} items"
              + (f", DUPLICATE NUMBERS {sorted(set(coll))}" if coll else ", numbering clean"))

    rc = 1 if (total_bad or total_missing or only_index or only_disk or dupes) else 0

    if run_checks:
        print("\n3. re-running each step's check script")
        for step in steps:
            scripts = sorted((step / "checks").glob("*verify*.py")) \
                if (step / "checks").exists() else []
            for s in scripts:
                r = subprocess.run([sys.executable, s.name], cwd=s.parent,
                                   capture_output=True, text=True)
                tail = [l for l in r.stdout.splitlines() if "PASS" in l or "FAIL" in l]
                summary = tail[-1] if tail else "(no summary)"
                mark = "ok  " if r.returncode == 0 else "FAIL"
                print(f"    [{mark}] {step.name}/{s.name}: {summary.strip()[:80]}")
                if r.returncode != 0:
                    rc = 1
    else:
        print("\n3. step check scripts not run (pass --checks to include them)")

    print(f"\n{'ALL VERIFIED' if rc == 0 else 'PROBLEMS FOUND — see above'}")
    return rc


if __name__ == "__main__":
    sys.exit(main())

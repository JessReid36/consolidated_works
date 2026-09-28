#!/usr/bin/env python3
"""
build_consolidated.py - assemble consolidated_works/ from the source repos.

Each step folder gets:
    scripts/        the script(s) that produced the outputs, copied verbatim
    outputs/        the committed outputs, or a .pointer file if oversized
    manifest.txt    repo / branch / commit / path / sha256 for every file
    interpretation.txt   fixed-skeleton note, stub unless already written

Nothing here is a source of truth. The manifest pins every copied file back
to its canonical home so duplicates can be culled mechanically later.

Usage:
    python3 build_consolidated.py \
        --code   ~/chorismate-thesis-code \
        --results ~/chorismate-thesis-results \
        --dest   ~/consolidated_works \
        [--steps 01_inputs_acquisition 02_terminal_repair ...]

Re-running is safe: existing interpretation.txt files are never overwritten.
"""
import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

# Files larger than this are recorded as a pointer, not copied.
# The results repo already keeps multi-hundred-MB trajectories out of git;
# this keeps consolidated_works clonable.
SIZE_LIMIT = 10 * 1024 * 1024  # 10 MB

# step id -> (title, [script specs], [output globs])
# script spec: (repo_key, path_within_repo)
# output glob: (repo_key, glob_within_repo)
STEPS = {
    "01_inputs_acquisition": (
        "Acquisition and provenance locking of raw inputs",
        [("code", "phase1_system_dev/step01_download_inputs.sh"),
         ("code", "phase1_system_dev/step01b_ligand_tsa_registration.sh")],
        [("results", "01_inputs/structures/*.pdb"),
         ("results", "01_inputs/ligands/*.mol2"),
         ("results", "00_admin/step01_download_manifest.tsv"),
         ("results", "00_admin/step01_provenance.txt"),
         ("results", "00_admin/sha256_step01_inputs.txt"),
         ("results", "00_admin/ligand_tsa_registration.tsv")],
    ),
    "02_terminal_repair": (
        "Terminal repair of the A/B/C trimer against 1DBF",
        [("code", "phase1_system_dev/step02_terminal_repair_abc.sh")],
        [("results", "02_preparation/protein_only/2cht_abc_raw_protein.pdb"),
         ("results", "02_preparation/protein_only/1dbf_abc_aligned_to_2cht.pdb"),
         ("results", "02_preparation/protein_only/abc_repaired.pdb"),
         ("results", "02_preparation/protein_only/abc_repaired_unadjusted.pdb"),
         ("results", "00_admin/step02_*")],
    ),
    "03_cleanup_validate": (
        "Hydrogen strip, altloc resolution and scaffold validation",
        [("code", "phase1_system_dev/step03_cleanup_validate.sh")],
        [("results", "02_preparation/protein_only/abc_repaired_clean.pdb"),
         ("results", "00_admin/step03_*")],
    ),
    "04_substrate_placement": (
        "Chorismate placement and derivation of the contact set",
        [("code", "phase1_system_dev/step04_place_chorismate_abc.sh"),
         ("code", "phase1_system_dev/step04b_contact_criterion.py")],
        [("results", "02_preparation/ligand_placement/*"),
         ("results", "00_admin/step04_*"),
         ("results", "00_admin/step04b_*")],
    ),
    "05_arg63_provenance": (
        "Arg63' rotamer provenance across the twelve crystallographic sites",
        [("code", "phase1_system_dev/step05_arg63_provenance.sh"),
         ("code", "phase1_system_dev/step05b_arg63_conformation_check.sh")],
        [("results", "00_admin/step05_*"),
         ("results", "00_admin/step05b_*")],
    ),
    "06_acceptance_freeze": (
        "Acceptance freeze of the pre-protonation complex",
        [("code", "phase1_system_dev/step06_freeze_preprotonation.sh")],
        [("results", "02_preparation/accepted_preprotonation/*"),
         ("results", "00_admin/step06_*"),
         ("results", "00_admin/sha256_step06_accepted_preprotonation.txt")],
    ),
    "07_protonation_hpp": (
        "H++ protonation-state assignment and acceptance audit",
        [("code", "phase1_system_dev/step07a_prepare_hpp_input.sh"),
         ("code", "phase1_system_dev/step07b_audit_repair_hpp.sh")],
        [("results", "03_amber/protonation_hpp/*"),
         ("results", "00_admin/step07a_*"),
         ("results", "00_admin/step07b_*")],
    ),
}

SKELETON = """\
# {step_id} — {title}

## What this step did
[TODO]

## Key numbers
[TODO]

## Acceptance threshold and verdict
[TODO — state the threshold, the measured value, and pass/fail]

## Literature comparison
[TODO — which paper, which quantity, agreement or disagreement. If none, say none.]

## Caveats and known issues
[TODO]

## Cited by
[TODO — write-up subsection(s) that reference this step. Blank means nothing
cites it yet, which is the signal for the cull.]
"""


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_info(repo):
    def run(*a):
        return subprocess.run(["git", "-C", str(repo), *a],
                              capture_output=True, text=True).stdout.strip()
    return run("rev-parse", "HEAD"), run("rev-parse", "--abbrev-ref", "HEAD")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--code", required=True, type=Path)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--dest", required=True, type=Path)
    ap.add_argument("--steps", nargs="*", default=None)
    args = ap.parse_args()

    repos = {"code": args.code, "results": args.results}
    for key, path in repos.items():
        if not (path / ".git").exists():
            sys.exit(f"not a git repo: {path}")
    info = {k: git_info(v) for k, v in repos.items()}
    for k, (commit, branch) in info.items():
        print(f"  {k:8s} {repos[k]}  {branch}  {commit[:7]}")

    wanted = args.steps or list(STEPS)
    for step_id in wanted:
        if step_id not in STEPS:
            sys.exit(f"unknown step: {step_id}")
        title, scripts, outputs = STEPS[step_id]
        root = args.dest / step_id
        (root / "scripts").mkdir(parents=True, exist_ok=True)
        (root / "outputs").mkdir(parents=True, exist_ok=True)

        lines = [f"step: {step_id}", f"title: {title}", ""]
        for kind, specs in (("script", scripts), ("output", outputs)):
            lines.append(f"## {kind}s")
            for repo_key, pattern in specs:
                repo = repos[repo_key]
                commit, branch = info[repo_key]
                matches = sorted(repo.glob(pattern))
                if not matches:
                    lines.append(f"  MISSING  {repo_key}:{pattern}"
                                 f"  (not present at {commit[:7]})")
                    print(f"  ! {step_id}: MISSING {repo_key}:{pattern}")
                    continue
                for src in matches:
                    if src.is_dir():
                        continue
                    rel = src.relative_to(repo)
                    dest_dir = root / ("scripts" if kind == "script" else "outputs")
                    dest = dest_dir / src.name
                    size = src.stat().st_size
                    digest = sha256(src)
                    if size > SIZE_LIMIT:
                        dest.with_suffix(dest.suffix + ".pointer").write_text(
                            f"canonical_repo: {repo_key}\n"
                            f"canonical_path: {rel}\n"
                            f"commit: {commit}\n"
                            f"size_bytes: {size}\n"
                            f"sha256: {digest}\n"
                            f"note: exceeds {SIZE_LIMIT} byte limit; not copied\n")
                        lines.append(f"  POINTER  {rel}  {size}  {digest}")
                    else:
                        shutil.copy2(src, dest)
                        lines.append(f"  {rel}  {size}  {digest}")
                    lines.append(f"    repo={repo_key} branch={branch} commit={commit}")
            lines.append("")

        (root / "manifest.txt").write_text("\n".join(lines))
        interp = root / "interpretation.txt"
        if not interp.exists():
            interp.write_text(SKELETON.format(step_id=step_id, title=title))
        print(f"  built {step_id}")

    print("\nDone. Review manifests, then commit.")


if __name__ == "__main__":
    main()

# consolidated_works

Assembled, per-step collection of the scripts, outputs and interpretation for the
chorismate mutase Phase 1 work.

**Nothing in this repository is a source of truth.** Canonical homes are:

- scripts  -> `JessReid36/chorismate-thesis-code`
- outputs  -> `JessReid36/chorismate-thesis-results`

Every file copied here is pinned back to its canonical home by repo, branch, commit and
SHA-256 in that step's `manifest.txt`. Duplicates are therefore safe: the cull is a
mechanical diff against the manifests, not a judgement call.

## Layout

    NN_step_name/
      scripts/            the script(s) that produced the outputs, copied verbatim
      outputs/            the committed outputs, or a .pointer file if oversized
      manifest.txt        repo / branch / commit / path / sha256 for every file
      interpretation.txt  fixed-skeleton note; see below

Files over 10 MB are not copied. They are recorded as `<name>.pointer` carrying the
canonical path, commit, size and SHA-256. The results repo already keeps multi-hundred-MB
trajectories out of git under `CHECKSUMS_large_files.txt`; this preserves that policy.

## interpretation.txt skeleton

Fixed headings, so that twenty-five of these written on different days remain consistent
and can be drafted from directly:

    What this step did
    Key numbers
    Acceptance threshold and verdict
    Literature comparison
    Caveats and known issues
    Cited by

`Cited by` is load-bearing. An empty `Cited by` is the signal that a step folder is a
candidate for the cull.

## Rebuilding

    python3 build_consolidated.py \
        --code    ~/chorismate-thesis-code \
        --results ~/chorismate-thesis-results \
        --dest    ~/consolidated_works

Re-running is safe: existing `interpretation.txt` files are never overwritten. Missing
sources are reported as `MISSING` in the manifest rather than silently skipped.

## Status

Steps 01-07 built and interpreted (2026-09-28), from code `tier1-realism` f0e059c and
results `main` aa131d5. Steps 08 onward not yet added.

Two known MISSING entries in `04_substrate_placement`:
`step04b_contact_criterion.py` and `step04b_contact_criterion_report.tsv` are not yet
committed to their canonical repos.

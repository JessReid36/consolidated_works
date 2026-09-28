# consolidated_works

Per-step collection of the scripts, outputs, verification and interpretation for the
chorismate mutase Phase 1 work, built one step at a time.

**Nothing here is a source of truth.** Canonical homes:

- scripts -> `JessReid36/chorismate-thesis-code`
- outputs -> `JessReid36/chorismate-thesis-results`

Every file is pinned back by repo, branch, commit and SHA-256 in that step's
`manifest.txt`, so duplicates are safe and the eventual cull is a mechanical diff
rather than a judgement call.

## Layout

    00_crosscutting/        claims that span steps, written as encountered
                            (comparator_decision, literature_position,
                             protocol_deviations, sourcing_gaps, open_items)
    INDEX.tsv               step -> subsection -> scripts -> outputs -> check status
    NN_step_name/
      scripts/              from the code repo, pinned
      scripts_as_run/       copies archived in results 00_admin; these CAN diverge
      outputs/              from the results repo, pinned
      checks/               re-runnable verification script + recorded output
      manifest.txt          provenance and checksums for every file, plus check environment
      interpretation.txt    fixed skeleton, see below
      writeup/              this step's Methods subsection, .md source + .docx

Files over 10 MB are not copied; they are recorded as `<name>.pointer` carrying canonical
path, commit, size and SHA-256, preserving the existing large-file policy in
`CHECKSUMS_large_files.txt`.

## Method, per step

1. Grep the script for what it actually writes. Never infer outputs from filenames.
2. Locate each declared output. Report ABSENT for declared-but-missing, UNDECLARED for
   present-but-unclaimed.
3. Read every output's contents.
4. Re-derive what is derivable; that becomes `checks/`.
5. Literature check: what did Claeyssens or Agbaglo do here, do we deviate, is it recorded.
6. Objective check: does this step's choice constrain or bias the Phase 2 charge design.
7. Write `interpretation.txt`, commit that folder alone, move on.

## interpretation.txt skeleton

    What this step did
    Key numbers
    Acceptance threshold and verdict
    Literature comparison
    Caveats and known issues
    Cited by

`Cited by` is load-bearing: an empty one flags a cull candidate, and the set of all of them
should reconcile against the write-up's reference list.

## Status

01 to 03 built, verified and drafted 2026-09-28 (code tier1-realism f0e059c, results main aa131d5).
Steps 04 onward not yet added.

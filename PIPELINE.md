# Pipeline map

Which consolidated_works folder covers which pipeline step, and which thesis
subsection it feeds. Read this before looking for a step by number: three folders
cover more than one step, and one step number appears in no folder name.

| Folder | Pipeline steps | Subsection |
|---|---|---|
| `01_inputs_acquisition` | 01, 01b | X.1.1 |
| `02_terminal_repair` | 02 | X.1.2 |
| `03_cleanup_validate` | 03 | X.1.3 |
| `04_substrate_placement` | 04 | X.1.4 |
| `05_arg63_provenance` | 05, 05b | X.1.4 (Arg63 paragraph) |
| `06_acceptance_freeze` | 06 | X.1.5 (freeze half) |
| `07_protonation_hpp` | 07a, 07b | X.1.5 |
| `08_ligand_gaff` | 08a, 08b | X.1.6 |
| `09_tleap_build` | 09a, 09b | X.1.7 |
| `10_md` | 10, 10a–10d | X.2 |
| `11_trajectory_analysis` | 11a–11f | X.2 (points to `10_md`) |
| `12_frame_selection` | 12a, 12b | X.2 (points to `10_md`) |
| `13_qmmm_setup` | **13a, 13b and 14** | X.3 |
| `15_17_path_endpoints_scan` | 15, 16, 17 | X.3 |
| `18_nebts_ts_irc` | 18, 18c–18f | X.3 and Results |
| `19_ensemble` | 19a–19e | Results |

## Two things the folder names do not say

**There is no folder named for step 14.** `13_qmmm_setup` covers steps 13a, 13b and 14
together: the bridge, the smoke test and the DFT single point are one operation and
neither half is interpretable alone. A folder `13_14_qmmm_bridge` existed briefly as a
duplicate and was removed; if a manifest or note still references it, that reference is
stale.

**Three folders cover ranges.** `13_qmmm_setup`, `15_17_path_endpoints_scan` and
`10_md` each cover several pipeline steps. The Subsection column, not the folder name,
is what maps to the thesis.

## Where the write-up lives

One subsection may draw on several steps. Where that happens, the prose lives in ONE
folder and the others hold a `writeup/README.txt` pointing to it.

| Subsection | Drafted in |
|---|---|
| X.1.1 – X.1.4 | the step's own folder |
| X.1.5 | `07_protonation_hpp/writeup/` (06 points there) |
| X.1.6 | `08_ligand_gaff/writeup/` |
| X.1.7 | `09_tleap_build/writeup/` |
| X.2 | `10_md/writeup/` (11 and 12 point there) |
| X.3 | not yet drafted |
| Results | not yet drafted |

## Verifying the repository

    python3 verify_all.py            # file integrity and bookkeeping
    python3 verify_all.py --checks   # also re-runs every step's check script

Expect one known failure at `11_trajectory_analysis`: two documentation discrepancies
in `step11_trajectory_analysis.md`, retained as failures until that file is corrected.

# X.1.2 Terminal repair of the A/B/C trimer

> **Drafting note (not part of the thesis text).** Follows the written specification in
> `WRITEUP_PROMPT.md` §2–§3; the house-style exemplar is not held. Source files for every
> statement below are in `consolidated_works/02_terminal_repair/`.

---

The 2CHT deposition resolves twelve chains, of which I retained A, B and C as the
biological trimer. Chains A, B and C each lack ordered residues at one or both termini, so
I repaired them from PDB entry 1DBF, which resolves the same protein at higher resolution
and without a bound ligand (Ladner *et al.*, 2000). `step02_terminal_repair_abc.sh`
extracted the A/B/C protein from `2cht_raw.pdb`, superposed 1DBF globally onto that trimer
over common backbone atoms, and transferred the unresolved terminal residues from the
aligned donor into the 2CHT frame. The superposition used 1396 common backbone atoms and
gave a root-mean-square deviation of 0.8037 Å.

I transferred residue 1 and residues 116 to 127 for chain A, residue 1 and residues 115 to
127 for chain B, and residues 1 to 2 and residues 115 to 127 for chain C. These ranges matched those used by
Agbaglo *et al.* (2024) in their preparation of the same deposition from the same donor.

Two grafted termini did not meet the existing backbone at a physical peptide C–N distance,
and I closed both by rigid translation of the terminal block to the canonical value of
1.329 Å, which conditions the join geometrically without altering the internal geometry of
the transferred block. The script permitted this closure only at a hard-coded set of
graft joins and treated any other abnormal join as a fatal error, preventing backbone
remote from the graft being moved without notice.

**Table X.1.2a.** Residues transferred from aligned 1DBF, per chain. Source:
`outputs/step02_residue_source_report.tsv`.

| Chain | Transferred residues | Count |
|---|---|---|
| A | 1, 116–127 | 13 |
| B | 1, 115–127 | 14 |
| C | 1–2, 115–127 | 15 |

**Table X.1.2b.** Graft joins closed by rigid translation. The shift in each case equals
the initial C–N distance less the canonical 1.329 Å. Source:
`outputs/step02_join_adjustment_report.tsv`.

| Chain | Join | Initial C–N (Å) | Shift (Å) | Block translated |
|---|---|---|---|---|
| A | 115–116 | 2.120 | 0.791 | A:116–127 |
| B | 1–2 | 2.649 | 1.320 | B:1 |

I checked peptide continuity across all 378 backbone links before and after closure,
accepting a link where the C–N distance fell between 1.15 and 1.70 Å. Two links fell
outside that window before closure and none after. The longest link in the repaired
structure was the graft join at C 114–115, at 1.647 Å.

Bond angles and dihedrals at the graft were not evaluated at this stage. The repaired
structure conditioned the join geometrically only, and the graft geometry was relaxed by
the minimisation and equilibration described in X.2, where it was re-checked.

## Pipeline of Operations

- `step02_terminal_repair_abc.sh` — extracts chains A, B and C from the 2CHT deposition,
  superposes 1DBF globally onto that trimer, transfers the unresolved terminal residues,
  closes permitted graft joins by rigid translation to the canonical peptide C–N distance,
  and reports the superposition, the residue provenance, the peptide continuity before and
  after closure, and the acceptance decision. Reads
  `01_inputs/structures/2cht_raw.pdb` and `01_inputs/structures/1dbf_raw.pdb`. Outputs
  `02_preparation/protein_only/2cht_abc_raw_protein.pdb`,
  `02_preparation/protein_only/1dbf_abc_aligned_to_2cht.pdb`,
  `02_preparation/protein_only/abc_repaired_unadjusted.pdb`,
  `02_preparation/protein_only/abc_repaired.pdb`,
  `00_admin/step02_alignment_report.tsv`,
  `00_admin/step02_residue_source_report.tsv`,
  `00_admin/step02_peptide_continuity_report.tsv`,
  `00_admin/step02_join_adjustment_report.tsv`,
  `00_admin/step02_repair_decision.txt`.
- `step02_verify.py` — confirms by SHA-256 that both inputs are the files produced at
  X.1.1, re-derives the superposition deviation from the committed structures, compares the
  grafted residue set in the repaired structure against both the script definition and the
  published ranges, recomputes every peptide C–N distance under the script's own acceptance
  window, and confirms that each recorded shift equals the initial distance less the
  canonical value. Outputs `checks/step02_verify.out`.

`[FIGURE X.1.2 — terminal repair schematic. Ribbon of the A/B/C trimer with the
transferred residues from Table X.1.2a highlighted per chain, and the two closed joins from
Table X.1.2b marked. An inset should show the C 114–115 join at 1.647 Å against the two
closed joins at 1.329 Å. Render from outputs/abc_repaired.pdb.]`

---

# Sourcing gaps

1. **`lategan_computational_2025.pdf`** — house-style exemplar, not held.

No new citation gaps arise at this step. Ladner *et al.* (2000) is already recorded as
outstanding at X.1.1; the DOI `10.1107/S0907444900004625` is held, the PDF is not.

# Inconsistencies found

1. **The acceptance window is wide, and one graft join exploits it.** A link is accepted
   between 1.15 and 1.70 Å against a canonical 1.329 Å. The graft join at C 114–115 sits at
   1.647 Å, 0.318 Å longer than canonical, and was therefore never flagged and never
   closed, while the equivalent joins on chains A and B were closed to 1.329 Å exactly. The
   three chains consequently enter minimisation with graft joins of unequal quality. The
   structure passes the gate as written; the observation is that the gate does not
   discriminate at this level.
2. **The A-chain closure translated twelve residues.** Closing the A 115–116 join moved the
   whole block A:116–127 rigidly by 0.791 Å relative to the crystallographic frame. The
   block is terminal and remote from the active site, and the displacement is recorded only
   in `step02_join_adjustment_report.tsv`.
3. **Graft bond angles and dihedrals remain unchecked.** The script defers them to
   post-minimisation and records that the affected joins must be re-checked there. That
   re-check has not been performed and is outstanding at X.2.

# Verification log

Performed by `checks/step02_verify.py` against the committed copies in
`02_terminal_repair/outputs/`. Nineteen checks, nineteen passed, none failed. Full output
in `checks/step02_verify.out`; the execution environment is recorded in `manifest.txt`.

- Continuity: `2cht_raw.pdb` and `1dbf_raw.pdb` confirmed byte-identical to the X.1.1
  outputs by SHA-256 against `sha256_step01_inputs.txt`.
- Chain extraction: `2cht_abc_raw_protein.pdb` holds chains A, B and C only.
- Superposition re-derived from `2cht_abc_raw_protein.pdb` and
  `1dbf_abc_aligned_to_2cht.pdb`: 1396 common backbone atoms, deviation 0.8037 Å, matching
  the report to four decimal places.
- Grafted residue set in `abc_repaired.pdb` matches the script patch definition for all
  three chains, and matches the ranges published by Agbaglo *et al.* (2024) for all three
  chains.
- Peptide C–N distances recomputed under the script's own window (1.15–1.70 Å): 378 links
  at both stages, 2 outside the window before closure, 0 after. Longest link after closure
  C 114–115 at 1.647 Å.
- Both recorded shifts equal the initial C–N distance less 1.329 Å to within 0.001 Å, and
  both joins are members of the script's permitted set.
- No barrier is quoted in this subsection, so no barrier re-derivation applies.

# References

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2024) 'The influence of model
building schemes and molecular dynamics sampling on QM-cluster models: the chorismate
mutase case study', *Physical Chemistry Chemical Physics*, 26, pp. 12467–12482.
doi:10.1039/d3cp06100k.

Ladner, J.E., Reddy, P., Davis, A., Tordova, M., Howard, A.J. and Gilliland, G.L. (2000)
'The 1.30 Å resolution structure of the *Bacillus subtilis* chorismate mutase catalytic
homotrimer', *Acta Crystallographica Section D*, 56, pp. 673–683.
doi:10.1107/S0907444900004625.

# X.1.3 Cleanup and structural validation

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Source files for every statement below are in
> `consolidated_works/03_cleanup_validate/`.

---

I stripped the crystallographic hydrogens from the repaired trimer, resolved its
alternate-location groups to a single conformer, and validated the resulting heavy-atom
scaffold, using `step03_cleanup_validate.sh`. The step alters no coordinate; it removes
and selects atoms only.

I detected hydrogens from the PDB element column rather than from atom names, since
name-based detection matches heavy atoms whose names begin with H. The element column was
populated for all 3638 input atoms, so the name-based fallback available in
`step03_cleanup_validate.sh` never engaged. Removal by element column took 601 hydrogens.

I resolved each alternate-location group to the conformer of highest occupancy, breaking
ties in favour of altloc A and then of the unflagged record. Seven groups required
resolution, distributed over two residues. In chain A, Thr122 carried conformers at
occupancies 0.47 and 0.53 across CB, OG1 and CG2, and conformer B was retained. In chain
C, Leu127 carried conformers at 0.64 and 0.36 across CB, CG, CD1 and CD2, and conformer A
was retained. Both residues lie in the C-terminal segments grafted from 1DBF in X.1.2, and
neither belongs to the contact set defined in X.1.4.

**Table X.1.3a.** Alternate-location groups resolved, with the occupancy margin between
conformers and the resulting displacement of the retained side chain. Source:
`outputs/step03_cleanup_report.tsv`, `checks/step03_verify.out`.

| Residue | Atoms | Occupancy A | Occupancy B | Retained | Conformer RMSD (Å) | Largest separation (Å) |
|---|---|---|---|---|---|---|
| A Thr122 | CB, OG1, CG2 | 0.47 | 0.53 | B | 2.123 | 2.673 |
| C Leu127 | CB, CG, CD1, CD2 | 0.64 | 0.36 | A | 1.941 | 3.246 |

The cleanup reduced 3638 input atoms to 3030 heavy atoms, by the removal of 601 hydrogens
and 7 alternate-location duplicates. The validation gates confirmed that no hydrogen, no
alternate-location flag, no duplicate atom name and no non-standard residue name remained,
that three C-terminal OXT atoms were present, and that all 378 peptide links fell within
the acceptance window. All twelve contact-set residues were present in all three chains,
and each chain comprised 127 residues numbered 1 to 127.

**Table X.1.3b.** Cleanup and validation audit. Source:
`outputs/step03_validation_audit.txt`.

| Quantity | Value |
|---|---|
| Input atoms | 3638 |
| Element column populated | 3638 |
| Hydrogens removed | 601 |
| Alternate-location groups resolved | 7 |
| Final heavy atoms | 3030 |
| Residual hydrogens | 0 |
| Residual alternate locations | 0 |
| C-terminal OXT | 3 |
| Non-standard residue names | none |
| Residues with duplicate atom names | 0 |
| Peptide links checked | 378 |
| Peptide links outside 1.15–1.70 Å | 0 |
| Contact-set residues absent | 0 |

## Pipeline of Operations

- `step03_cleanup_validate.sh` — removes crystallographic hydrogens by element column,
  resolves alternate-location groups by occupancy, renumbers and rewrites the structure,
  then re-parses the written file and validates it against six gates. Outputs
  `02_preparation/protein_only/abc_repaired_clean.pdb`,
  `00_admin/step03_cleanup_report.tsv`, `00_admin/step03_validation_audit.txt`.
- `step03_verify.py` — confirms that the input is the committed step 02 output, reproduces
  the atom-count arithmetic from the structures, recomputes peptide continuity, reproduces
  the occupancy rule and reports the occupancy margin and conformer separation for each
  resolved group, and confirms that no retained atom changed coordinates. Outputs
  `checks/step03_verify.out`.

`[FIGURE X.1.3 — the two resolved alternate-location groups. Render A Thr122 and C Leu127
with both conformers overlaid, the retained one solid and the discarded one transparent,
annotated with the occupancies from Table X.1.3a. Source structure
02_preparation/protein_only/abc_repaired.pdb, which still carries both conformers. Caption
should state that both residues lie in grafted C-terminal segments and outside the contact
set.]`

---

# Sourcing gaps

1. **No literature comparator exists for this step.** Neither Claeyssens *et al.* (2011)
   nor Agbaglo *et al.* (2024) report their treatment of alternate locations or of
   crystallographic hydrogens. The write-up should state that no precedent is available
   rather than imply one.

# Inconsistencies found

1. **One alternate-location decision rests on a 0.06 occupancy margin and is
   structurally real.** A Thr122 was resolved in favour of conformer B at occupancy 0.53
   against 0.47. The retained and discarded conformers differ by 2.123 Å RMSD over the
   three affected atoms, with OG1 and CG2 displaced by 2.5 and 2.7 Å. The occupancy rule
   is applied consistently and is defensible, but at this margin the outcome is close to
   arbitrary. Neither the margin nor the conformer separation is recorded by
   `step03_cleanup_validate.sh`; both are produced by `checks/step03_verify.py`. Both
   affected residues lie outside the contact set, in the C-terminal segments that X.2
   subsequently identifies as the only flexible region of the trimer, so the choice does
   not propagate to the reaction.
2. **The alternate locations are inherited from the donor structure.** Both Thr122 and
   Leu127 were grafted from 1DBF rather than taken from 2CHT, per
   `step02_residue_source_report.tsv`. The disorder is a property of the 1.30 Å donor,
   which resolves alternate conformations that the scaffold does not.
3. **The peptide-continuity window is unchanged from X.1.2.** The graft join at C 114–115
   at 1.647 Å passes here for the same reason it passed there, and is not re-flagged.

# Verification log

Performed by `checks/step03_verify.py`. Twenty-eight checks, twenty-eight passed, none
failed. Full output in `checks/step03_verify.out`; environment recorded in `manifest.txt`.

- Input continuity confirmed against the committed `abc_repaired.pdb` from X.1.2.
- Atom-count arithmetic reproduced from the structures rather than read from the audit:
  3638 − 601 − 7 = 3030, matching the 3030 heavy atoms written.
- Output purity confirmed: no hydrogens, no alternate-location flags, no duplicate atom
  names, no non-standard residue names, three OXT.
- Peptide continuity recomputed: 378 links checked, none outside 1.15–1.70 Å. Worst link
  C 114–115 at 1.647 Å.
- Occupancy rule reproduced independently for all seven groups; the retained conformer in
  every group is the higher-occupancy one.
- Coordinate invariance confirmed: maximum displacement of any retained atom, 0.000 Å.
- No barrier is quoted in this subsection, so no barrier re-derivation applies.

# References

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2024) 'The influence of model
building schemes and molecular dynamics sampling on QM-cluster models: the chorismate
mutase case study', *Physical Chemistry Chemical Physics*, 26, pp. 12467–12482.
doi:10.1039/d3cp06100k.

Claeyssens, F., Ranaghan, K.E., Lawan, N., Macrae, S.J., Manby, F.R., Harvey, J.N. and
Mulholland, A.J. (2011) 'Analysis of chorismate mutase catalysis by QM/MM modelling of
enzyme-catalysed and uncatalysed reactions', *Organic & Biomolecular Chemistry*, 9,
pp. 1578–1590. doi:10.1039/c0ob00691b.

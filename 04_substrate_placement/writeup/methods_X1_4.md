# X.1.4 Substrate placement and the active-site contact set

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Source files for every statement below are in
> `consolidated_works/04_substrate_placement/`.

---

I placed one chorismate into each of the three active sites and derived the set of
residues forming those sites from the placed geometry, using
`step04_place_chorismate_abc.sh`.

Each substrate geometry had been aligned to a specific crystallographic analogue site in
X.1.1. I superposed the analogue of that source site onto the analogue of the target site
and applied the resulting transform to the substrate, carrying it into the target site in
the analogue pose. The fits used 17 common analogue atoms and gave RMSDs of 0.1267,
0.0797 and 0.1518 Å for sites A, B and C. I renamed the duplicate atom names inherited
from the reference geometries to unique labels for tleap safety, preserving coordinates,
atom types, charges and bonds; each placed substrate retained 24 atoms, 16 of them heavy,
at net charge −2.

**Table X.1.4a.** Placement transforms. Source: `outputs/step04_placement_transform_report.tsv`.

| Substrate | Source analogue | Target analogue | Common atoms | Fit RMSD (Å) | Centre offset (Å) |
|---|---|---|---|---|---|
| `cha_a` | K210 | A203 | 17 | 0.1267 | 0.161 |
| `cha_b` | J212 | B201 | 17 | 0.0797 | 0.160 |
| `cha_c` | L211 | C202 | 17 | 0.1518 | 0.203 |

I assigned a residue to the contact set if any of its heavy atoms lay within 4.0 Å of any
heavy atom of the placed substrate in at least one of the three sites. Distances were
measured against every protein heavy atom in the trimer rather than against a predefined
residue list, so the membership of the set and the chain origin of each contact were
outcomes of the measurement. The chorismate mutase active site forms at the interface of
two adjacent monomers, and measuring across all chains tests that arrangement rather than
assuming it.

The criterion returned twelve residues: Arg7, Glu78, Arg90, Tyr108 and Leu115 from the
substrate's own chain, and Phe57′, Ala59′, Lys60′, Arg63′, Val73′, Thr74′ and Cys75′ from
the adjacent chain. Twelve residues met the criterion at sites A and B and eleven at site
C, where Arg63′ lay at 6.380 Å; the set is the union over the three sites. The assignment
proved stable with respect to the cutoff, returning the same twelve residues at 3.5 Å and
at 4.0 Å, collapsing to six below 3.5 Å, and admitting a thirteenth only at 4.5 Å. The
nearest excluded residue was Ala9, at 4.449 Å.

**Table X.1.4b.** Closest heavy-atom distance (Å) between each contact-set residue and the
placed substrate, per site. Primed residues are contributed by the adjacent chain. Source:
`outputs/step04_placement_contact_report.txt`, `outputs/step04_catalytic_contact_report.tsv`.

| Residue | Origin | Site A | Site B | Site C |
|---|---|---|---|---|
| Arg7 | same | 2.782 | 2.782 | 2.709 |
| Arg90 | same | 2.845 | 2.865 | 2.997 |
| Glu78 | same | 3.362 | 2.900 | 3.098 |
| Tyr108 | same | 2.956 | 2.586 | 2.785 |
| Leu115 | same | 3.548 | 3.764 | 3.301 |
| Phe57′ | cross | 3.620 | 3.684 | 3.390 |
| Ala59′ | cross | 3.358 | 3.222 | 2.964 |
| Lys60′ | cross | 3.024 | 3.129 | 3.080 |
| Arg63′ | cross | 2.899 | 3.060 | 6.380 |
| Val73′ | cross | 3.594 | 3.525 | 3.492 |
| Thr74′ | cross | 3.448 | 3.266 | 3.612 |
| Cys75′ | cross | 3.384 | 3.231 | 3.188 |

**Table X.1.4c.** Residues excluded by the criterion, by closest approach in any site.
Source: `checks/step04_verify.out`.

| Residue | Closest approach (Å) |
|---|---|
| Ala9 | 4.449 |
| Pro58′ | 5.164 |
| Gly8 | 5.209 |
| Arg116 | 5.501 |
| Pro117 | 5.515 |

Every expected chain-origin assignment was recovered by measurement, across 39 tested
assignments in all three sites, with no mismatch. The closest protein-substrate contact
was 2.782 Å at site A, 2.586 Å at site B and 2.709 Å at site C, against a severe-overlap
threshold of 0.80 Å.

The twelve residues returned by this criterion were the twelve that Agbaglo *et al.* (2024)
obtained for this enzyme by an independent route, as the union of residues appearing in
their crystal-structure quantum-cluster models. Two sets derived by different methods from
different starting structures coinciding exactly is a check on the placement, and I report
it as such rather than as the source of the set. A correction to that work was published
in 2026, restating the reference structure used for the deviations in their frame
selection and revising four of their reported free energies; neither amendment touches the
residue composition used here (Agbaglo *et al.*, 2026).

Arg116 lay outside the contact set at 5.501–6.834 Å. Claeyssens *et al.* (2011) report it
as contributing 1.0 ± 0.9 kcal mol⁻¹ of electrostatic transition-state stabilisation and
describe it as an active-site residue. The two statements are compatible, since a contact
criterion selects residues in van der Waals proximity whereas an electrostatic
contribution does not require proximity. I use the contact set for the geometric purpose
it was constructed for, confirming that each substrate copy occupies the crystallographic
pose, and do not carry it forward as a definition of the quantum region or as the residue
basis for the electrostatic decomposition. Those partitions are defined on their own
criteria in X.3 and in the corresponding results section.

`[DECISION REQUIRED — the residue basis for the per-residue electrostatic decomposition is
left open here. Claeyssens et al. (2011, Table 3) decompose over eight residues, including
Arg116 and excluding Lys60′, Val73′, Thr74′, Ala59′ and Leu115, all of which are contacts
in this structure. Their selection reflects sixteen reaction paths; the ensemble here is
larger, so the selection step could be removed and contributions computed over every
residue within a stated radius. ensemble_barriers.tsv currently carries only
whole-environment interaction terms with no per-residue breakdown, so this requires work
that does not yet exist and may fall outside the Phase 1 scope.]`

## Pipeline of Operations

- `step04_place_chorismate_abc.sh` — superposes each source analogue onto its target
  analogue, applies the transform to the substrate, renames duplicate atom names, writes
  the combined structure, and measures protein-substrate contacts both as an unbiased
  4.0 Å scan over all residues and as a chain-origin test over a named residue set.
  Outputs `02_preparation/ligand_placement/abc_with_chorismate_unprotonated.pdb`,
  `02_preparation/ligand_placement/cha_a_placed.mol2`,
  `02_preparation/ligand_placement/cha_b_placed.mol2`,
  `02_preparation/ligand_placement/cha_c_placed.mol2`,
  `00_admin/step04_placement_transform_report.tsv`,
  `00_admin/step04_placement_contact_report.txt`,
  `00_admin/step04_placement_audit.txt`,
  `00_admin/step04_catalytic_contact_report.tsv`.
- `step04_verify.py` — confirms input continuity against X.1.1 and X.1.3, reproduces each
  analogue-to-analogue transform from the raw deposition, confirms substrate integrity and
  net charge after placement, reproduces the unbiased contact scan entry for entry,
  compares the resulting set against the published set, and reports the cutoff sensitivity
  and exclusion margin that the placement script does not produce. Outputs
  `checks/step04_verify.out`.

`[FIGURE X.1.4a — one active site, rendered from abc_with_chorismate_unprotonated.pdb,
showing the placed substrate with the twelve contact-set residues, coloured by chain of
origin to make the inter-subunit arrangement visible, and Arg116 shown in outline as the
nearest excluded residue.]`

`[FIGURE X.1.4b — cutoff-sensitivity plot. Number of residues in the contact set against
cutoff over 3.0–6.0 Å, from checks/step04_verify.out, marking the 3.5–4.0 Å plateau at
n = 12 and the 4.0 Å working value. Caption should state that the nearest excluded residue
is Ala9 at 4.449 Å.]`

---

# Sourcing gaps

1. **Ranaghan *et al.* (2003, 2004)** — the substitution procedure cited by
   Claeyssens *et al.* (2011) as their references 9 and 10. Neither PDF is held and no DOI
   is available for either.
2. **Bartlett and Johnson (1985)** — the transition-state analogue, required if the
   analogue is cited by name. No DOI held.

# Inconsistencies found

1. **The archived copy of the placement script is materially stale.** The copy in
   `00_admin/` uses a residue list containing Arg116 and omitting Val73′. The committed
   `step04_catalytic_contact_report.tsv` contains Val73′ and labels Arg116 as second shell,
   so the report was produced by the working version. This is the only step so far where
   the archived and working copies differ in substance rather than in environment
   settings.
2. **The script comment misattributes the published residue set.** It credits the
   supporting information of Agbaglo *et al.* (2024); the set is in the main text, p. 8,
   in the description of frame-selection scheme S3. The comment also repeats the
   pre-correction wording of that scheme, which the 2026 correction revises.
3. **Cutoff sensitivity is not recorded by the placement script.** The 4.0 Å cutoff is
   applied without any report of how the set would change at other values. The ladder and
   the exclusion margin are produced by `checks/step04_verify.py`.

# Verification log

Performed by `checks/step04_verify.py`. Thirty-two checks, thirty-two passed, none failed.

- Input continuity confirmed against X.1.3 and X.1.1 by SHA-256.
- Source analogue sites hard-coded in the placement script match the registration measured
  in X.1.1 for all three substrates.
- All three analogue-to-analogue transforms reproduced from the raw deposition: 17 common
  atoms and RMSD 0.1267, 0.0797, 0.1518 Å, matching to four decimal places.
- Substrate integrity after placement: 24 atoms, 16 heavy, net charge −2.0000, for all
  three.
- Unbiased 4.0 Å contact scan reproduced entry for entry, maximum distance discrepancy
  0.0009 Å.
- Union of the three site sets equals the twelve residues published by Agbaglo *et al.*
- Cutoff ladder: 6, 12, 12, 13, 13, 21 residues at 3.0, 3.5, 4.0, 4.5, 5.0 and 6.0 Å.
- Chain-origin assignments: 39 rows tested, 0 mismatches; Arg116 labelled second shell.
- Minimum protein-substrate contacts 2.782, 2.586 and 2.709 Å, all above the 0.80 Å
  overlap threshold.
- No barrier is quoted in this subsection, so no barrier re-derivation applies.

# References

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2024) 'The influence of model
building schemes and molecular dynamics sampling on QM-cluster models: the chorismate
mutase case study', *Physical Chemistry Chemical Physics*, 26, pp. 12467–12482.
doi:10.1039/d3cp06100k.

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2026) 'Correction: The
influence of model building schemes and molecular dynamics sampling on QM-cluster models:
the chorismate mutase case study', *Physical Chemistry Chemical Physics*, 28,
pp. 8283–8284. doi:10.1039/d6cp90040b.

Claeyssens, F., Ranaghan, K.E., Lawan, N., Macrae, S.J., Manby, F.R., Harvey, J.N. and
Mulholland, A.J. (2011) 'Analysis of chorismate mutase catalysis by QM/MM modelling of
enzyme-catalysed and uncatalysed reactions', *Organic & Biomolecular Chemistry*, 9,
pp. 1578–1590. doi:10.1039/c0ob00691b.

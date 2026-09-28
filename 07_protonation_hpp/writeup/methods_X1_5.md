# X.1.5 Acceptance freeze and protonation-state assignment

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Assembled from the freeze fragment in `06_acceptance_freeze/writeup/` and the
> protonation half drafted here. **The H++ run parameters are left blank for completion
> from the lab workbook**; see Sourcing gaps.

---

I froze the accepted pre-protonation complex using `step06_freeze_preprotonation.sh`. The
substrate geometries inherited repeated generic atom names from their source, which tleap
cannot resolve, so I renamed every substrate atom uniquely before freezing. The renaming
preserved coordinates, atom types, partial charges and bond connectivity; each substrate
retained 24 atoms in 24 bonds at a net charge of exactly −2, the value assigned to the
quantum region in X.3. The renaming tables were recorded per substrate in
`cha_a_renaming.tsv`, `cha_b_renaming.tsv` and `cha_c_renaming.tsv`. The combined complex
comprised 3030 protein atoms and 72 substrate heteroatoms across three substrate residues,
with no residue carrying duplicate atom names. I checksummed the four accepted files in
`sha256_step06_accepted_preprotonation.txt`.

**Table X.1.5a.** Substrate audit after renaming. Source:
`06_acceptance_freeze/outputs/step06_ligand_mol2_audit.tsv`.

| Substrate | Atoms | Heavy atoms | Bonds | Net charge |
|---|---|---|---|---|
| `cha_a` | 24 | 16 | 24 | −2.000000 |
| `cha_b` | 24 | 16 | 24 | −2.000000 |
| `cha_c` | 24 | 16 | 24 | −2.000000 |

I assigned protonation states with the H++ server. `step07a_prepare_hpp_input.sh` wrote the
protein-only, TER-delimited structure for submission, comprising three chains of 1010 atoms
over residues 1 to 127, carrying no heteroatoms and no hydrogens. The server was run at
pH 7.0, salinity \_\_\_\_\_\_ M, internal dielectric \_\_\_\_\_\_ and external dielectric
\_\_\_\_\_\_.

The server returned 381 residues carrying 3153 added hydrogens, with all 3030 template
heavy atoms present and none added. Nine histidines were assigned: His36 and His54 as HIP
and His106 as HIE, identically in each of the three chains. I used only the protonation and
orientation decisions from this output, discarding the H++ topology and coordinates and
applying the assignments to the frozen template.

`step07b_audit_repair_hpp.sh` classified every heavy-atom displacement between the template
and the accepted structure and applied four acceptance thresholds. Backbone drift was
0.0000 Å over 1527 atoms against a threshold of 0.05 Å, non-flippable side-chain drift
0.0000 Å over 2676 atoms against 0.20 Å, and active-site backbone drift 0.0000 Å over 144
atoms against 0.05 Å. No atom of the contact set moved at all. Twenty-two heavy atoms
moved beyond 1.0 Å, to a maximum of 2.2452 Å, and all twenty-two lay in glutamine and
histidine side chains outside the contact set, which is the expected signature of
terminal-amide and imidazole flips.

**Table X.1.5b.** Heavy-atom displacement between the frozen template and the accepted
protonated structure, by category. Source:
`outputs/step07b_audit_summary.tsv`.

| Category | Atoms | Maximum (Å) | RMS (Å) | Beyond 1.0 Å |
|---|---|---|---|---|
| Backbone | 1527 | 0.0000 | 0.0000 | 0 |
| Non-flippable residues | 2676 | 0.0000 | 0.0000 | 0 |
| Contact-set residues | 321 | 0.0000 | 0.0000 | 0 |
| Contact-set backbone | 144 | 0.0000 | 0.0000 | 0 |
| Flippable residues | 354 | 2.2452 | 0.4959 | 22 |
| All heavy atoms | 3030 | 2.2452 | 0.1695 | 22 |

No residue name fell outside ff14SB, no heavy atom was missing or extra, all 378 peptide
links were normal, no residue carried duplicate atom names, and no disulfide was detected.
The structure was accepted for parameterisation.

## Pipeline of Operations

- `step06_freeze_preprotonation.sh` — renames substrate atoms uniquely, rebuilds the
  combined complex, audits substrate charge and name uniqueness, and freezes the accepted
  files with checksums. Outputs
  `02_preparation/accepted_preprotonation/abc_cha_preprotonation.pdb`,
  `02_preparation/accepted_preprotonation/cha_{a,b,c}.mol2`,
  `02_preparation/accepted_preprotonation/cha_{a,b,c}_renaming.tsv`,
  `00_admin/step06_acceptance_decision.txt`, `00_admin/step06_combined_pdb_audit.txt`,
  `00_admin/step06_ligand_mol2_audit.tsv`,
  `00_admin/sha256_step06_accepted_preprotonation.txt`.
- `step07a_prepare_hpp_input.sh` — writes the protein-only, TER-delimited structure for
  submission to H++. Outputs
  `03_amber/protonation_hpp/hpp_input_protein_only_with_TER.pdb`,
  `00_admin/step07a_hpp_input_report.tsv`.
- `step07b_audit_repair_hpp.sh` — audits the returned structure against the frozen
  template, classifies heavy-atom displacement by backbone, flippability and contact-set
  membership, derives the histidine assignment from ring-hydrogen presence, and records an
  explicit accept decision against four thresholds. Outputs
  `03_amber/protonation_hpp/abc_protonated_hpp_accepted.pdb`,
  `00_admin/step07b_audit_summary.tsv`, `00_admin/step07b_his_assignment.tsv`,
  `00_admin/step07b_drift_by_atom.tsv`, `00_admin/step07b_drift_by_residue.tsv`,
  `00_admin/step07b_acceptance_decision.txt`.
- `step06_verify.py`, `step07_verify.py` — reproduce the renaming, the composition audits,
  every displacement category, the acceptance thresholds and the histidine assignment, and
  confirm that neither the freeze nor the protonation moved any template coordinate.
  Outputs `checks/step06_verify.out`, `checks/step07_verify.out`.

`[FIGURE X.1.5 — flowchart from the cleaned trimer and placed substrates through the
acceptance freeze, submission to H++, and the audited acceptance, terminating at
abc_protonated_hpp_accepted.pdb as the input to X.1.6. Mark the external boundary at the
server submission, since that step is not reproducible from the repository. Annotate the
acceptance gate with its four thresholds and the achieved values.]`

---

# Sourcing gaps

1. **H++ run parameters.** pH 7.0 is agreed across all copies of the script comments.
   Salinity and the internal and external dielectric constants appear only in those
   comments, and the copies contradict one another: the archived copy of
   `step07a_prepare_hpp_input.sh` records internal 80 and external 10, while the working
   copy and `step07b_audit_repair_hpp.sh` both record internal 10 and external 80.
   `hpp_output_from_server.pdb` carries no parameter record, only a URL remark. These
   values are a record of a past interactive submission and cannot be recovered by
   re-running anything. **Blanks left in the text for completion from the lab workbook.**
2. **H++ method paper.** No DOI held. Search terms: *Anandakrishnan, Aguilar, Onufriev —
   H++ 3.0 automating pK prediction and the preparation of biomolecular structures —
   Nucleic Acids Research 2012*.
3. **No literature comparator for protonation assignment.** Claeyssens *et al.* (2011)
   state only that the structure was hydrogenated and solvated. Agbaglo *et al.* (2024)
   report no protonation protocol. The write-up should state that no precedent is
   available rather than imply one.

# Inconsistencies found

1. **The two committed copies of `step07a_prepare_hpp_input.sh` disagree on the dielectric
   constants**, as above. Two records favour internal 10 and external 80, which is the
   conventional assignment for a solute in water, but no committed output establishes it.
2. **An undeclared output is present.** `step07b_heavy_atom_drift.tsv` exists in the
   results repository but is written by no script in the pipeline. Its summary figures
   agree with the recomputed values, so it is consistent, but its provenance is unknown
   and it duplicates `step07b_drift_by_residue.tsv` in a different format.
3. **The as-run archive convention changes partway through the pipeline.**
   `step07b_audit_repair_hpp.sh` is archived under `phase1_system_dev/` in the results
   repository rather than in `00_admin/` with its siblings. It is identical to the working
   copy.
4. **`step06_acceptance_decision.txt` names a superseded residue list.** Generated as a
   literal block inside `step06_freeze_preprotonation.sh`, it enumerates the adjacent-chain
   contributors as Arg63′, Lys60′, Thr74′, Cys75′, Phe57′ and Ala59′, omitting Val73′,
   which X.1.4 establishes as a contact at 3.49–3.59 Å in all three sites. The frozen
   structure is unaffected.

# Verification log

Performed by `checks/step06_verify.py` and `checks/step07_verify.py`. Forty-three and
forty-two checks respectively, all passed.

- Renaming confirmed as a pure relabelling: coordinates identical to 0.00e+00 Å, atom types
  and charges unchanged, bond blocks identical, renaming tables consistent with input and
  output. Net charge exactly −2.000000 in all three substrates.
- Combined complex reproduced: 3030 + 72 = 3102 records, 3 substrate residues, none with
  duplicate names. Protein coordinates unchanged from X.1.3 at 0.00e+00 Å.
- Submission structure reproduced: three chains of 1010 atoms, residues 1–127, three TER
  records, no heteroatoms, no hydrogens.
- Conservation confirmed: all 3030 template heavy atoms present in the return, none added,
  3153 hydrogens.
- All nine displacement categories reproduced independently and matched the recorded
  summary to four decimal places.
- All four acceptance thresholds reproduced. Every displacement beyond 1.0 Å resolves to a
  glutamine or histidine side chain.
- Histidine assignment re-derived from ring-hydrogen presence in the accepted structure and
  matched the recorded assignment for all nine, with identical assignment across the three
  chains.

# References

Claeyssens, F., Ranaghan, K.E., Lawan, N., Macrae, S.J., Manby, F.R., Harvey, J.N. and
Mulholland, A.J. (2011) 'Analysis of chorismate mutase catalysis by QM/MM modelling of
enzyme-catalysed and uncatalysed reactions', *Organic & Biomolecular Chemistry*, 9,
pp. 1578–1590. doi:10.1039/c0ob00691b.

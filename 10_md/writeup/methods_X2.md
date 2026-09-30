# X.2 Conformational sampling and reaction-path frame selection

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Source files are in `consolidated_works/10_md/`, `11_trajectory_analysis/` and
> `12_frame_selection/`. Covers pipeline steps 10, 11 and 12.

---

## X.2.1 Molecular dynamics

I minimised, equilibrated and then ran unrestrained dynamics on the solvated complex of
X.1.7, following the protocol of Agbaglo *et al.* (2024) for the same enzyme. The
minimisation restrained the solute heavy atoms harmonically at 200 kcal mol⁻¹ Å⁻² over
10000 cycles, leaving water and hydrogens free. I then relaxed the restraint over five
consecutive 20 ps stages of constant-pressure dynamics at 200, 100, 50, 10 and
0 kcal mol⁻¹ Å⁻², using Langevin dynamics at 300 K and 1 atm with a collision frequency
of 2.0 ps⁻¹, the SHAKE algorithm on all bonds involving hydrogen, and a 9 Å cut-off for
non-bonded interactions.

Production ran unrestrained under the same thermostat and barostat for 20 ns with a 2 fs
timestep, writing coordinates every 1 ps to give 20000 frames. I extended production by a
further 40 ns under identical conditions, giving 60 ns in total.

**Table X.2.1a.** Simulation protocol, with the corresponding specification of Agbaglo
*et al.* (2024). Source: `outputs/10a_min/min10a.in`, `outputs/10b_equil/npt{1..5}.in`,
`outputs/10c_production/prod.in`, `outputs/10d_production_extend/prod_ext.in`.

| Stage | This work | Agbaglo *et al.* (2024) |
|---|---|---|
| Minimisation restraint | 200 kcal mol⁻¹ Å⁻², solute heavy atoms | 200 kcal mol⁻¹ Å⁻², protein heavy atoms |
| Equilibration | 5 × 20 ps, 200 → 100 → 50 → 10 → 0 | five 20 ps simulations, iteratively relaxed |
| Ensemble | NPT, 300 K, 1 atm, Langevin | NPT, 300 K, 1 atm, Langevin |
| Constraints | SHAKE on bonds to hydrogen | SHAKE on bonds to hydrogen |
| Non-bonded cut-off | 9 Å | 9 Å |
| Production | 20 ns, 1 ps per frame, 20000 frames | 20 ns, 1 ps per frame, 20000 frames |
| Extension | **+40 ns, 60 ns total** | none |

The extension departed from the published protocol, which ran production for 20 ns. Two
thirds of the reaction-path frames selected in X.2.3 came from beyond 20 ns, so the
sampling underlying the reaction-path results covered a range the published protocol did
not.

## X.2.2 Trajectory quality

I assessed the trajectory before selecting any frame for reaction-path calculation,
measuring backbone deviation, per-residue flexibility, active-site contact persistence and
the reaction coordinate itself. The analysis used a standalone reader rather than cpptraj,
which was unavailable across the cluster.

Backbone root-mean-square deviation from the first production frame averaged 3.69 Å over
the full 60 ns with a maximum of 4.53 Å. Successive 10 ns blocks gave means of 2.65, 3.67,
3.89, 4.03, 3.83 and 4.03 Å, so the deviation rose over the first 20 ns and then held.
Agbaglo *et al.* report a mean backbone deviation of 2.66 Å with a standard deviation of
0.76 Å for the same enzyme, measured against their initial equilibrated structure.

Per-residue fluctuation over 381 residues averaged 0.69 Å across the lower four-fifths of
the distribution. The flexible residues were confined to the chain termini, reaching
13.9 Å at B253. The thirteen residues forming the substrate contact set of X.1.4 all fell
below 1.0 Å, with a mean of 0.70 Å and a maximum of 0.97 Å at Arg63 of chain A.

The reaction coordinate, defined as the difference between the breaking C4–O3 and forming
C6–C1 distances, distinguished the three substrate copies sharply. The breaking bond
remained within 1.324–1.581 Å across all copies and all frames, as expected for a bond
held by a classical force field.

**Table X.2.2a.** Substrate geometry by binding site, means over 20000 production frames,
and the proportion of frames within the near-attack distance criterion over the full
60 ns. Source: `outputs/rxn_coord_per_frame.dat`,
`outputs/rxn_coord_per_frame_full.dat`.

| Site | Breaking C4–O3 (Å) | Forming C6–C1 (Å) | *r* (Å) | Frames with forming ≤ 3.7 Å |
|---|---|---|---|---|
| A | 1.45 | 4.15 | −2.70 | 1.8 % |
| **B** | 1.46 | **3.29** | **−1.83** | **99.4 %** |
| C | 1.45 | 3.83 | −2.38 | 14.4 % |

## X.2.3 Frame selection

I selected reaction-path starting structures on two criteria applied in sequence. The
first excluded any binding site whose substrate was only transiently in a reactive
arrangement, requiring the site to satisfy both an intact Arg90–substrate contact below
3.2 Å and a forming distance below 3.7 Å in at least half of the frames sampled. Of the
three sites, one met the threshold and two did not: site A was competent in 2 % of frames,
site B in 95 %, and site C in 22 %. All subsequent reaction-path work therefore used
site B alone.

The second criterion was the near-attack conformer definition of Hur and Bruice (2003),
applied in full. Their definition combines a forming-bond distance below 3.7 Å with two
approach-angle conditions, and I tested all three rather than the distance alone. Of the
twelve frames selected on distance, nine satisfied the complete definition. The three that
did not failed on the second approach angle by 21.4°, 24.7° and 23.0° against a 20°
tolerance, where the unreactive conformers Hur and Bruice identify in water deviate by
approximately 110°. I retained these three, labelled them separately, and selected three
additional frames satisfying all three conditions.

Applying a near-attack criterion as a selection filter carries no commitment to the
near-attack hypothesis as an explanation of catalysis. Ranaghan *et al.* (2004) set out
the case against that hypothesis, noting that the proposed dominance of near-attack
conformer formation has been criticised as unreliable and that free-energy perturbation
work finds no contribution from near-attack conformer formation to catalysis; the same
body of work attributes catalysis in this enzyme to transition-state stabilisation.
The filter here identifies reactant-state structures competent to react, and nothing more.

`[DECISION REQUIRED — whether the near-attack frames are reported separately in the
Results chapter or excluded. Retaining them permits a test: if near-attack geometry
governs reactivity, their barriers should sit systematically higher than those of the
complete-definition frames. Excluding them removes three of the original set and loses
that test.]`

**Table X.2.3a.** Selection criteria and outcomes. Source:
`outputs/step12a_selection.log`, `outputs/12b_nac_select.log`,
`outputs/selection_manifest.tsv`.

| Criterion | Threshold | Outcome |
|---|---|---|
| Site competence | ≥ 50 % of frames | A 2 %, B 95 %, C 22 % — site B only |
| Forming distance | ≤ 3.7 Å | applied per frame |
| Approach angles | within 20° | 9 of 12 initial frames satisfied both |
| Final manifest | — | 45 frames, all site B, 820–59999 ps |

## Pipeline of Operations

- `step10a_min_submit.sh`, `step10b_equil_submit.sh`, `step10c_production_submit.sh`,
  `step10d_production_extend_submit.sh` — minimisation, the five-stage equilibration
  ladder, production and its extension. Outputs `04_amber_md/10a_min/10a_min.rst7`,
  `04_amber_md/10b_equil/npt{1..5}.rst7`, `04_amber_md/10c_production/prod.nc`,
  `04_amber_md/10d_production_extend/prod_ext.nc`, with the corresponding `.out` and
  `.mdinfo` records.
- `step10a_audit.sh`, `step10b_audit.sh`, `step10c_audit.sh` — per-stage audits of
  completion, temperature, pressure and density.
- `step11a_rmsd.sh`, `step11a_rmsd_full.sh` — backbone deviation over production and over
  the full trajectory. Outputs `04_amber_md/11_analysis/rmsd_vs_time.dat`,
  `rmsd_vs_time_ext.dat`, `rmsd_vs_time_full.dat`.
- `step11b_rmsf.sh` — per-residue fluctuation. Outputs
  `04_amber_md/11_analysis/rmsf_per_residue.dat`.
- `step11c_rxn_coord.sh`, `step11c_rxn_coord_full.sh` — reaction coordinate per frame for
  all three substrate copies. Outputs `04_amber_md/11_analysis/rxn_coord_per_frame.dat`,
  `rxn_coord_per_frame_full.dat`.
- `step11d_contacts.sh` — persistence of the catalytic contacts of X.1.4.
- `step11e_autocorr.py`, `theta2_autocorr.py` — autocorrelation of the reaction coordinate
  and of the second approach angle, establishing the interval between independent frames.
- `step11f_nac_check.py` — the complete three-condition near-attack test of Hur and Bruice
  (2003). Outputs `05_qmmm/11f_nac_check.log`.
- `step12a_select.sh` — site competence filter and distance-clause selection. Outputs
  `05_qmmm/12_frame_selection/selection_manifest.tsv`,
  `05_qmmm/12_frame_selection/step12a_selection.log`.
- `step12b_nac_select.py`, `step12b_extract.sh` — supplementary selection against the
  complete definition, and extraction of the selected frames. Outputs
  `05_qmmm/selection_manifest.tsv`, `05_qmmm/12_frame_selection/frames/`.
- `step10_verify.py`, `step11_verify.py`, `step12_verify.py` — verification scripts.
  Outputs the corresponding `.out` files under each step's `checks/`.

`[FIGURE X.2.1 — sampling and selection flowchart. Solvated topology from X.1.7 into
minimisation, the five-stage ladder, production and extension, then the analysis branch
(deviation, fluctuation, reaction coordinate, near-attack test) feeding the two selection
filters and the final manifest. Annotate the extension arrow with 20 ns → 60 ns and the
site filter with 2 % / 95 % / 22 %.]`

`[FIGURE X.2.2 — trajectory quality, three panels. (a) backbone deviation against time
over 60 ns from rmsd_vs_time_full.dat, with 10 ns block means overlaid and a horizontal
line at the 2.66 Å of Agbaglo et al.; (b) per-residue fluctuation from
rmsf_per_residue.dat, contact-set residues marked; (c) forming C6–C1 distance
distributions for the three sites from rxn_coord_per_frame_full.dat, with the 3.7 Å
criterion marked. Caption should state that the flexible residues are confined to the
chain termini and that the three sites differ by a factor of 55 in near-attack
occupancy.]`

---

# Sourcing gaps

None new.

# Inconsistencies found

1. **The trajectory analysis document states a selection criterion that was not applied.**
   `step11_trajectory_analysis.md` closes by specifying selection that "draws from all
   three sites, weighted toward CHA#2". Every one of the 45 selected frames is site B. The
   executed criterion is the site competence threshold, which excludes sites A and C by
   construction, and the document requires correction.
2. **The fluctuation statistics in that document mix two populations.** It quotes a mean
   of 0.69 Å beside a median of 0.74 Å; 0.686 Å is the mean of the lower four-fifths of
   the distribution and 0.737 Å is the median of all 381 residues. The lower-four-fifths
   median is 0.685 Å. Both figures are individually correct and describe different
   populations.
3. **The claim that active-site residues are absent from the flexible set is imprecise.**
   On a strict eightieth-percentile threshold of 0.91 Å, Arg63 and Arg116 of chain A sit
   marginally above it at 0.97 and 0.94 Å. The substantive claim holds, since all thirteen
   contact-set residues fall below 1.0 Å while the termini reach 13.9 Å.
4. **The document covers 20 ns of a 60 ns trajectory.** The full-range statistics reported
   in X.2.2 were computed from the committed extended data files for this write-up.
5. **No as-run script archive exists** for steps 10, 11 or 12, as for X.1.6 and X.1.7.

# Verification log

Performed by `checks/step10_verify.py`, `step11_verify.py` and `step12_verify.py`.
Twenty-nine, twenty-five and eighteen checks respectively; two of the step 11 checks fail,
and those two are inconsistencies 2 and 3 above, retained as failures until the source
document is corrected.

- Protocol read from the committed AMBER inputs and confirmed against Agbaglo *et al.*
  element by element, including the five 20 ps stages and the monotone restraint ladder.
- Peptide continuity: all 378 links within the acceptance window after both minimisation
  and equilibration. The grafted join at C 114–115 relaxed from 1.647 Å to 1.362 Å and is
  no longer the most deviant link in the structure.
- The seven close contacts reported at build: the substrate contact relaxed from 1.375 Å
  to 3.080 Å; the three inter-residue contacts at the graft relaxed to 3.3–4.2 Å; the
  three remaining pairs are geminal hydrogens on single methyl carbons, whose separation
  of 1.78 Å is fixed by covalent geometry.
- All deviation, fluctuation and reaction-coordinate figures reproduced from the committed
  data files.
- Near-attack labelling confirmed traceable: the three frames labelled separately are
  exactly those the complete test rejected.
- Site competence reproduced from the selection log; exactly one site met the threshold.
- All 45 manifest entries confirmed to originate from a single substrate copy.

# References

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2024) 'The influence of model
building schemes and molecular dynamics sampling on QM-cluster models: the chorismate
mutase case study', *Physical Chemistry Chemical Physics*, 26, pp. 12467–12482.
doi:10.1039/d3cp06100k.

Hur, S. and Bruice, T.C. (2003) 'The near attack conformation approach to the study of the
chorismate to prephenate reaction', *Proceedings of the National Academy of Sciences*,
100, pp. 12015–12020. doi:10.1073/pnas.1534873100.

Ranaghan, K.E., Ridder, L., Szefczyk, B., Sokalski, W.A., Hermann, J.C. and Mulholland,
A.J. (2004) 'Transition state stabilization and substrate strain in enzyme catalysis: ab
initio QM/MM modelling of the chorismate mutase reaction', *Organic & Biomolecular
Chemistry*, 2, pp. 968–980. doi:10.1039/b313759g.

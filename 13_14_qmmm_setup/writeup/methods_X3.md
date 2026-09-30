# X.3 Reaction-path calculation and transition-state characterisation

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Source files are in `consolidated_works/13_14_qmmm_setup/`,
> `15_17_path_endpoints_scan/` and `18_nebts_ts_irc/`. Covers pipeline steps 13, 14,
> 15, 16, 17 and 18. **No barrier is reported in this subsection; barriers belong to the
> Results chapter.**

---

## X.3.1 Quantum region and level of theory

I partitioned the system with the substrate as the quantum region and the remainder of
the protein and solvent as the classical region, coupled by electrostatic embedding. The
quantum region comprised the 24 atoms of one chorismate molecule. The substrate is not
covalently bonded to the protein, so the partition cuts no bond and requires no link
atoms.

The partition matches that of Claeyssens *et al.* (2011), who treat the
chorismate, transition state and prephenate as the quantum region at 24 atoms and
likewise require no link atoms.

I converted the AMBER topology of X.1.7 into the force-field format required by the
quantum chemistry package, then confirmed the bridge on a single frame with a low-cost
composite Hartree–Fock method before running the production level of theory on the same
geometry. Both single points terminated normally and converged in 14 cycles.

The production level of theory was the B3LYP hybrid density functional with the D3
dispersion correction and Becke–Johnson damping, the def2-SVP basis set, the def2/J
auxiliary basis, and the resolution-of-identity chain-of-spheres approximation to the
Coulomb and exchange integrals. Claeyssens *et al.* use B3LYP with the 6-31G(d) basis set
and no dispersion correction.

`[DECISION REQUIRED — whether to quantify the basis-set difference. The functional is
common to both, but def2-SVP and 6-31G(d) are different basis sets and the dispersion
correction is additional here. A single point at B3LYP/6-31G(d) on an existing optimised
geometry would give the offset directly. The alternative is to rest the comparison on
differential stabilisation, where a systematic basis-set offset largely cancels, and to
state the absolute offset as unquantified.]`

**Table X.3.1a.** Quantum-chemical treatment, with the corresponding specification of
Claeyssens *et al.* (2011). Source: `outputs/dft_sp.inp`, `outputs/smoke.inp`.

| | This work | Claeyssens *et al.* (2011) |
|---|---|---|
| Quantum region | chorismate, 24 atoms | chorismate/TS/prephenate, 24 atoms |
| Link atoms | none required | none required |
| Functional | B3LYP | B3LYP |
| Dispersion | D3 with Becke–Johnson damping | none |
| Basis set | def2-SVP | 6-31G(d) |
| Integral approximation | RIJCOSX with def2/J | none stated |
| Embedding | electrostatic | electrostatic |
| Classical region | periodic box, 55680 atoms | 25 Å sphere, outer 5 Å fixed |

## X.3.2 Reaction path

I optimised the reactant and product under the quantum-mechanical/molecular-mechanical
treatment of X.3.1, then drove a restrained scan between them along the reaction
coordinate. The coordinate was defined as the difference between the breaking C4–O3 and
forming C6–C1 distances, following the definition used by Claeyssens *et al.* (2011) for
the corresponding atoms. The scan ran over 20 windows, restraining both distances toward
target values from 1.464 and 3.251 Å to 2.400 and 1.540 Å.

The scan did not resolve the transition-state region. Between the eleventh and twelfth
windows the achieved coordinate moved from −0.947 to +0.416 Å, a discontinuity of
1.363 Å against a mean window step of 0.084 Å on the reactant branch, and no structure
was obtained between those values. Claeyssens *et al.* locate the transition state at a coordinate of −0.5 Å,
which falls inside the unsampled interval. The achieved coordinate also trailed its target
progressively on the reactant branch, reaching a lag of 0.554 Å at the eleventh window,
because the scan restrained two interatomic distances and the coordinate is a derived
quantity that no restraint acted on directly.

The scan therefore served to establish that the endpoints are connected and that a
one-dimensional restrained scan cannot locate the saddle for this reaction. Transition
states were located by the nudged elastic band method described in X.3.3.

This differs from the approach of Claeyssens *et al.*, who generate profiles by restrained
optimisation along the reaction coordinate in steps of 0.2 Å, reduced to 0.1 Å near the
transition state, with the quantum and classical systems fully optimised at each step over
a range of −2.2 to +2.2 Å, and who take the highest point of the resulting profile as an
approximate transition-state structure.

## X.3.3 Transition-state location and characterisation

I located the transition state in two stages. A climbing-image nudged elastic band was
first converged over an active region of 1959 atoms, comprising the substrate, the
surrounding protein residues and 139 mobile water molecules, with the remainder of the
system held fixed. The converged climbing image was then refined by a transition-state
optimisation over a reduced active region of 102 atoms, which permits an analytical
Hessian. The quantum region was the same 24 atoms in both stages.

**The two active regions are distinct and their energies are not interchangeable.** An
energy calculated over the 1959-atom region and one calculated over the 102-atom region
differ by the contribution of the atoms free to relax in one and fixed in the other. Every
energy reported in the Results chapter states the active region it belongs to.

Vibrational analysis of the refined structure returned two imaginary frequencies, at
−313.30 and −18.65 cm⁻¹. The first corresponds to motion along the reaction coordinate.
The second lies below the 50 cm⁻¹ threshold conventionally applied to distinguish a
reaction-coordinate mode from numerical artefact. Agbaglo *et al.* (2024) attribute small
imaginary frequencies of this kind to partial Hessians computed with frozen Cartesian
coordinates, and describe a modified Hessian treatment in which the phenomenon does not
arise; the reduced-region optimisation used here holds all atoms outside 102 fixed, which
is the condition they identify. The structure is therefore a first-order saddle point with
respect to the reaction coordinate, carrying one additional low-frequency artefact of the
partial Hessian.

I confirmed that the saddle connects the intended minima by following the intrinsic
reaction coordinate in both directions from the converged structure. I then optimised the
reactant and product within the same 102-atom active region, so that all three stationary
points are defined over a common region.

`[FIGURE X.3.1 — reaction-path and transition-state flowchart. Topology from X.1.7 and
selected frames from X.2.3 into the bridge conversion, then endpoint optimisation and the
restrained scan, the scan branch terminating at the unsampled interval; then the
climbing-image band over 1959 atoms, the reduced-region refinement over 102 atoms, the
vibrational analysis and the intrinsic reaction coordinate. Annotate the scan branch with
the 1.363 Å discontinuity and the two band boxes with their active-region sizes.]`

`[FIGURE X.3.2 — restrained scan profile. Energy against achieved reaction coordinate from
scan_energies.tsv, with the target coordinate from step16_scan_targets.txt overlaid to
show the lag, the unsampled interval shaded, and a vertical line at the −0.5 Å transition
state of Claeyssens et al. Caption should state that the discontinuity spans the
literature transition-state position and that no barrier is read from this profile.]`

## Pipeline of Operations

- `step13a_bridge.sh` — converts the AMBER topology of X.1.7 into the ORCA force-field
  format. Outputs `05_qmmm/13_bridge/complex_solvated.ORCAFF.prms`,
  `05_qmmm/13_bridge/convff_A_single.log`.
- `step13b_smoketest.sh` — single point at a low-cost composite Hartree–Fock method,
  confirming the bridge. Outputs `05_qmmm/13_smoketest/smoke.inp`, `smoke.out`.
- `step14_dft_sp.sh` — single point at the production level of theory on the same
  geometry. Outputs `05_qmmm/14_dft_sp/dft_sp.inp`, `dft_sp.out`.
- `step15b_reactant_opt.sh`, `step17_product_opt.pbs` — endpoint optimisations. Outputs
  `05_qmmm/15_reactant_opt/reactant_opt.pdb`, `.property.txt`, and the corresponding
  product files under `05_qmmm/17_product_opt/`.
- `step16_restraint_scan_driver.sh`, `step16_restraint_scan.pbs`,
  `step16_scan_targets.txt` — the 20-window restrained scan and its target list. Outputs
  `05_qmmm/16_scan/scan_energies.tsv`, `05_qmmm/16_scan/scan_bonds_all.tsv`.
- `stage_tsguess.sh`, `band_status.sh`, `reference_inputs/neb_nebci.inp` — transition-state
  guess staging, band monitoring and the band input template.
- `step18` scripts — the climbing-image band over the 1959-atom active region, the
  reduced-region transition-state optimisation over 102 atoms, the vibrational analysis,
  the intrinsic reaction coordinate, and the reduced-region endpoint optimisations.
  Outputs `05_qmmm/18_nebts/nebts_NEB-CI_converged.activeRegion.xyz`,
  `05_qmmm/18c_reduced_region/ts_reduced.out`, `ts_reduced.pdb`, `active_reduced.txt`,
  `ts_frequencies_summary.txt`, `05_qmmm/18d_irc/irc_IRC_F.pdb`, `irc_IRC_B.pdb`,
  `05_qmmm/18e_reactant_reduced/reactant_reduced.out`,
  `05_qmmm/18f_product_reduced/product_reduced.out`.
- `step13_verify.py`, `step15_17_verify.py`, `step18_verify.py` — verification scripts.
  Outputs the corresponding `.out` files under each step's `checks/`.

---

# Sourcing gaps

None new.

# Inconsistencies found

1. **The frequency summary file is internally inconsistent and carries a superseded
   barrier.** `ts_frequencies_summary.txt` states one reaction-coordinate imaginary mode
   in its header while listing two in its body, and its footer quotes a barrier computed
   by differencing a 102-atom transition state against a reactant optimised over the band
   region. It must not be cited for either the mode count or the barrier until corrected.
2. **No written justification exists for the basis-set choice.** The functional matches
   Claeyssens *et al.*; the basis set and the dispersion correction do not, and neither
   repository records a reason.
3. **The intrinsic-reaction-coordinate output is not committed.** The forward and backward
   trajectories and structures are present, so the path is inspectable, but the
   calculation's own output file is absent and its convergence cannot be checked from the
   repository.
4. **The band energies are not committed**, so the band profile cannot be re-derived from
   the archived files. The converged climbing image and the trajectories are present.
5. **Both validation single points used frame 00820**, which is the pilot frame excluded
   from the reaction-path ensemble. This is correct for a validation calculation, from
   which nothing is quoted, and the Results chapter does not draw on it.

# Verification log

Performed by `checks/step13_verify.py`, `step15_17_verify.py` and `step18_verify.py`.
Twenty-seven, sixteen and twenty-nine checks respectively; all passed.

- Quantum-region indices resolved against the topology and confirmed to select exactly one
  substrate copy, the one the reaction-path frames are drawn from. The neighbouring copies
  abut it in the atom ordering, so an indexing error of one would select a region spanning
  two substrates while still returning 24 atoms, the correct charge and a converged
  self-consistent field.
- Charge and multiplicity confirmed as the substrate dianion and a closed-shell singlet.
- Zero link atoms confirmed from the output, consistent with a quantum region that cuts no
  bond.
- The three energy quantities a QM/MM calculation reports — the classical subsystem, the
  quantum subsystem and the total — confirmed present and internally consistent, with the
  total equal to the sum of the other two. The harvesting scripts were confirmed to select
  the total rather than the quantum subsystem alone.
- Scan completeness confirmed at 20 windows, all energies on the QM/MM total scale.
- The discontinuity and the coordinate lag reproduced from the committed scan output.
- The reduced active region confirmed at 102 atoms with the quantum region entirely
  inside it, and the band active region confirmed larger at 1959 atoms.
- All three reduced-region stationary points confirmed to have terminated normally and to
  lie in a common active region.
- Vibrational analysis confirmed to report exactly one imaginary frequency above the
  ghost-mode threshold and one below it.

# References

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2024) 'The influence of model
building schemes and molecular dynamics sampling on QM-cluster models: the chorismate
mutase case study', *Physical Chemistry Chemical Physics*, 26, pp. 12467–12482.
doi:10.1039/d3cp06100k.

Claeyssens, F., Ranaghan, K.E., Lawan, N., Macrae, S.J., Manby, F.R., Harvey, J.N. and
Mulholland, A.J. (2011) 'Analysis of chorismate mutase catalysis by QM/MM modelling of
enzyme-catalysed and uncatalysed reactions', *Organic & Biomolecular Chemistry*, 9,
pp. 1578–1590. doi:10.1039/c0ob00691b.

# X.1.7 System assembly and parameterisation

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Source files for every statement below are in `consolidated_works/09_tleap_build/`.

---

I combined the protonated protein with the three parameterised substrates into a single
structure and built the simulation topology, using `step09a_combine_pretleap.sh` and
`step09b_tleap_build.sh`. The assembly derived every count from its inputs, so
an upstream change surfaced here.

The combined structure comprised 6183 protein atoms and 72 substrate atoms across 384
residues. I estimated the formal charge before building, obtaining −3 for the protein
from its ionisable residues and the histidine assignment of X.1.5, and −6 for the three
substrates at −2 each, giving −9 for the complex.

I built the topology with tleap, sourcing the ff14SB protein force field, the general
AMBER force field for the substrates, and the TIP3P water model. The substrates entered
through the residue template and the additional parameters generated in X.1.6. The dry
topology carried 6255 atoms in 384 residues at a net charge of −9, matching the
pre-build estimate. I solvated the complex in a rectangular box with a 10 Å buffer and
neutralised it with sodium ions, giving 55680 atoms in 16865 residues, comprising 16472
water molecules and 9 sodium ions, at a net charge of zero in a box of
92.915 × 76.664 × 95.443 Å.

**Table X.1.7a.** Composition through assembly and solvation. Source:
`outputs/step09a_pretleap_audit.txt`, `outputs/step09b_topology_audit.txt`.

| Quantity | Value |
|---|---|
| Protein atoms | 6183 |
| Substrate atoms (3 × 24) | 72 |
| Combined atoms | 6255 |
| Residues | 384 |
| Protein formal charge | −3 |
| Substrate charge (3 × −2) | −6 |
| Complex charge | −9 |
| Sodium ions added | 9 |
| Chloride ions added | 0 |
| Solvated atoms | 55680 |
| Water molecules | 16472 |
| Solvated net charge | 0.000 |
| Box dimensions (Å) | 92.915 × 76.664 × 95.443 |

The build recipe follows Agbaglo *et al.* (2024), who parameterise the same enzyme with
ff14SB and the general AMBER force field, solvate in a cubic 10 Å box of TIP3P water, and
neutralise with nine sodium ions. The agreement extends to the ion count.

tleap reported 28 warnings and no errors. I classified each against a defined set rather
than inspecting them individually: 12 were terminal residue renames, 2 reported the
non-zero charge of the unsolvated complex, 1 was a unit-check note, and 14 reported close
contacts over 7 distinct atom pairs. No warning fell outside the classification. The close
contacts are examined in X.2, where the minimisation and equilibration that relieve them
are described.

**Table X.1.7b.** tleap diagnostic classification. Source:
`outputs/step09b_warning_classification.tsv`.

| Class | Count |
|---|---|
| Terminal residue renaming | 12 |
| Non-zero charge of the dry complex | 2 |
| Unit-check note | 1 |
| Close contact reports | 14 |
| Unclassified | 0 |

## Pipeline of Operations

- `step09a_combine_pretleap.sh` — merges the protonated protein of X.1.5 with the three
  frozen substrates of X.1.6, estimates the formal charge of each component, and confirms
  that every substrate copy matches the parameterised template in atom count and order.
  Outputs `03_amber/tleap_build/complex_for_tleap.pdb`,
  `00_admin/step09a_pretleap_audit.txt`, `00_admin/step09a_charge_estimate.tsv`,
  `00_admin/step09a_cha_template_check.tsv`, `00_admin/step09a_residue_inventory.tsv`.
- `step09b_tleap_build.sh` — writes the tleap input, builds the dry and solvated
  topologies, and audits both against the tleap log. Outputs
  `03_amber/tleap_build/tleap_complex.in`, `03_amber/tleap_build/complex_dry.prmtop`,
  `03_amber/tleap_build/complex_dry.inpcrd`, `03_amber/tleap_build/complex_solvated.prmtop`,
  `03_amber/tleap_build/complex_solvated.inpcrd`,
  `03_amber/tleap_build/complex_solvated.pdb`,
  `00_admin/step09b_topology_audit.txt`, `00_admin/step09b_warning_classification.tsv`,
  `00_admin/step09b_prmtop_residue_counts.tsv`, `00_admin/step09_fingerprints.txt`,
  `00_admin/sha256_step09_tleap_build.txt`.
- `step09_verify.py` — reproduces the charge arithmetic end to end, confirms atom-count
  conservation through assembly, reads the force fields actually sourced from the tleap
  input, and resolves each close contact against the graft provenance of X.1.2 and the
  contact set of X.1.4. Outputs `checks/step09_verify.out`.

`[FIGURE X.1.7 — assembly and solvation schematic. Three input nodes (protonated protein
from X.1.5, three parameterised substrates from X.1.6, force-field parameters from X.1.6)
converging on step09a, then step09b branching to the dry and solvated topologies.
Annotate the assembly arrow with 6183 + 72 = 6255 atoms and the solvation arrow with
6255 → 55680 atoms and −9 → 0 charge. Render from the file list in the Pipeline of
Operations; no data file required.]`

---

# Sourcing gaps

None new. The force-field and solvation choices follow Agbaglo *et al.* (2024), already
held.

# Inconsistencies found

1. **No as-run script archive exists for this step**, as for X.1.6. Neither script is
   archived alongside its outputs, unlike those of X.1.1 to X.1.5.
2. **The bucket counts sum to 29 against a reported total of 28 warnings.** The extra is
   the unit-check note, which the classifier counts separately from the warning stream.
   Not an error, but the two numbers count slightly different things.
3. **The solvated PDB saturates at 9999 residues**, so water and ion counts must be read
   from the topology rather than from the structure file. The build script states this and
   the audit follows it.
4. **Topology files embed a wall-clock date**, so their raw checksums differ between
   otherwise identical builds. `step09_fingerprints.txt` provides a date-insensitive
   fingerprint for comparison.

# Verification log

Performed by `checks/step09_verify.py`. Thirty-six checks, thirty-six passed.

- Charge arithmetic reproduced end to end: −3 protein plus −6 substrate gives −9, matching
  the dry topology, and the sodium count equals the magnitude of that charge.
- Atom-count conservation confirmed: 6183 + 72 = 6255, equal to the dry topology atom
  count and to the pre-tleap structure.
- All three substrate copies confirmed to match the parameterised template in atom count
  and order, with no duplicate, missing or extra atoms.
- Force fields read directly from the committed tleap input: ff14SB, the general AMBER
  force field, and TIP3P, with a 10 Å solvation buffer and sodium-only neutralisation.
- Solvated system confirmed: topology and coordinates agree on atom count, a periodic box
  is present, the net charge is integral and zero, and all three substrates survive
  solvation.
- All 28 tleap warnings classified, none unclassified, no errors and no fatal errors.
- No barrier is quoted in this subsection, so no barrier re-derivation applies.

# References

Agbaglo, D.A., Summers, T.J., Cheng, Q. and DeYonker, N.J. (2024) 'The influence of model
building schemes and molecular dynamics sampling on QM-cluster models: the chorismate
mutase case study', *Physical Chemistry Chemical Physics*, 26, pp. 12467–12482.
doi:10.1039/d3cp06100k.

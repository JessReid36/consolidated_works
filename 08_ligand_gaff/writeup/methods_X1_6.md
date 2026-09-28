# X.1.6 Substrate parameterisation

> **Drafting note (not part of the thesis text).** Follows `WRITEUP_PROMPT.md` §2–§3.
> Source files for every statement below are in `consolidated_works/08_ligand_gaff/`.

---

I parameterised chorismate with the general AMBER force field, deriving partial charges by
the AM1-BCC scheme and assigning atom types and missing parameters with Antechamber and
parmchk2. The derivation was split across two machines. AM1-BCC requires `sqm`, which does
not run on the HPC, so `step08a_am1bcc_charges.sh` derived the charges on a local
workstation under a self-contained AmberTools installation, and only the charge file was
transferred. `step08b_ligand_gaff.sh` then performed GAFF typing and parameter checking on
the HPC in read-charge mode, which requires no `sqm`.

I parameterised one canonical substrate template and applied it to all three copies. The
three frozen substrates share atom names and bond topology exactly, each carrying 24 atoms
in 24 bonds at a net charge of −2, which licenses the single parameter set.

The AM1-BCC charges as written summed to −1.997000 at the four decimal places of the mol2
format, and I renormalised them to exactly −2 by distributing the −0.003 residual evenly
over the 24 atoms, adding −0.000125 e to each. The renormalised charges were carried into
the GAFF typing unchanged.

Antechamber assigned ten GAFF types across the molecule — c, c2, c3, ce, h1, ha, ho, o, oh
and os — leaving no Tripos types in place and preserving the atom names and connectivity of
the template. It rewrote the bond orders of the two carboxylate groups from aromatic to
single; GAFF represents carboxylate delocalisation through equal partial charges and
identical oxygen types rather than through bond order, so the rewrite is the expected
treatment and the connectivity is unchanged.

parmchk2 generated the missing parameters and reported no parameters requiring attention
and none missing. Four parameters were supplied by analogy, and their penalty scores span
two orders of magnitude.

**Table X.1.6a.** Parameters supplied by analogy, with the penalty score parmchk2 assigned
to each substitution. Source: `outputs/cha.frcmod`, `checks/step08_verify.out`.

| Parameter | Type | Substituted by | Penalty |
|---|---|---|---|
| `c2-ce-os-c3` | dihedral | `X-c2-os-X` | 237.0 |
| `c -ce-os-c3` | dihedral | `X-c2-os-X` | 237.0 |
| `ce-o -c -o` | improper | general improper | 3.0 |
| `c2-ce-os` | angle | `c2-ce-oh` | 2.0 |

The two dihedrals carrying a penalty of 237.0 both span the `os`–`c3` pair. In the typed
substrate, O3 is the only atom of type `os` and C4 is of type `c3`, and the reaction
coordinate defined in X.2 identifies C4–O3 as the bond that breaks during the
rearrangement. The least reliable parameters in the substrate force field are therefore the
torsions about the reacting bond. The substrate occupies the quantum region throughout the
reaction-path work of X.3 onward, so these parameters do not enter any barrier. They do
govern the classical molecular dynamics of X.2, from which the reaction-path starting
structures are drawn, and the near-attack geometry that selection depends on is sensitive
to this torsion.

`[DECISION REQUIRED — whether to test the effect of these parameters on the sampled
ensemble. The test is to compare the C4–O3 torsion distribution over the production
trajectory against the same torsion in the QM/MM optimised reactants, which requires only
analysis of existing trajectories. Agreement would establish the parameterisation as
adequate for sampling, which is all it is required to do. Disagreement would not
invalidate any barrier but would constitute a sampling caveat for X.2 and for the ensemble
results.]`

## Pipeline of Operations

- `step08a_am1bcc_charges.sh` — derives AM1-BCC charges on a local workstation, extracts
  them in canonical atom order, renormalises to exactly −2, and records the derivation.
  Held in `local_workstation/` in the code repository. Outputs `cha_am1bcc.mol2`,
  `charges_am1bcc.dat`, `antechamber_am1bcc.log`, `00_admin/step08a_provenance.txt`.
- `step08b_ligand_gaff.sh` — verifies that the three substrate copies share names and
  topology, runs Antechamber in read-charge mode for GAFF typing, runs parmchk2 for the
  missing parameters, audits the result and checksums it. Outputs
  `03_amber/ligand_gaff/cha_gaff.mol2`, `03_amber/ligand_gaff/cha.frcmod`,
  `03_amber/ligand_gaff/antechamber_cha_gaff.log`,
  `00_admin/step08b_ligand_identity_check.tsv`, `00_admin/step08b_cha_gaff_audit.tsv`,
  `00_admin/step08_method_note.txt`, `00_admin/sha256_step08_ligand_gaff.txt`.
- `step08_verify.py` — confirms the charge template is the frozen substrate from X.1.5,
  reproduces the renormalisation atom by atom, confirms read-charge mode preserved the
  charges exactly, and ranks every substituted parameter by penalty score. Outputs
  `checks/step08_verify.out`.

`[FIGURE X.1.6 — the typed substrate, rendered from cha_gaff.mol2, atoms labelled with
their GAFF type and coloured by partial charge. Highlight O3 and C4 and the torsions listed
in Table X.1.6a. Caption should state that these are the reacting bond and the parameters
supplied by weakest analogy.]`

---

# Sourcing gaps

None new. The reference direction for Antechamber and GAFF is Agbaglo *et al.* (2024),
already held.

# Inconsistencies found

1. **Two declared run logs are absent from the repository.** `parmchk2_cha_gaff.log` and
   `antechamber_am1bcc.log` are written by the scripts but not committed. Neither carries
   a result; the parmchk2 outcome is fully recorded in `cha.frcmod` and the audit.
2. **No as-run script archive exists for step 08.** Neither half is archived alongside its
   outputs, unlike every script in X.1.1 to X.1.5.
3. **The toolchain is split and the split is untested.** Charges were derived under
   AmberTools 19 and typing performed under Antechamber 22.0. Both versions are recorded
   in the provenance, but no committed output tests whether the two agree on typing for
   this molecule.
4. **Zero parameters requiring attention is not the same as zero concern.** parmchk2
   reported no attention flags, which the audit records faithfully. The penalty scores,
   which are not summarised anywhere in the committed outputs, show two substitutions two
   orders of magnitude weaker than the rest.

# Verification log

Performed by `checks/step08_verify.py`. Twenty-eight checks, twenty-eight passed.

- Charge template confirmed as the frozen substrate from X.1.5 by SHA-256.
- Charge file carries 24 values summing to −2.00000000.
- Renormalisation reproduced: raw mol2 charges sum to −1.997000, per-atom shift
  −0.00012500 e, spread across atoms 1.11e-16 e.
- Read-charge mode confirmed: typed charges equal the charge file to 0.00e+00 e.
- Atom names and bond connectivity preserved; four carboxylate bond orders rewritten.
- Ten GAFF types assigned, no Tripos types surviving, type set matching the audit.
- Four substituted parameters ranked by penalty: 237.0, 237.0, 3.0, 2.0.
- Three-copy identity confirmed for names, topology and net charge.
- All three frozen checksums reproduce.
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

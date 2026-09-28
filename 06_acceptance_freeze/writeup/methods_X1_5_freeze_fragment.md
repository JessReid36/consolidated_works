# Fragment — Methods X.1.5, acceptance freeze

> Not a standalone subsection. X.1.5 covers the acceptance freeze and the protonation-state
> assignment together, and is assembled in `07_protonation_hpp/writeup/`. This fragment is
> the freeze half.

---

I froze the accepted pre-protonation complex using `step06_freeze_preprotonation.sh`. The
substrate geometries inherited repeated generic atom names from their source, which tleap
cannot resolve, so I renamed every substrate atom uniquely before freezing. The renaming
preserved coordinates, atom types, partial charges and bond connectivity; each substrate
retained 24 atoms in 24 bonds at a net charge of exactly −2, the value assigned to the
quantum region in X.3. The renaming tables are recorded per substrate in
`cha_a_renaming.tsv`, `cha_b_renaming.tsv` and `cha_c_renaming.tsv`.

The combined complex comprised 3030 protein atoms and 72 substrate heteroatoms across
three substrate residues, with no residue carrying duplicate atom names. I checksummed the
four accepted files — the combined structure and the three renamed substrates — in
`sha256_step06_accepted_preprotonation.txt`.

**Table X.1.5a.** Substrate audit after renaming. Source:
`outputs/step06_ligand_mol2_audit.tsv`.

| Substrate | Atoms | Heavy atoms | Bonds | Net charge |
|---|---|---|---|---|
| `cha_a` | 24 | 16 | 24 | −2.000000 |
| `cha_b` | 24 | 16 | 24 | −2.000000 |
| `cha_c` | 24 | 16 | 24 | −2.000000 |

---

# Verification log (freeze half)

Performed by `checks/step06_verify.py`. Forty-three checks, forty-three passed.

- Renaming confirmed as a pure relabelling for all three substrates: coordinates identical
  to 0.00e+00 Å, atom types and charges unchanged, bond blocks identical, renaming tables
  consistent with both input and output.
- Atom names unique after renaming: 24 distinct names per substrate.
- Net charge exactly −2.000000 in all three.
- Combined complex reproduced: 3030 + 72 = 3102 records, 3 substrate residues, 0 with
  duplicate names.
- Protein coordinates unchanged from X.1.3 at 0.00e+00 Å; substrate coordinates unchanged
  from X.1.4 to within PDB three-decimal rounding.
- All four frozen checksums reproduce.

# Inconsistencies found (freeze half)

1. **The acceptance decision text names a superseded residue list.**
   `step06_acceptance_decision.txt` is generated as a literal block inside
   `step06_freeze_preprotonation.sh`. It enumerates the adjacent-chain contributors to the
   active site as Arg63′, Lys60′, Thr74′, Cys75′, Phe57′ and Ala59′, omitting Val73′,
   which X.1.4 establishes as a contact at 3.49–3.59 Å in all three sites. The frozen
   structure is unaffected; only the text is out of date.

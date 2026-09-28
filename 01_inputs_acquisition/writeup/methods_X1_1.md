# X.1.1 Acquisition and provenance locking of raw inputs

> **Drafting note (not part of the thesis text).** The house-style exemplar
> `lategan_computational_2025.pdf` is not held. This subsection follows the written
> specification in `WRITEUP_PROMPT.md` §2–§3. Software names and versions are given here
> only where they identify a source file; per §2 they are consolidated in the Computational
> Details section rather than repeated in each methods subsection. Source files for every
> statement below are in `consolidated_works/01_inputs_acquisition/`.

---

I retrieved the protein scaffold as PDB entry 2CHT, the crystal structure of the
monofunctional chorismate mutase from *Bacillus subtilis* in complex with a
transition-state analogue (Chook, Ke and Lipscomb, 1993). I retrieved PDB entry 1DBF, the
1.30 Å structure of the same catalytic homotrimer, as the donor for the terminal repair
described in X.1.2 (Ladner *et al.*, 2000). I took the chorismate substrate coordinates as
three geometries pre-aligned to crystallographic analogue sites, from the CP2K examples
repository at commit `d71583a411a0a7d8557e94a24b87a3ceb87ec34e`. Both retrievals were
performed by `step01_download_inputs.sh`, which recorded the retrieval timestamp and the
upstream commit in `step01_provenance.txt` and the role of each file in
`step01_download_manifest.tsv`.

2CHT resolved twelve chains, A to L, carrying 13529 protein atoms and 1815 heteroatoms.
The only non-water heteroatom species present was the analogue, at 204 atoms distributed
one copy per chain across twelve sites. 1DBF resolved three chains, A to C, carrying 3156
protein atoms. Each of the three substrate geometries carried 24 atoms in 24 bonds, of
elemental composition C₁₀H₈O₆ and net formal charge −2, which set the charge assigned to
the quantum region in X.3.

`step01b_ligand_tsa_registration.sh` registered each substrate geometry to the analogue
site against which it had been aligned, by nearest molecular centre, assigning `liga.mol2`
to site K210 at a centre offset of 0.161 Å, `ligb.mol2` to site J212 at 0.160 Å, and
`ligc.mol2` to site L211 at 0.203 Å. The nearest alternative site lay at 22.272, 22.207 and
22.184 Å respectively.

I locked the identity of all five raw inputs by SHA-256 checksum in
`sha256_step01_inputs.txt`, since deposited structures may be re-versioned and repository
commits may be rewritten.

**Table X.1.1a.** Raw inputs, their role, and their recorded SHA-256 (first sixteen
characters). Source: `outputs/step01_download_manifest.tsv`,
`outputs/sha256_step01_inputs.txt`, `outputs/step01_provenance.txt`.

| File | Role | Source | SHA-256 (truncated) |
|---|---|---|---|
| `2cht_raw.pdb` | Protein scaffold | RCSB, PDB 2CHT | `b938020bbbcbe301` |
| `1dbf_raw.pdb` | Terminal-repair donor | RCSB, PDB 1DBF | `59183f711ea35c1b` |
| `liga.mol2` | Chorismate reference A | CP2K examples, `d71583a4` | `fd19c40914739e57` |
| `ligb.mol2` | Chorismate reference B | CP2K examples, `d71583a4` | `31800dab8780e66a` |
| `ligc.mol2` | Chorismate reference C | CP2K examples, `d71583a4` | `a4b8584c16f4237c` |

**Table X.1.1b.** Registration of each substrate geometry to its crystallographic analogue
site, with the distance to the nearest alternative site. Source:
`outputs/ligand_tsa_registration.tsv`, `checks/step01_verify.out`.

| Ligand | Registered site | Centre offset (Å) | Nearest alternative (Å) |
|---|---|---|---|
| `liga.mol2` | K210 | 0.161 | 22.272 |
| `ligb.mol2` | J212 | 0.160 | 22.207 |
| `ligc.mol2` | L211 | 0.203 | 22.184 |

Claeyssens *et al.* (2011) modelled the enzyme-catalysed reaction from the same
deposition with the analogue bound, substituting chorismate for the analogue. The
procedure adopted here departs from theirs in two respects. The terminal repair described
in X.1.2 has no counterpart in their reported procedure, and follows Agbaglo *et al.*
(2024), who repaired the same deposition from the same donor structure. The substrate
geometries used here were pre-aligned reference structures, whereas Claeyssens *et al.*
described their substitution by reference to earlier work
`[PLACEHOLDER — needs sourcing: Claeyssens et al. (2011) cite their references 9 and 10
for the substitution procedure, being Ranaghan, Ridder, Szefczyk, Sokalski, Hermann and
Mulholland, Mol. Phys. 2003, 101, 2695–2714 and Org. Biomol. Chem. 2004, 2, 968–980.
Neither is held and no DOI is available for either. See Sourcing gaps.]`. The analogue
carried seventeen heavy atoms against sixteen in the chorismate dianion, so the
superposition performed in X.1.4 was not an atom-for-atom correspondence under either
procedure.

## Pipeline of Operations

- `step01_download_inputs.sh` — retrieves PDB entries 2CHT and 1DBF from the RCSB and the
  three chorismate reference geometries from the CP2K examples repository at a pinned
  commit, reports a structural inventory of both depositions, and records role assignments,
  checksums and retrieval provenance. Outputs `01_inputs/structures/2cht_raw.pdb`,
  `01_inputs/structures/1dbf_raw.pdb`, `01_inputs/ligands/liga.mol2`,
  `01_inputs/ligands/ligb.mol2`, `01_inputs/ligands/ligc.mol2`,
  `00_admin/step01_download_manifest.tsv`, `00_admin/sha256_step01_inputs.txt`,
  `00_admin/step01_provenance.txt`.
- `step01b_ligand_tsa_registration.sh` — registers each chorismate reference geometry to
  the crystallographic analogue site against which it was aligned, by nearest molecular
  centre. Outputs `00_admin/ligand_tsa_registration.tsv`.
- `step01_verify.py` — recomputes all five checksums, reproduces the structural inventory
  and the analogue site list, confirms the elemental composition and net charge of each
  substrate geometry, and measures each registration against the nearest alternative site.
  Outputs `checks/step01_verify.out`.

`[FIGURE X.1.1 — schematic of input acquisition and provenance locking. Three source nodes
(RCSB 2CHT, RCSB 1DBF, CP2K examples at commit d71583a4) feeding step01, with step01b as a
downstream branch producing the site registration. Annotate each source with its role and
truncated SHA-256 from Table X.1.1a, and the registration branch with the three centre
offsets from Table X.1.1b. Render from the file list in the Pipeline of Operations; no data
file required.]`

---

# Sourcing gaps

1. **Ranaghan *et al.* (2003, 2004)** — the substitution procedure cited by
   Claeyssens *et al.* (2011) as their references 9 and 10. Neither PDF is held and no DOI
   is available for either. *Mol. Phys.* **101**, 2695–2714 and *Org. Biomol. Chem.* **2**,
   968–980.
2. **Ladner *et al.* (2000)** — DOI `10.1107/S0907444900004625`, read from the `JRNL DOI`
   record of the committed `1dbf_raw.pdb`. The DOI is held; the PDF is not.
3. **Bartlett and Johnson (1985)** — the transition-state analogue. Required only if the
   analogue is cited by name. Given by Claeyssens *et al.* (2011) as *J. Am. Chem. Soc.*
   **107**, 7792–7793; no DOI held.
4. **`lategan_computational_2025.pdf`** — house-style exemplar, not held.

# Inconsistencies found

1. **The acquisition script promises a check it does not perform.** The comment in
   `step01_download_inputs.sh` directs that elemental composition and net charge be
   confirmed from the mol2 before that charge is used to set the quantum-region charge. The
   code reports only the molecule name, atom count and bond count. No committed output
   contained the elemental or charge confirmation. The confirmation is implemented in
   `checks/step01_verify.py` and passes for all three geometries.
2. **The registration has no ambiguity guard.** `step01b_ligand_tsa_registration.sh`
   accepts the nearest analogue centre unconditionally, without testing the margin to the
   next nearest. The margin is measured by `checks/step01_verify.py` rather than by the
   acquisition script.
3. **The archived and working copies of both scripts differ.** The copies held in
   `00_admin/` differ from those in the code repository by a single added line setting
   thread limits for OMP, OpenBLAS and MKL. The difference cannot affect any output, and is
   recorded because the same archive has been found materially stale at a later step.
4. **Retrieval provenance does not record a deposition version.**
   `step01_provenance.txt` records a retrieval timestamp and the upstream repository
   commit, but no RCSB file version. Identity is pinned by the recorded checksums.

# Verification log

Performed by `checks/step01_verify.py` against the committed copies in
`01_inputs_acquisition/outputs/`. Twenty checks, twenty passed, none failed. Full output in
`checks/step01_verify.out`; the execution environment is recorded in `manifest.txt`.

- SHA-256 recomputed for all five raw inputs and compared against
  `sha256_step01_inputs.txt`: five of five match.
- 2CHT inventory reproduced: 13529 ATOM, 1815 HETATM, chains A–L, sole non-water
  heteroatom species TSA at 204 atoms.
- Twelve analogue sites enumerated: A203, B201, C202, D206, E204, F205, G209, H207, I208,
  J212, K210, L211.
- 1DBF inventory reproduced: 3156 ATOM, 499 HETATM, chains A–C.
- Substrate geometry confirmed for all three files: 24 atoms, 24 bo
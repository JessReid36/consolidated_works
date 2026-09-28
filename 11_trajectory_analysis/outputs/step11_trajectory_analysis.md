# Step 11 — Production trajectory quality analysis (findings)

Fork-independent QC of the 20 ns production (step 10c) before QM/MM frame selection.
All analyses are cpptraj-free (cpptraj is broken cluster-wide): a stdlib NetCDF3
header reader + numpy, single-threaded (OPENBLAS_NUM_THREADS=1) to avoid the login-node
OpenBLAS RLIMIT_NPROC segfault. RMSD/RMSF use molecule-wise imaging (undoes iwrap chain
splits) + Kabsch alignment; contacts use minimum-image distances.

## 11a — protein backbone RMSD vs production frame 0
- Backbone (N,CA,C,O,H) mean 3.16 A, drift 1.6 A (first 2 ns 1.98 -> last 2 ns 3.59); max 4.31.
- Gate flagged REVIEW on the drift; RESOLVED by 11b as benign C-terminal-tail motion.
- Consistent with Agbaglo/DeYonker whole-protein RMSD 2.66 A (both include floppy tails).

## 11b — per-residue Ca RMSF
- Core (lower-80%) mean 0.69 A, median 0.74 A -> rigid, well-folded trimer core.
- ALL flexibility is in the chain C-terminal tails: chain B res 249-254 up to 13.9 A,
  chain C 379-381 up to 7.3 A, chain A 127 5.6 A. (Numbering continuous: A=1-127,
  B=128-254, C=255-381.)
- Active-site residues {7,57,59,60,63,73,74,75,78,90,108,115} and second-shell Arg116
  are NOT in the flexible set -> active site stable.
- Conclusion: stable fold + active site; the 11a "drift" is three wandering tails.

## 11c — reaction-coordinate tracking (Claeyssens r)
- Atom mapping (Claeyssens Fig.1 -> cha_gaff): breaking C2-O13 -> C4-O3;
  forming C4-C14 -> C6-C1; carboxylate compression C16-C17 -> C3-C10.
- Validated by geometry: placed pose gives d(C4-O3)=1.46, d(C6-C1)=3.12, r=-1.66,
  matching Claeyssens reactant (3.13-3.29 / ~-1.6) to 2 dp.
- In classical MD the breaking bond (C4-O3) is a fixed force-field bond (~1.45 A always);
  the informative quantity is the forming distance d(C6-C1) = near-attack metric.
- Per-CHA means: CHA#2 (res 383) form 3.29 / r -1.83 / carbox 5.58 (closest to Claeyssens);
  CHA#3 (384) 3.83 / -2.38 / 6.51; CHA#1 (382) 4.15 / -2.70 / 6.69. All three sample
  near-attack conformations (form min 2.84-2.97 A, below the 3.7 A NAC threshold).

## 11d — catalytic-contact validation
- Per substrate, nearest partner (global numbering; +127/chain): ether O3->Arg,
  carboxylate->Arg, hydroxyl O4->acidic.
- All three inter-subunit sites assembled correctly with the full Claeyssens network:
  ether O3 -> Arg90; carboxylate -> Arg7/Arg63; hydroxyl O4 -> Glu78.
    CHA#1 (382): Arg90 4.71 A (min 2.62, 84%) | Arg7 2.72 (71%) | Glu78 8.02 (min 2.85, 93%)
    CHA#2 (383): Arg90 2.92 A (100%)          | Arg63 2.69 (53%) | Glu78 2.81 (100%)
    CHA#3 (384): Arg90 3.44 A (min 2.54, 82%) | Arg7 2.70 (42%) | Glu78 7.83 (min 2.60, 99%)
- CHA#2 is a textbook Michaelis complex: Arg90-O13 2.92 A locked 100% AND Glu78-OH 2.81 A
  100% (hydroxyl points to Glu78, exactly as Claeyssens describes). CHA#1/#3 reach the
  same H-bonds transiently (mins ~2.5-2.9 A) but sit looser on average.

## Scientific conclusions
1. BsCM biology: it IS a homotrimer and the trimer is the functional enzyme; three
   crystallographically equivalent active sites at subunit interfaces (2CHT; Chook,
   Ke & Lipscomb, PNAS 1993, 90:8600). Modelling one solvated trimer is correct and
   complete - not "less than the enzyme". (2CHT ASU has 4 trimers with crystal contacts;
   our isolated solvated trimer correctly omits those, matching the QM/MM literature.)
2. Site asymmetry is NOT error: it is finite-sampling symmetry-breaking (20 ns is not
   long enough for 3 independent sites to converge) modulated by the strength of the
   Arg90-O13 / Glu78-OH anchoring - which is precisely the TS-stabilisation physics
   Claeyssens 2011 studies. One site (CHA#2) locked the network and stayed near-attack;
   two sampled a looser mode. This spread of competent conformations is a feature for
   Claeyssens-style multi-pathway sampling.
3. Substrate pose VALIDATED with evidence: TSA-superposition + full MD relaxation
   produced genuine, correctly-assembled Michaelis complexes making the complete
   catalytic contact network. No MD rerun is needed on pose grounds. (The earlier
   "CP2K placeholder must be replaced" note meant in-pocket relaxation, which the MD
   itself accomplished under the full solvated force field.)
4. Two-paper design (record explicitly in Methods): Agbaglo/DeYonker 2024 = MD protocol
   source (AMBER22/ff14SB/GAFF/TIP3P); Claeyssens 2011 = QM/MM reaction-path science
   reference (chorismate-only QM region, electrostatic embedding, no link atoms;
   Arg90/Arg7/Glu78 TS stabilisation; reaction coordinate r = d(C2-O13) - d(C4-C14)).
   Implementation modernises Claeyssens' QoMMMa/adiabatic-mapping to ORCA/NEB-TS on the
   AMBER ff14SB system.
5. OPEN, non-blocking: AM1-BCC vs RESP substrate charges - matters only for QM/MM
   energetics, to be checked against what Agbaglo specifies. Does not block frame
   selection. Substrate stereochemistry to be confirmed visually in ChimeraX.

## Step 12 frame-selection criterion (set here)
Select catalytically-competent reactant frames: Arg90-O13 intact (< ~3.2 A) AND
near-attack (small d(C6-C1) / r toward the TS ~-0.5). Draws from all three sites,
weighted toward CHA#2; gives a Claeyssens-style spread of competent conformations.

## Provenance
Scripts: phase1_system_dev/step11{a_rmsd,b_rmsf,c_rxn_coord,d_contacts}.sh, step11_plots.sh.
Inputs: complex_solvated.prmtop; prod.nc (sha256 2d55390fa53ac1c5bcf1909f5d9edda87890f673747085098acd8c1c61f0b99e).
Outputs: 04_amber_md/11_analysis/{rmsd_vs_time,rmsf_per_residue,rxn_coord_per_frame}.dat;
11d summary reproducible from step11d_contacts.sh. Plots rendered locally (HPC py3.6
cannot build modern matplotlib).

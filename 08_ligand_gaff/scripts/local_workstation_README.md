# Local-workstation steps (NOT run on the HPC)

Scripts in this folder run on a local Ubuntu workstation, not the cluster, because
they need tools the HPC lacks. Convention: anything directly under
phase1_system_dev/ runs on the HPC; anything in local_workstation/ runs on your PC.

- step08a_am1bcc_charges.sh - derives AM1-BCC charges for chorismate. AM1-BCC needs
  a working `sqm`, which the HPC lacks (missing libopenblas.so.0). Produces
  charges_am1bcc.dat, which is then scp'd to the HPC for step08b_ligand_gaff.sh.
  Needs a one-time conda AmberTools (micromamba) install; the script stops with
  guidance if run without it (e.g. on the HPC).

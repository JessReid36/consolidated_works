Stage manifests. NONE of these is the ensemble set.

The canonical manifest is ../selection_manifest.tsv, which lives at the TOP LEVEL of
05_qmmm in the results repository and holds 45 rows: the 44 ensemble frames plus frame
820, the excluded pilot.

The files here live under 05_qmmm/12_frame_selection/ and are intermediate products:

  step12a_selection_manifest.tsv      12 rows, the step12a distance-clause selection
  selection_manifest_nac_extra.tsv     3 rows, the three full-NAC frames added by step12b
  selection_manifest_eq30.tsv         30 rows, equilibration-cut selection
  selection_manifest_ext_only.tsv     30 rows, extension-trajectory selection
  selection_manifest_ext_renumbered.tsv  30 rows, the same renumbered
  *_angles.tsv                        the same sets with the NAC angle columns added

NOTE THE FILENAME COLLISION. The step12a output is also called selection_manifest.tsv in
its own directory. It is renamed here to step12a_selection_manifest.tsv so it cannot be
confused with the canonical file. Copying 05_qmmm/12_frame_selection/*.tsv into a flat
directory alongside the canonical manifest silently overwrites it. That happened once
while this step folder was being built.

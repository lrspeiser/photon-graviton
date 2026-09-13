# Legacy targeting-code provenance check

The preceding turn traced the infrared exceptions to historical target lists. This turn located the public apogeetarget source and tested whether a visible missing-value branch explains those exceptions.

## Source evidence

The public [apogeetarget repository](https://svn.sdss.org/public/repo/apogee/apogeetarget/) exposes tags dated 2011–2013 in its current index (revision 86573). Those tags do not establish which program/version produced our 2017 southern designs. Four relevant routines and the tag index have been cached with exact hashes and URLs.

`check_data_quality.pro` applies the near-infrared quality flags, proximity, contamination and extended-source cuts, explicitly deferring the mid-infrared error threshold to science selection. This is consistent with treating the two filtering stages separately.

`select_science.pro` checks whether the mean infrared error is finite while ignoring NaNs; `select_science2.pro` performs that guard without ignoring NaNs. Each omits the infrared-error threshold in its alternative branch. For an illustrative error array [0.05, missing, 0.15], ordinary averaging selects the branch without the threshold, whereas NaN-aware averaging selects the branch with it. The diagnostic reproduces those guard decisions in NumPy; it does not execute IDL or reproduce the survey pipeline.

## Tested against the available parent

Among all 44,957 parent entries with finite, nonsentinel targeting extinction, **none has a nonfinite selected infrared error**. Of these entries, 24,082 have errors above 0.1 mag. Thus the missing-error guard does not explain disabling the cut on this available subset. We have not proved that this subset is identical to either historical routine's input.

The 42 observed exceptions therefore remain unexplained by this code audit. Do not adopt the NaN branch as their cause, claim the 2017 pipeline was buggy, or relax the unseen parent selection based on it. The historical criteria parameters and actual southern implementation remain missing evidence.

## Research consequence

The prior candidate parent counts remain sensitivity cases, not calibrated inclusion probabilities. This result narrows a proposed explanation rather than settling selection. The next physical comparison can retain both selection cases as explicitly conditional sensitivity analyses while the actual selection and camera/field geometry are reconstructed. No stellar-density or companion-gravity conclusion is justified by this audit alone.

Run `python research_work/results/selection-code-provenance/run.py`. The result records code hashes, index revision, tag names, the illustrative branch check and the actual parent error counts. No held-out kinematics were read, no selection weights or physical parameters changed, and all six scientific objectives remain open. These are software/data-provenance checks, not new physical equations.

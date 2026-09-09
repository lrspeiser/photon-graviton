# Browser recovery result

Recovered original files from the ChatGPT Library associated with Galactic Star Distances on 9 September 2026.

The full-project download, photon-graviton-all-files.zip (25,593,665 bytes), has no ZIP central directory/end record. Its 261 consecutive local-file records nevertheless decompress completely and pass their individual CRC-32 and declared-size checks. These intact entries were extracted in staging, then compared with the included original SNAPSHOT_MANIFEST.json. No corrupt or partial payload was accepted.

Combining this recovery with earlier archive contents restores 271 of the 276 files listed in the original snapshot. All 271 match the recorded SHA-256; there are no mismatches. Five original snapshot files remain absent. The manifest describes the available research checkpoints, not every file ever created.

Restored into the local repository. The full archive added or updated 155 paths; the separate time_theory_revision_data.zip restored six more files with matching snapshot hashes. Conflicting README.md, reconstructed requirements.txt and the previous causal report were backed up under the local recovery workspace before replacing them with the downloaded originals. The original report is preserved unchanged at [archive/original-uploads/report(3).md](archive/original-uploads/report%283%29.md). This report describes the recovery stage; the later repository update includes the recovered sources and subsequent research.

## Major recoveries

- companion_causal_test/run.py, followup.py and circulation.py.
- DustPedia THEMIS/DL14 tables, protocol, source hashes, per-object energy budgets, cross-validation predictions, bootstrap/follow-up results and circulation results.
- companion_wave_test source, data and cluster transcription, plus companion_deposition_fit code and results.
- option3_test/run_test.py and protocol.md, and missing early branch sources.
- minimal_clock_interaction/derivation.md and other previously absent research folders.
- redshift_paper/temporal_redshift_paper.docx, authenticated against the snapshot described as latest paper v9.
- Original requirements.txt, restore_data.py and SNAPSHOT_MANIFEST.json.

## Verification

All three latest Python scripts parse. circulation.py ran successfully in a scratch copy and its parsed JSON output exactly matched the archived result during recovery. All three were subsequently reproduced in isolation; see the [baseline reproduction report](research_work/results/baseline/reproduction-report.md) and [current verification](research_work/publication-verification.json) for the documented numerical differences. Existing archived results remain intact.

## Remaining original snapshot files

- `unified_paper_v9/build.py`
- `unified_paper_v9/design_tokens.json`
- `unified_paper_v9/math_format.py`
- `unified_paper_v9/verification.json`
- `unified_paper_v9/version8_original.docx`

The remaining gaps are historical version-9 paper-building files; the finished version-9 paper is recovered. They do not include the three requested companion scripts. See [current original-file hashes](research_work/results/checkpoint/current-original-hashes.json) and [baseline manifest audit](research_work/results/baseline/manifest-audit.json). Earlier local recovery inventories are superseded by these checks and the explicit missing-file list above.

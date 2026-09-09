Temporal redshift paper: scientific data supplement

all_164_groups.csv contains every selected object, all three splits, observed labels, frozen primary predictions, and residuals. It is a disclosed data release, not a blinded input.

The original experiment is preserved unchanged in redshift_paper_sources/time-redshift-expanded.zip. Extract that archive to rerun its pipeline as described in the paper.

To regenerate this paper from the supplement root, run python redshift_paper/build_paper.py with numpy, scipy, matplotlib, python-docx and lxml installed. The script extracts the preserved source archive if needed.

The temporal optical mechanism is a conditional hypothesis developed after the empirical fit. No physical validation of that mechanism is claimed.

# Second force-source refinement in progress

Previous goal turn was progress: it completed the source and force comparisons and recorded failed gates. This turn diagnoses cancellation and verifies the local close-passage formulas, while launching the next force quadrature at unchanged physical parameters.

Active integrations: refine12.py1 uses session37174; refine12.py3 uses session21683. Each must contain72 completed records in orbits12-R1.json or orbits12-R3.json and reach terminal process status. Recheck these handles or authoritative process command lines before considering restart. Partial files are not completion evidence. The partial raw files are excluded from the current commit.

After both finish, run export12.py and document12.py. Inspect every orbit/age gate, retain failures, compare8/12 with unchanged potential/force thresholds, and then commit raw files and summaries. A failed age gate requires numerical repair before relying on its force value. Do not relabel a source-convergence failure as physical rejection or hide it by tuning normalization.

The contribution diagnostic is in force-conditioning.md and concentration.json. It does not prove that an individual large contribution is wrong. Analytic passage checks are in passage-checks.json. All nine goals remain incomplete; no observed-star fit, supply closure, self-gravity or successful joint redshift model is claimed.

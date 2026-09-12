# Angular refinement in progress

The preceding checkpoint changed authoritative results and was progress. This checkpoint adds the population-to-gravity derivation, a fixed refinement protocol, and two active integrations. It does not establish angular convergence.

Runs launched with run.py 8 and run.py 12. Tool session handles: 30388 (8) and 50979 (12). Authoritative completion requires a terminal process and exactly64/144 records respectively in grid8-L40.json and grid12-L40.json. These files are rewritten after each trajectory and are partial until then. Recheck live processes before restarting; a missing observation handle alone is not evidence of failure if the process is still running.

After both complete, run export.py, inspect every numerical gate and adjacent-grid comparison, report relative residence differences as well as absolute gates, and retain failures. Commit final raw data and summaries only after this inspection. The currently running raw JSON files are intentionally not part of this checkpoint commit.

Read population-to-gravity.md for the conditional derivation linking continuous moving injection to density, potential and force. It supplies the next calculation target, not a completed galaxy or redshift theory. All nine goals remain active.

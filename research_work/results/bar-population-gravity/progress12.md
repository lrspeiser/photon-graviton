# Second force-source refinement completed

Both processes are terminal: sessions37174 and21683 each completed72 trajectories. The complete raw files are orbits12-R1.json and orbits12-R3.json. Do not restart those jobs.

export12.py and document12.py evaluated the source8/source12 comparison. All208 combined orbit/age status checks pass;2/30 potential gates and19/30 force gates fail. Final central force changes remain large. No direction, failed comparison or physical parameter was removed or tuned.

The independent benchmark in ../stream-force-quadrature/report.md demonstrates the same potential/force discrepancy for a smooth straight-stream source and verifies an unsoftened integration remedy. Next adapt local close-passage/source-coordinate integration to the curved bar trajectories; benchmark success is not a completed bar solver. Full energy supply, self-gravity, redshift/timing and all nine goals remain incomplete.

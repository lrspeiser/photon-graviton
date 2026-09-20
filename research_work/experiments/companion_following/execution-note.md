# GF-1 execution note: unavailable ring balance

The first screen completed all 600 phenomenological chains, all 600
phenomenological rings and all 18 conservative chains. It stopped while
initializing the batch of conservative rings:
`ValueError: No positive circular balance for ['C-K5-E8']`.

This is a declared physical-fixture failure, not evidence for or against a
different radius or an elliptical orbit. The protocol explicitly requires
recording it as unavailable. No failed output is removed. The original
screen-v1 source manifest and three completed trajectory archives stay intact;
its empty conservative-ring directory records the interrupted stage.

A completion invocation checks each conservative circular balance separately
and runs only cases with a positive balance at the original radius. Its
combined index references the original results and records the unavailable
case with behavioral_pass=false and conservation_pass=false. An unavailable
initial condition is counted separately from an executed trajectory.

The audit and refinement readers use this combined index. Model equations,
timesteps, selection rules and pass thresholds are unchanged. No parameter is
retuned after seeing a failure. The initial command was source commit d45b775
and returned exit code 1; the completion source has its own committed manifest.

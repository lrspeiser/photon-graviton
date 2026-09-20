# MS-1 execution and refinement note

The first run completed every simulation and saved the audit and completion
marker, but returned exit1 when its console-only json.dumps encountered a
NumPy Boolean. Fix the print path without changing or rerunning that evidence.

Eight of60 first timestep comparisons fail. Keep them. Refine the same six
selected laws at dt0.005 and0.0025, plus the highest and lowest mean inward
positive-coupling cases. Keep the original absolute angle/position and relative
exposure tolerances. This is an exposed numerical diagnosis of sharp transport
boundaries, not a new physical parameter fit or new holdout.

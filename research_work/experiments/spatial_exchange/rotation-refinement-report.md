# SE-RR1: declared n64 rotation comparison passes

20 September 2026. Protocolcc5ba79; executablea3f976a. The n64 rotated
evolution passes its individual energy/angular/boundary gates and both
declared radius1.5 circulation comparisons. Independent endpoint, polynomial
interpolation and decision reconstruction passes34 checks; largest endpoint
discrepancy1.09e-19. The earlier n48 failure remains unchanged.

| Method | n48 rotation difference | n64 rotation difference | Declared gate |
|---|---:|---:|---|
|Linear,512 loop points|2.52575%|1.61535%|<2% and decreases: pass|
|Four-point cubic,512 points|2.26686%|1.31560%|<2% and decreases: pass|

Projected A-sector angular momentum changes4.75e-10 relative under this
rotation. Relative total-energy drift is5.28e-10 and angular-vector drift
5.90e-7; near-edge energy stays below1.51e-19 of initial energy. The same
equations, coefficients, excitation, box and timestep were used. Only grid
resolution changed from the failed rotated case; unrotated n64 was reused.

This is not accuracy at every radius. At radius2 the very small circulation
still changes approximately13.5% with linear and14.4% with cubic sampling
under rotation. There was no preregistered radius2 accuracy gate, and its
observed error must not be hidden by the radius1.5 pass. Further outer-field
resolution would be needed before interpreting that tail or tracing light
through it with a claimed precision. One rotation and one short excitation
also do not prove general rotational invariance, source-size independence or
long-term stability.

The current model has finite prescribed internal excitation and measured
field transfer. It still lacks ordinary-matter production/fuel and a fully
coupled matter/light evolution. No galaxy or cluster prediction is established.
The audit independently reconstructs endpoint ledgers and both interpolants;
intermediate extrema come from saved traces. Full initial/final arrays and
source/input hashes remain in rotation-refinement-v1.

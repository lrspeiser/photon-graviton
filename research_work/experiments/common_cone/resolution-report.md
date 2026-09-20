# SR-2 final report: consistency passes, accuracy gate narrowly fails

20 September 2026. All nine declared runs finish and pass their individual
numerical gates. The campaign fails the cross-discretization accuracy gate:
the second/fourth-order bends at n70 differ1.012759%, exceeding1%. Preserve
this result without rounding it into a pass or altering the threshold.

| Accuracy test | Measured | Requirement | Result |
|---|---:|---:|---|
| Fourth-order n56/n70 bend difference |0.621867%|<1%|Pass|
| Fourth-order refinement difference decreases |0.00025734 to0.00003382 rad|Decrease|Pass|
| Second/fourth order, n70 |1.012759%|<1%|Fail|
| Timestep refinement |1.783e-7%|<0.1%|Pass|
| Rotation direction error / bend |0.303080%|<1%|Pass|
| Earlier second-order replay maximum error |5.20e-18|<1e-10|Pass|

Independent reconstruction passes214 checks across all nine archives; maximum
independent energy discrepancy is1.78e-15. Maximum recorded ledger drift is
6.95e-11. The original declarations cover gradient/compact-deposition and
Fourier controls as well; these are numerical checks, not physical stability
or observational validation.

| Run | Signed bend (model radians) | Matter energy decrease | Angular residual | Wall / CPU seconds |
|---|---:|---:|---:|---:|
|Second order n28|-0.006540922|0.02419849|1.64e-4|8.50 /8.36|
|Fourth order n28|-0.006344769|0.02373587|1.61e-4|9.51 /9.17|
|Fourth order n40|-0.005728925|0.02311838|3.01e-5|205.15 /198.98|
|Fourth order n56|-0.005471586|0.02286095|1.01e-5|613.25 /595.67|
|Fourth order n70|-0.005437770|0.02282111|6.28e-6|1298.35 /1258.34|
|Second order n56|-0.005554723|0.02304862|1.19e-5|381.80 /371.19|
|Second order n70|-0.005492842|0.02295052|5.18e-6|803.72 /778.63|
|Fourth order n56, half dt|-0.005471586|0.02286095|1.01e-5|1439.48 /1399.83|
|Fourth order n56, rotated|-0.005474020|0.02283497|1.37e-5|755.06 /735.38|

The matter decrease is the archived initial-minus-final particle Hamiltonian
energy, including the photon, not an inferred rest-mass loss or a physical
stellar fuel rate. These timings were taken on the working machine during
other jobs and are not controlled performance benchmarks. Angular residuals
are the runner's canonical-vector drift diagnostic.

Signed bend is the change in the probe angle measured in the source plane
relative to its rotated initial direction. The unsigned3D bend is
atan2(|v_initial cross v_final|,v_initial dot v_final). It agrees with the
signed magnitude to printed precision in the unrotated cases. In the rotated
case it is0.005474044893 radians, while the signed planar bend is
-0.005474020316. The rotation comparison transforms the final probe velocity
back before comparing directions; its error is1.65833e-5 radians.

This is the earlier common-cone model with a finite-radius probe/source
treatment, not the later shrinking point-probe or local-rotor model. Grid
agreement cannot establish a point-ray interpretation. The wider fourth-order
stencil and pre-boundary clearance do not validate a strict numerical causal
front or long-time outgoing boundary. The largest grid still fails the stated
cross-method target. Further resolution or an independently declared numerical
change is needed before claiming that target precision.

No observational fit, sustained galaxy orbit or cluster-lensing prediction is
established. No old gravity curve is used as the success criterion. Established
finite differences, Hamiltonian adjoints and Fourier analysis are mathematical
tools credited in resolution-protocol.md. Raw states, traces, source hashes
and final failed summary remain in resolution-v1; the independent evidence
is resolution-audit.json. All twelve broader goal requirements remain active.

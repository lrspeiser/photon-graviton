# Closing the angular-kernel and transfer-rate calculation

13 September 2026. Consistency test of the earlier Gaussian overlap branch.

## Outcome

The very small observer moments from a kernel that tracks apparent source size cannot be combined freely with a constant transfer rate. In the existing Gaussian receiving-overlap model, making the kernel narrower reduces the rate and eventually strongly suppresses the forward transition through longitudinal momentum mismatch. At unchanged target density and coupling, accumulation saturates near the source rather than growing across millions of parsecs.

This rejects that particular combination of assumptions. It is not a universal exclusion of narrow scattering, collective interactions or companion transport. It also does not change the empirical one-third reference.

## Matched kernel and provenance

The earlier spatial-response calculation postulated a Gaussian transition-density overlap, with natural-unit length a (in inverse eV):

\[
|F(q)|^2=e^{-a^2q^2},\quad
q^2=\delta^2+2E(E-\delta)(1-\cos\theta).
\]

This is the previous heavy-store, negligible-recoil approximation. The Gaussian profile and operator are hypotheses; Fourier overlap and angular integration are known mathematics. In its narrow-angle limit,

\[
\langle\theta^2\rangle\simeq\frac1{a^2E(E-\delta)},\qquad
S\simeq\frac{3e^{-(a\delta)^2}}{8a^2E(E-\delta)}.
\]

S is the rate relative to the point-overlap version at the same photon energy. Match the angular second moment of the previous source-size prescription, theta_rms=d R_sun/s, using

\[
a(s)=\frac{s}{dR_{sun}\sqrt{E(E-\delta)}}.
\]

This is only moment matching: a Gaussian angular distribution is not the fixed-radius kick distribution used to derive d=0.157241775. Thus it does not prove the same empty-channel fraction or loading bias. Its purpose is to test the rate cost of comparable narrowing.

With E held at 2 eV, the resulting rate shape is

\[
S(s)\propto s^{-2}e^{-(s/s_c)^2},\qquad
s_c=\frac{dR_{sun}\sqrt{E(E-\delta)}}{\delta}.
\]

These equations are derived consequences of our earlier overlap ansatz, not a new fundamental law. Holding E fixed isolates the kernel effect; this is not a coupled energy-evolution simulation.

## Numerical scale

| Transfer delta | Coherence size at 1 pc | s_c | log10 rate suppression at 1 pc | log10 rate suppression at 1000 pc |
|---:|---:|---:|---:|---:|
| 1e-8 eV | 27.830 m | 0.70904 pc | -18.191 | -863887.64 |
| 2e-10 eV | 27.830 m | 35.4519 pc | -17.327 | -368.872 |

The coherence size scales linearly with s: about 27.8 km at 1000 pc. For a nonzero energy transfer even perfectly forward scattering has momentum mismatch delta/c in this approximation. Once a delta exceeds one, the Gaussian overlap strongly suppresses that transition. Increasing coherence cannot therefore narrow the angular kernel at no rate cost.

## Integrated accumulation diagnostic

Normalize the fractional loss coefficient at s_ref=1 pc to the previous reference alpha_ref=0.0002488993286 per Mpc, and keep coupling and target density unchanged farther out. In the frozen-energy probe,

\[
\alpha(s)=\alpha_{ref}(s_{ref}/s)^2
\exp[-(s^2-s_{ref}^2)/s_c^2].
\]

Its infinite effective path length is

\[
L_{eff}=s_{ref}\left[1-\sqrt\pi x\,\operatorname{erfcx}(x)\right],\quad
x=s_{ref}/s_c.
\]

| Transfer | Effective path after 1 pc | Diagnostic redshift exp(alpha_ref L_eff)-1 |
|---:|---:|---:|
| 1e-8 eV | 0.157845 pc | 3.929e-11 |
| 2e-10 eV | 0.951556 pc | 2.368e-10 |

The uniform-reference relation over the earlier 30.660139 Mpc distance instead gives about 0.00766048. Merely multiplying the local normalization would need factors about 1.94e8 or 3.22e7 for that one path, and would still leave a source-localized effect that saturates with distance. It would not restore the desired common distance dependence.

These redshift numbers are conditional frozen-energy diagnostics, not observations or fitted predictions. The overall point-model coupling is not established. Other source environments and evolving photon energies would need a full transport calculation.

## What could change

Maintaining a constant rate within this Gaussian prescription would require a compensating coupling-squared times target-density factor proportional to s^2 exp[(s/s_c)^2]. No such compensation has been derived, and source distance is not automatically a locally measurable interaction parameter. Even eliminating the longitudinal Gaussian factor would leave the inverse-square rate cost under this normalization.

Other momentum-space responses, coherent frequency conversion, different carriers or a kernel set by a local field remain candidates. Their angular profile, rate, energy and momentum exchange must be derived together. The next useful decision is to compare such complete interactions, rather than combining the favorable angular moments of one prescription with the rate of another. The previous source-size-tracking observer result should be read as a kinematic benchmark only.

## Verification

Ten distance/energy cases are evaluated. Four non-underflowing cases independently call the original exact angular integrator and check both suppression and angular moment against the narrow-angle expression to 1e-10. Two integrated rate shapes agree with the analytic erfcx expression to relative tolerance 1e-9. Logarithmic rates preserve extremely small values without silently interpreting underflow as exact zero. These checks establish consistency inside this branch, not its physical existence.

Files: `kernel-rate-closure.py`, `kernel-rate-closure-results.json`. Related: [spatial response](spatial-response-report.md), [observer transport](observer-transport-report.md).

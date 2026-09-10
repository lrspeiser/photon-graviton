# Does photon-to-gravity conversion produce our redshift?

2026-09-10. This tests a possible cause behind the existing redshift calculator; it does not replace that calculator. No galaxy rate was refitted, no new observations were scored and no holdouts were opened. [Protocol](protocol.md), [calculation](run.py), [comparison](comparison.csv), [results](results.json).

**Result:** stationary, phase-matched linear photon/graviton mixing can transfer exactly the required total energy while leaving surviving photons at their original frequencies. It also preserves the separation between shifted copies of a source event. That specialization is unsuitable as the sole cause of our accumulated redshift and event stretching. It remains a possible local gravitational-radiation source, subject to its own physical limitations.

## The candidate tested

The Gertsenshtein mechanism describes electromagnetic/gravitational mode mixing in a background magnetic field. Its published treatment includes coherence requirements and propagation effects; simply increasing the background field does not guarantee efficient conversion. Here we use the ideal phase-matched, stationary two-mode limit as a deliberately favorable spectral test, not a full neutron-star plasma calculation. [Primary derivation](https://arxiv.org/html/2301.02072v3)

**Known coupled-mode/unitary evolution mathematics:** with normalized photon and gravitational amplitudes, a constant real effective mixing coefficient mu gives

\[
i\frac{d}{ds}\begin{pmatrix}A_\gamma\\A_g\end{pmatrix}
=\begin{pmatrix}0&\mu\\\mu&0\end{pmatrix}
\begin{pmatrix}A_\gamma\\A_g\end{pmatrix}.
\]

For a pure incoming photon mode, the conversion fraction is sin^2(theta), where theta=mu*L. No project novelty is claimed for this matrix, and mu is not identified with alpha or a measured magnetic field. Equal mode frequencies allow the chosen normalization to track transferred energy.

**Diagnostic choice, not a physical fit:** set sin^2(theta)=1-exp(-alpha D), so the converted energy equals the amount required by our existing exponential photon-energy-loss rule. This removes inadequate conversion strength as an explanation for failure in this test. The rate has not been derived from a real environment.

## Same energy transfer, different spectrum

For the 100-million-light-year example, required transfer is about 0.7602% of incoming photon energy:

| Quantity | Stationary mode conversion | Desired partial-energy redshift |
|---|---:|---:|
| Photon energy remaining | 99.2398% | 99.2398% |
| Photon number remaining | 99.2398% | 100% |
| Mean frequency of surviving photons, relative to input | 100% | 99.2398% |
| Spectral-line redshift | 0 | 0.0076605 |
| Separation of otherwise identical source events | Unchanged | Target: multiplied by 1.0076605 |

The desired column specifies the requirements of the proposed number-preserving redshift channel; it is not a derived interaction. Both columns can balance the photon-energy ledger. The spectrum distinguishes them.

The same conclusion holds in four synthetic distance examples, including larger extrapolated transfer fractions. Two input spectral lines are propagated explicitly; photon/gravity energy integrals, photon counts and spectral means are tracked separately. These probes are not measured galactic spectra. Counts represent ideal mode occupations, not detected individual gravitons.

## Why a different static void profile does not repair this channel

**Known consequence of time-translation invariance:** any linear stationary converter acts diagonally in temporal frequency,

\[
A_{\gamma,\rm out}(\omega)=t(\omega)A_{\gamma,\rm in}(\omega),\qquad
A_{g,\rm out}(\omega)=r(\omega)A_{\gamma,\rm in}(\omega).
\]

Changing its strength with position, or with a fixed void map, changes the coefficients but does not create new frequencies from an exactly monochromatic input. Frequency-selective conversion can reshape a finite-width line or move its centroid by selectively removing photons; that is not an arbitrary input spectrum uniformly translated to lower frequencies. The constant-mixing example has no such selective reshaping and preserves every line center exactly.

Similarly, if input A(t) becomes A(t-delta), its output is the correspondingly shifted original output. A stationary converter can impose a delay or broaden a pulse, but it cannot universally multiply the separation of arbitrary otherwise identical input events. This timing conclusion is analytic, not a new measured supernova score or a separately simulated full light curve.

In the ideal common-speed limit used here, photon and converted-wave flight adds no relative delay along the same ray. Real dispersion or phase mismatch can change arrival behavior and conversion probability; zero relative delay is not a prediction for every magnetic environment. Ordinary endpoint clocks are assumed in this approximation, not derived from it. Consequently a successful relative-lag check cannot compensate for its missing redshift.

## What a viable conversion mechanism must add

Our desired interaction must change the frequencies of surviving photons. Possibilities to investigate include an inelastic process with an explicit receiving system, a time-dependent coupling, or nonlinear collective dynamics. Each must specify where energy and momentum go, preserve the observed line quality, and derive event timing and clock response from the same process. A static change in mixing strength alone is insufficient.

The earlier quantum receiver and time-dependent field candidates were attempts in that direction, with their documented coherence, resource, clock and momentum gaps. This result does not fill those gaps or exclude every such extension. It narrows the local photon-generated-wave branch: keep ordinary stationary mixing as a possible source contribution, not as the redshift mechanism.

## Verification and goal status

Matrix exponentiation agrees with the analytic solution within 1.7e-16; unitarity error is at most 1.2e-16. Integrated energy balances and a doubled frequency grid pass their stated tolerances. The observed amount of redshift was never used as evidence for this conversion mechanism: only the old calibrated target energy fraction fixed a diagnostic mixing angle.

The redshift calculator remains z=exp(alpha D)-1. This audit excludes one candidate cause under its stated assumptions. Independent physical void mapping, a complete surviving interaction/clock model and a frozen full-model withheld test remain outstanding; all four goal stages stay active.

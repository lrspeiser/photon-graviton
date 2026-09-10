# Theory-level priority review

**The current evidence supports useful empirical redshift and extra-gravity relations. It does not yet establish that photons create the gravitational source.** The largest remaining opportunity is to derive and test the connection between those two effects, together with event timing and energy accounting. More precise stellar fits alone cannot supply that connection.

This review implements the user's request to periodically step back from small technical issues. The new bar-field calculation is saved as a provisional foundation; further numerical refinement is not allowed to substitute for resolving the main physical questions.

## A decisive check of the current connection

The drafted empirical gravity relation is

\[
g_c=Aa_*\left(\frac{g_b}{a_*}\right)^p
   =A a_*^{1-p}g_b^p,\qquad a_*=c^2\alpha.
\]

**Status:** this is the project's empirical power-law template written in known mathematical form. It is not a first-principles derivation or a verified unique formula. The proposed identification of its scale with photon loss does not by itself make the connection predictive.

Because A is adjustable, the transformation

\[
\alpha\rightarrow q\alpha,\qquad A\rightarrow A q^{p-1}
\]

leaves every gravitational prediction unchanged. **This is an exact algebraic consequence of the drafted formula**, not a proposed new law.

We executed the transformation on all **3,150 retained rotation measurements from 149 galaxies**, preserving their original model predictions and sample roles. For q=0.01, 0.1, 1, 10 and 100, the largest numerical change in any predicted rotation speed was 0.000000000000057 km/s. The rotation-only parameter-derivative matrix has rank two for the three quantities log(A), log(alpha), and p, confirming the redundant direction.

In plain language: **we can change the photon-loss rate enormously and compensate with the separately adjustable gravity strength. The rotation fit cannot distinguish those cases.** The redshift observations still constrain alpha; this is not a claim that the combined data cannot measure it. It means that rotation is not an independent prediction of the photon-loss mechanism while A remains free.

For example, the factor-0.01 and factor-100 cases predict radically different redshifts over 100 million light-years, approximately 0.000076 and 1.145, yet give the same rotation curves after adjusting A. The actual redshift data would reject inappropriate values. They do not explain why the gravitational coupling takes the compensating value. See `bridge-identifiability.json` and `bridge_test.py` for the calculation and input hashes.

The local-deposit formulation has a related issue:

\[
\Delta\Phi_c=-\chi\int u_d(\mathbf x')K(|\mathbf x-\mathbf x'|)\,d^3x'.
\]

**Status:** proposed effective response using known kernel/superposition mathematics. Gravity alone constrains the product of response strength chi and deposited-energy distribution, subject to the kernel assumptions. Scaling chi up and the deposited energy down by the same factor gives the same field. An independent transfer/capture/storage/response law is needed to break that freedom. Calling the carriers gravitons does not set that law.

## The five major questions

| Priority question | What the evidence currently says | Why it matters / next decisive outcome |
|---|---|---|
| Does one process stretch both wavelengths and whole events? | Stationary energy loss fails the retained spectral-aging comparison; the time-dependent candidate fixes kinematics but lacks a completed physical source. | Derive a common propagation/interaction rule that predicts spectrum, duration and brightness. Appending a timing factor is not an explanation. |
| Does that process conserve energy and momentum across every participating field? | Photon-plus-companion packet bookkeeping balances; the evolving time field's work and loss-free companion propagation remain unresolved. | Identify the donor, receiver and driver explicitly, including forces and inverse processes. A balanced two-reservoir equation does not account for an omitted driver. |
| Does photon loss independently predict the deposited gravitational field? | The numerical test above exposes an exact adjustable-amplitude degeneracy. Capture shape and storage support also remain unspecified. | Derive a response/capture rule or constrain it with independent observables, then freeze it before predicting new systems. |
| Does one gravitational field fit motions and lensing? | Rotation benchmarks improve substantially; the bulge distribution and lensing have not been jointly predicted. | Use one field equation and the same source, not separate correction factors for stellar motion and light deflection. |
| Does it survive the other observations and genuine predictive tests? | Radio-line agreement is encouraging under limited assumptions; timing, thermal-background and other gates remain unresolved. Existing sample exposure and selection matter. | Maintain the full coverage register and use declared holdouts/external tests. Better performance on one convenient subset cannot certify the theory. |

These conclusions are grounded in the current [electromagnetic audit](../electromagnetic-audit/report.md), [joint galaxy analysis](../joint-galaxy-audit/report.md), [stellar holdout audit](../stellar-holdout-audit/report.md), and the new degeneracy calculation. Their different datasets, assumptions and residual scores must not be pooled into an invented global significance.

## Work order following this review

1. **First resolve the propagation and receiving-sector accounting at the equation level.** Compare explicitly stationary conversion with the drafted evolving time-field candidate. Identify which interaction could actually shift photon frequencies rather than merely remove photons. Require the same rule to predict event stretching and identify where the energy goes. Preserve the user's desired normal local clocks and companion behavior, and state clearly if a candidate cannot meet them.
2. **Derive the connection to gravitational response before promoting it as confirmed.** A common fitted A can remain a phenomenological benchmark, but it must not be presented as evidence for photon origin. Specify what independent observation or field equation fixes chi, the capture fraction, support and spatial deposition law. The user-deferred total source-energy calculation remains deferred, not silently passed.
3. **Then resume the orbital and lensing predictions with declared parameters.** Use the saved bar-field foundation and stellar uncertainties. Finish numerical accuracy only to the level needed by those physical predictions. Keep the reserved stellar scores unopened until the model/selection likelihood is frozen.
4. **Repeat this broad review at each major model change or expensive new analysis.** State the proposed physical claim, the evidence that could refute it, and whether the work reduces a principal uncertainty. Keep failures and parameter degeneracies visible. Do not optimize the pipeline simply because it is easier than the underlying physics.

The goal is a strong, defensible case if the theory earns one. If a shared mechanism fails, report that result and test a clearly labeled revision; do not make agreement inevitable by giving every observable its own adjustable rule. The goal remains active and incomplete.

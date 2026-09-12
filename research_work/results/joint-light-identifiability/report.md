# What makes the joint redshift, timing and brightness test predictive?

This audit derives an exact ambiguity in the declared bolometric model and checks the local DES metadata. It does not fit observations, rule out companion transport, or require a microscopic origin of time before proceeding. It identifies what must be independently constrained so that an apparent agreement has predictive content.

## Conditional model and formula provenance

**Proposed assumptions:** static geometric dilution at physical distance D, photon-number conservation, measured spectral stretch S=1+z, and event stretch S^b. The current shared-stretch candidate has b=1; energy-only transport has b=0. This parameterization does not derive the cause of stretching.

**Statistical nuisance assumptions, not new physics:** allow intrinsic source width and peak bolometric luminosity to have trends

\[
W_{\rm int}=W_0 S^{e_W},\qquad L_{\rm peak}=L_0 S^{e_L}.
\]

These powers describe possible sample/source differences across redshift. They do not assume expansion, a cosmic age, or that such trends actually exist. A frozen hypothesis can set them to zero, provided that assumption is stated and tested rather than hidden.

**Known radiometric bookkeeping, conditional on the assumptions:**

\[
W_{\rm obs}=W_0S^{b+e_W},\qquad
F_{\rm peak}=\frac{L_0}{4\pi D^2}S^{e_L-1-b}.
\]

One energy factor 1/S comes from photon redshift; another factor S^-b comes from the arrival-rate change. For a fixed normalized source shape f,

\[
F(t)=\frac{L_0}{4\pi D^2}S^{e_L-1-b}
 f\!\left(\frac{t}{W_0S^{b+e_W}}\right).
\]

No empirical success follows merely from writing these relations.

## Exact ambiguity

**Derived consequence of this model, using ordinary algebra:** for any constant delta,

\[
b\mapsto b+\delta,\quad e_W\mapsto e_W-\delta,\quad e_L\mapsto e_L+\delta
\]

leaves both width and the entire bolometric light curve unchanged. For example, extra propagation stretching can be offset by intrinsically shorter, brighter sources. This remains true even if every distance is independently known perfectly. Unknown distances introduce further freedom.

The two observable redshift slopes constrain b+e_W and e_L-1-b: two combinations of three parameters. Their Jacobian has rank 2 and null direction (1,-1,1). Adding an independent constraint on either evolution exponent raises the rank to 3. This is conditional identifiability, not proof that a usable independent constraint has been acquired.

## Why total event energy is not a third independent answer

Let C_f be the fixed dimensionless integral of f. The observed bolometric fluence is

\[
\mathcal F=\int F(t)dt
 =\frac{C_f L_0 W_0}{4\pi D^2} S^{e_L+e_W-1}.
\]

The timing exponent b cancels. In this model fluence is the product of peak flux, duration and a shape factor, so its slope adds no independent equation. The transformation above also preserves intrinsic total event energy L_peak W_int C_f. Consequently, an energy-budget anchor alone does not break this particular ambiguity. Additional information about emission duration or peak luminosity can. Real changing spectral/temporal shapes require a richer model; this calculation makes no claim that full spectroscopy has no additional information.

## Executed checks

The script checks six parameter transformations across five redshifts (0.001 to 1.2) and 1001 time samples. Width, peak flux, full artificial bolometric curves and analytic fluence remain unchanged to a maximum relative numerical difference of 1.96e-16. Matrix-rank checks verify the ambiguity and the effect of an additional independent evolution constraint. These are numerical checks of algebra, not astronomical tests.

## What the local DES inputs do and do not supply

The pinned HEAD file has 99 columns, recorded in results.json. It supplies redshifts, positions, host associations, approximate peak dates and host properties. It does not contain an explicit independently measured physical-distance or intrinsic-supernova-luminosity field that closes the joint comparison. Host magnitudes are not supernova absolute luminosities. This statement concerns the audited product, not every potentially available astronomical catalog.

The pinned release README says that AB and Fragilistic offsets are not applied directly to these data; they are defined in the separate `kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.input` calibration file. No filename matching calib/kcor/filter/fragil was found in the acquired timing-input cache. The documented calibration, passband and covariance dependencies must be resolved before using the fluxes for precision physical brightness predictions. Their absence does not invalidate the earlier artificial-flux timing calculations.

Primary input provenance: [pinned DES release](https://github.com/des-science/DES-SN5YR/tree/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/0_DATA). The earlier input audit and this script verify the retained SHA-256 hashes. No fresh flux outcomes or independent distances were acquired in this pass.

## Decision for the joint-light test

1. Retain the user-stipulated galaxy distances as fictional-universe inputs, with their provenance. Do not call a distance inferred from the same supernova brightness or redshift law independent evidence for that law.
2. Acquire and apply the documented passband/zero-point calibration before physical-flux inference.
3. Freeze source-population restrictions, including the allowed duration and luminosity trends, before fitting the joint transport parameters. Local calibration fixes normalizations; by itself it does not determine high-redshift evolution slopes.
4. Constrain at least one relevant source-evolution direction independently, or report the family of equivalent solutions and an explicitly assumption-dependent fit. A fully specified source-emission model is another possible route; it must predict the constraint rather than fit it to this same discrepancy.

This supports a predictive effective theory without demanding ultimate microphysics. It prevents source freedom from being mistaken for evidence that a particular time/companion mechanism is the cause. The active timing calibration continues separately; none of the six scientific demonstrations is complete.

Reproduce from the repository root:

```powershell
python research_work/results/joint-light-identifiability/run.py
```

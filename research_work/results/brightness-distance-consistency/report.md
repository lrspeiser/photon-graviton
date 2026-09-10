# Does model-induced dimming repair the redshift comparison?

2026-09-10. Exposed-data sensitivity connecting stages 2 and 3. No parameters fitted, new data acquired or genuinely withheld observations evaluated. The original stipulated distances remain unchanged in the main analysis. [Protocol](protocol.md), [calculation](run.py), [all predictions](predictions.csv), [results and hashes](results.json).

**Result:** under an ideal bolometric brightness-distance reinterpretation, the shared spectral/event-stretch model lowers median inferred path distance by 0.714%, but all-sample RMS redshift residual changes only from 447.05 to 444.45 km/s. Some individual predictions change substantially more; the large overall scatter persists. This neither establishes a physical correction to published SBF distances nor validates the proposed interaction.

## The conditional calculation

The existing empirical conversion rate is held at alpha=0.0002488993286382367 per Mpc. Under the proposed spectral stretch S=exp(alpha D), event stretch S^b, conserved photon number and static geometric dilution, known radiometric bookkeeping gives bolometric flux attenuation S^-(1+b). For an ideal luminosity-calibrated distance indicator,

\[
D_{\rm app}=D\exp[p\alpha(D-D_a)],\qquad p=(1+b)/2.
\]

D_a is a declared physical calibration-anchor distance. The zero-anchor scenario means unaffected calibration, not an actual source at zero distance. b=1 is the shared spectral/event-stretch candidate; b=0 is the earlier energy-only control that lacks event stretching. Neither case supplies a new matter-clock coupling.

**Known inverse-function mathematics, not a new physical law:**

\[
D=\frac{W(p\alpha D_{\rm app}e^{p\alpha D_a})}{p\alpha}.
\]

For the unmodified baseline p=0, D=D_app. No observed redshift enters this inversion. However, alpha was calibrated from these same previously exposed galaxy data, so the exercise is not independent validation. A joint photometric fit could require a different alpha, and would then need an appropriately frozen new observational test.

## All declared scenarios

Residuals below are c times predicted-minus-observed redshift, expressed in km/s for compatibility with the existing catalog diagnostic. They are not inferred galaxy motions or travel speeds. Every original row is retained.

| Scenario | Anchor | All 164 rows RMS | Historical 104 training rows RMS | Historical 35 validation rows RMS | Historical 25 test rows RMS |
|---|---:|---:|---:|---:|---:|
| Original stipulated distances | — | 447.05 | 457.58 | 437.07 | 415.41 |
| Energy-only brightness sensitivity | 0 Mpc | 445.09 | 457.82 | 431.42 | 408.73 |
| Energy-only brightness sensitivity | 1 Mpc | 445.11 | 457.80 | 431.51 | 408.82 |
| Energy-only brightness sensitivity | 10 Mpc | 445.28 | 457.64 | 432.30 | 409.62 |
| Shared spectral/event stretching | 0 Mpc | 444.45 | 459.31 | 427.11 | 403.51 |
| Shared spectral/event stretching | 1 Mpc | 444.45 | 459.23 | 427.25 | 403.64 |
| Shared spectral/event stretching | 10 Mpc | 444.47 | 458.62 | 428.50 | 404.88 |

All partitions were exposed before this audit. The historical test improvement is not a new blind result; the training RMS worsens in the same scenario. No anchor is selected by these scores, and no significance or full-likelihood preference is claimed.

With shared stretching and unaffected calibration, distance decreases range from 0.253% to 2.242%, median 0.714%. Predicted redshift shifts correspond to -1.93 through -159.54 km/s, median -15.56 km/s. The all-sample bias changes from +14.26 to -18.29 km/s. Thus the correction changes both shape and mean offset; it cannot simply be described as removing a uniform bias. For the 10-Mpc anchor, median distance change is -0.468% and all-sample RMS is 444.47 km/s.

## Uncertainty and verification

**Known derivative algebra conditional on the same model:**

\[
\frac{dD}{d\mu}=\frac{\ln10}{5}\frac{D}{1+p\alpha D}.
\]

The output includes this derivative times each quoted modulus error as a local sensitivity. It omits common calibration, source populations, peculiar motions, endpoint effects and survey selection, so it is not a complete predictive uncertainty and is not used to claim a chi-square success.

Independent bracketed roots match the Lambert W inversion within 2.2e-13 Mpc. Finite differences match the analytic modulus derivative within 4.2e-10 relative error. Forward distance round trips pass, and the p=0 control reproduces the archived baseline scores. These numerical checks do not establish that a real SBF reduction should receive the ideal bolometric correction.

Actual SBF data use specific passbands, colors, population calibration and observational corrections. Existing reductions may already include some spectral transformations; blindly applying another full factor risks double counting. The result is therefore a sensitivity scenario alongside the user's stipulated published distances, not an adopted catalog rewrite.

## Consequences for the goal

Brightness and independent-distance inference must share one explicit calibration model. This check shows that the ideal dimming effect is too small to eliminate the present overall residual scatter, despite changes in particular objects. The void comparison still needs physical environment/selection information; the shared-transport model still needs clock and interaction closure. The local photon-generated gravitational-wave branch must separately derive its spectral and waveform predictions. No full candidate is ready for a genuinely withheld joint test yet. All four goal stages remain incomplete.

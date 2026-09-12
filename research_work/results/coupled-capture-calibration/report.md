# Recalibrating capture with gravitational feedback included

## Outcome

Reducing the capture exposure improves the feedback model's amplitude, but does not recover the earlier no-feedback fit. The recalibrated p=2 andp=3 training RMS discrepancies are23.67 and22.62km/s, respectively. Both curves are too low in the inner bins and too high in the outer bins.

This distinguishes a distribution problem from simply reusing an overly large capture amount. It does not exclude all capture laws, outer boundaries or retention dynamics. The comparison is still a training calibration of a specified100kpc, thin-capture, in-place-retention model.

## What was fitted

**Unchanged project growth postulate:**

\[
\partial_s\rho_D=C\left[\frac{W_{\rm total}}{W_{\rm total}+W_*}\right]^p,
\quad W_{\rm total}=\max[-\Phi_b-\Phi_D(s),0],
\]

with zero initial deposits and exposure s=0..1. Poisson's equation supplies Phi_D and the force from the same density. This is chronological growth rather than assigning a final-density equilibrium. The gravity mathematics is established; the capture, source and permanent-retention prescription is hypothetical.

Keep p=2 or3, Wstar=40000(km/s) squared and the100kpc boundary. Only C changes. It combines incident intensity, capture normalization and duration in consistent units. Recalibrating that combination does not establish that the available photon/companion supply can fund it.

The objective uses the same12 Jeans-proxy bins from542 already processed Cepheid training stars. It averages predicted squared circular speeds at the actual training radii before taking a square root. No additional stellar sample, reserved validation/test stars, vertical inference or lensing observation enters the fit. Existing distance, equilibrium, tracer-profile and ordinary-mass-model assumptions remain. RMS values are unweighted discrepancies, not statistical confidence levels.

## Full curve, not only its best-fitting point

| Model | Training RMS, km/s |
|---|---:|
| Ordinary matter alone | 75.24 |
| p=2, previous no-feedback calibration | 13.78 |
| p=2, feedback at that old amount | 193.15 |
| p=2, recalibrated with feedback | **23.67** |
| p=3, previous no-feedback calibration | 8.76 |
| p=3, feedback at that old amount | 370.33 |
| p=3, recalibrated with feedback | **22.62** |

The selected amplitudes are C=1.16994e7 and1.93629e7 solar masses/kpc cubed for p=2 andp=3. They are several times lower than the previous no-feedback calibrations.

| Training bin | Inferred Jeans proxy | Recalibrated p=2 | Recalibrated p=3 |
|---|---:|---:|---:|
| 6–7kpc | 242.0km/s | 199.7km/s | 201.0km/s |
| 8–9kpc | 235.6km/s | 205.6km/s | 207.0km/s |
| 11–12kpc | 232.6km/s | 216.2km/s | 217.2km/s |
| 17–18kpc | 226.2km/s | 257.3km/s | 255.6km/s |

The whole twelve-bin comparison is in [radial-comparison.csv](radial-comparison.csv). For p=3 the inner deficit remains41.0km/s and the outer excess29.4km/s. An improved overall RMS cannot conceal those systematic residuals. They indicate an overly rising contribution in this calibrated family; individual-star speeds are not equated to these circular-speed proxies.

## Mass and feedback have not become negligible

The fitted models still require deposited effective masses of3.35e13 solar masses for p=2 and2.98e13 for p=3 inside100kpc. Their solar deposited potential depths are about2.18e6 and1.95e6(km/s) squared. These are model outputs, not independently measured Galactic masses. The distant material contributes to binding depth even where its contribution to inner orbital acceleration is small.

Final density divided by C spans0.659–0.880 for p=2 and0.346–0.743 for p=3. Feedback remains substantial. Saturation and permanent retention still reshape the deposited mass distribution; a scalar reduction of incoming exposure cannot independently repair every radius.

No incident-flux budget, binding/support energy or finite-speed formation dynamics was added. The current quasistatic gravity solve and thin-capture assumption remain conditional. A favorable local fit would not by itself settle those requirements.

## Search and numerical verification

The [protocol](protocol.md) fixes a13-point scan in log10 C from4 to8.5, retaining all points. Each interior scan minimum is refined within its neighboring interval. Neither selected minimum is at the imposed boundary, and both optimizers succeed. This is a declared one-dimensional numerical search, not a theorem excluding every unobserved narrow extremum or every parameterized theory.

Calibration uses400 radial nodes,64 angular nodes, even multipoles through24 and128 midpoint exposure steps. At the selected C, separate time refinement to256 steps and combined spatial/time refinement to800 nodes,128 angles, multipoles through48 and256 steps change every predicted bin speed by less than0.105km/s. Both declared0.5km/s gates pass. The quoted final RMS values evaluate the calibrated C on that finer grid, without claiming an independently optimized fine-grid parameter value.

The archived plus/minus1-percent C sensitivities show how exposure changes the predictions; they are not parameter uncertainty estimates. Input hashes and all scan, optimization and refinement results are retained in [results.json](results.json). A partial calibration checkpoint is in [calibration.json](calibration.json).

Reproduce:

```text
python research_work/results/coupled-capture-calibration/run.py
python research_work/results/coupled-capture-calibration/export.py
```

## What needs to change next

Test the physically unresolved spatial assumptions: the outer capture boundary and whether captured energy truly stays where it first lands, or redistributes into supported bound states. Those alternatives can change the mass profile, unlike a single amplitude adjustment. They must have a common rule and an energy/transport accounting; they cannot be independent density corrections chosen for each observed bin.

First quantify boundary dependence of the coupled calibration. Then use an explicit retention/redistribution model if the fixed-location assumption is changed. Keep the full radial predictions, candidate lensing and vertical-force consequences, and required supply visible. Persistent redshift, admissible vertical observations and genuinely unseen validation remain unresolved. The astronomical alpha and reserved-data roles are unchanged; all nine goals remain active.

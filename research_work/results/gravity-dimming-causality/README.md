# JR-3: can low light be caused by conversion into gravity?

**21 September 2026, America/Los_Angeles. Executed exploratory diagnostic, not a causal detection.**

## Theory first

Matter and radiation may exchange energy with a companion state that contributes to gravity. Low observed light might therefore be an effect of exchange, rather than evidence of a smaller source. This experiment tests that interpretation using the archived 149 SPARC galaxies, 3,152 rotation measurements, six SLACS lenses and 40 stellar-motion bins. No clusters were analyzed.

The exact R10 force/source law and universal parameters stay fixed. Missing light is propagated into both ordinary gravity and the source-linked companion, without adding the same companion energy twice. R10 remains a phenomenological implementation, not the overall theory. Failure here does not refute other source, memory, spectral or formation mechanisms.

Baseline main: `41ede5342a0ddabdff7dbea2a4f1951eceec65c6`. Primary protocol: `2841a36a9d443cb4ca6d513056e3def5b9d1ca08`, committed before the primary run. AMENDMENT-1 was written locally after primary results and before its reversible follow-up; it was published afterward and is not retroactively called preregistered. All computations ran locally, not in GitHub Actions.

## 1. Hidden-light calculation

Define band transmission T by L_obs=T*L_intrinsic, u=1/T and net dimming f=1-T. At fixed intrinsic stellar mass-to-light ratio, a 50% light loss requires twice the stellar source, not 50% more. Distances, inclinations, light-profile shapes, velocities, error bars and lens geometry remain unchanged.

This tests uniform post-emission band transfer at fixed stellar populations. It does not simulate gravity suppressing star formation, changing the luminosity of each star, redistributing light spatially, or creating a frequency-dependent spectrum. Missing light and a larger stellar mass-to-light ratio are algebraically degenerate in this test. The old R10 coefficients were themselves inferred using the original photometry, so the result is not an unbiased opacity limit for the entire companion hypothesis.

| Calculation | Curves below 20% fractional RMS | Median fractional RMS | Mean per-galaxy RMSE |
|---|---:|---:|---:|
| Original R10 | 95/149 | 15.05% | 17.71 km/s |
| All galaxies lose 25% of starlight | 83/149 | 17.21% | 20.47 km/s |
| All galaxies lose 50% of starlight | 57/149 | 25.33% | 31.91 km/s |
| Individually optimized 0–90% missing starlight | 107/149 | 10.87% | 14.50 km/s |
| Individually optimized photon gain or loss | 117/149 | 9.59% | 11.16 km/s |

The individual fits minimize fractional RMS; original-error chi-square is also recorded, and a separately labeled chi-square optimization is preserved. No fitted error floor. These free-transmission rows are inverse feasibility calculations, not a universal causal law or independently measured opacities.

The dim-only fit assigns nonzero loss to 62 galaxies and zero to 87. The signed diagnostic gives 63 net losses and 86 gains; 34 reach the allowed tenfold-brightening bound. Those gains require a supplied incoming energy channel or a different population, not energy from nothing.

### The difficult subset separates by sign

Of the 54 original greater-than-20% outliers, 38 have overpredicted motion and 16 underpredicted motion. **Dimming repairs 12 of the 16 underpredicted cases, but none of the 38 overpredicted cases.** All 38 prefer zero loss. A scan of 101 source corrections from u=1 to u=10 gives no decreasing rotation step across the 3,152 measured radii. This is a measured property of the fixed implementation, not a theorem about every gravity mechanism.

The underpredicted subset's median error falls from 27.05% to 12.69%, with median missing-light fraction 68.1%. The twelve repaired cases require approximately 60–88% hidden starlight. UGC07399 changes from 32.6% to 8.6% RMS with 81.7% loss; UGC07125 remains at 55.8% and prefers zero loss.

JR-2's low-surface-brightness outliers are gas-rich relative to inferred stars, not gas-poor. A grey attenuation sensitivity applied equally to stellar and gas flux gives 106/149 curves below 20% after individual dim-only fits; universal 25% and 50% losses give 81/149 and 44/149. Infrared and 21-cm conversion rates are not thereby derived to be equal. Reducing line brightness would not itself destroy gas atoms.

Two shared laws were fitted to the 89 training galaxies: constant log u=tau, and log u=tau*H(u), where H is the companion fraction of modeled squared circular speed at Re, solved self-consistently on the hidden-light-corrected source. H does not directly use an observed rotation speed, but inherits the previously fitted model. **Both choose tau=0.** The existing 29/31 comparison groups then retain their original predictions. All these samples were previously exposed; this is not fresh blind evidence.

## 2. Conversion back into light

The follow-up solves an explicit equal-rate reversible channel. With incident photon energy 1, companion energy X and total exchange depth a:

    T = 1 + (X-1)*(1-exp(-a))/2
    companion_out = 1+X-T.

Both output energies sum to the input. A return signal needs incoming companion energy. The chosen X=0 equal-rate incoherent limit can remove at most 50%; that is not a bound on coherent or unequal-rate conversion.

For one shared model, set X=R*H(u), a=tau*H(u), and solve u*T(u)=1. Only R and tau are fitted, with six declared starts. The influx prescription is phenomenological, not a measured or derived transport law.

The best executed training solution has tau=1.04943 and **R=9, at its upper bound**. It prefers **net brightening in 148 galaxies and dimming in one**, with median photon multiplier 2.53. Training mean fractional-square loss decreases from 0.13705 to 0.10927, but comparison performance worsens: 29/31-group mean RMSE changes from 17.68/15.09 to **22.66/17.39 km/s**. Across all 149, only 76 remain below 20%, with median fractional RMS 19.63%.

The fit trades large outliers against other objects; its improved training objective is not general success. This candidate is not promoted. Its return-energy supply and replenishment are not established. Matrix-exponential comparison error is 1.78e-15, self-consistency error 5.33e-15, and channel-energy error 8.88e-16. Each final galaxy solution has one bracketed root. Numerical success does not validate the influx hypothesis.

## 3. Lenses

At fixed stellar populations, the missing-light fractions needed to match the observed ring are:

| Lens | Missing light | Original stellar RMS | Stellar RMS after correction | After also adjusting orbital anisotropy |
|---|---:|---:|---:|---:|
| J1112+0826 | 17.00% | 1.69% | 6.95% | 4.75% |
| J1621+3931 | 21.31% | 3.50% | 6.05% | 5.19% |
| J1630+4520 | 20.20% | 1.82% | 8.17% | 6.40% |

The other three rings require about 5% brightening rather than dimming. Rings are matched by construction, but the same mass correction worsens the independently scored stellar profiles. The inherited conditional distances, spherical deprojection, constant population scaling and restricted orbital family remain assumptions. No new spectra, full lens-image likelihood or causal light-loss measurement was obtained.

## 4. Stored energy is not the same as present energy loss

A separate reversible energy model uses radiation E_gamma, companion energy E_chi, input luminosity L_in and nonnegative rates:

    E_gamma_dot = L_in - k_escape*E_gamma - k_plus*E_gamma + k_minus*E_chi
    E_chi_dot = k_plus*E_gamma - (k_minus+k_chi_escape)*E_chi
    L_out = k_escape*E_gamma.

Therefore

    d(E_gamma+E_chi)/dt = L_in-L_out-k_chi_escape*E_chi.

At steady state, L_in-L_out=k_chi_escape*E_chi. **If companion energy cannot escape and can return, a nonzero stationary stored state coexists with L_out=L_in.** A large stored gravitational state does not automatically imply continual dimming today.

Persistent dimming instead requires a growing state, energy escaping through another channel, or energy returning outside the observed band/direction/aperture. A different source-formation history is another mechanism, not calculated here. This conservation statement concerns all-direction bolometric power; it does not equate measured 3.6-micron brightness with total power.

Analytic stationary solutions agree with independent linear solves across 81 rate fixtures to 1.42e-15 relative error. Four initially empty evolutions were also compared to exact exponentiation; the largest recorded discrepancy is 3.13e-9. The outgoing-energy ledger is included. These are dimensionless controls, not observed galaxy rates or ages.

## 5. Conditional energy supply

Only if the original R10 reservoir gravitates as ordinary stored energy, M_chi=A*r_t/G requires E=M_chi*c^2. Let Lcat be dimensionless catalog band luminosity and DEFINE unknown B by L_bol,obs=B*Lcat*L_sun,bol. For a grey, fixed-history example, loss is B*(u-1)*Lcat*L_sun,bol. With retention eta and effective history T:

    B*eta*T*(u-1) >= M_chi*c^2/(Lcat*L_sun,bol).

The right side has median 2.33e14 years for the extrapolated envelope, or 5.48e13 years for mass within the last measured radius. These are energy/power units, **not inferred ages**; B is unknown, and no universe age is imposed. Returning the same energy cannot count as new accumulated inventory on every cycle. This is a lower supply target for the existing reservoir, not an added gravitational component or a complete budget for larger hidden-source corrections.

A different energy-to-gravity response, incident energy from elsewhere, or another historical source needs an explicit relation and accounting. This conditional calculation does not refute those branches. One infrared band is not silently converted into measured bolometric watts.

## Conclusion and reproduction

There is a compatible minority, but these data do not establish dimness caused by companion conversion. The strongest unresolved physical distinction is growing storage versus light returned elsewhere versus genuinely altered stellar production. A spatial/spectral source-and-return law must determine which happens; assigning an opacity to each galaxy cannot establish the cause.

The original observations, uncertainties and formula hashes remain unchanged. Disk replay error is 5.68e-14 km/s. Independent arithmetic readback gives **6,201 checks, zero failures**. The complete originating-conversation bundle retains original inputs, all positive/negative cases, exact executed code, full arrays, optimizer starts, logs and the longer report; the compact results are not a replacement for the full arrays.

Use the extracted JR-2 package (containing original-JR1 and code), with a fresh output directory:

    python test_dimming.py --jr2 Photon-Graviton-JR2 --output results-replay
    python reversible_followup.py --jr2 Photon-Graviton-JR2 --output results-replay
    python readback_audit.py --jr2 Photon-Graviton-JR2 --output results-replay

Primary catalog provenance: Lelli, McGaugh & Schombert (2016), SPARC, https://arxiv.org/abs/1606.09251. Lens, population and TDCOSMO inputs retain the exact pinned sources and assumptions of JR-1/JR-2. No conventional dark-matter or MOND substitute was used for the tested companion branch.

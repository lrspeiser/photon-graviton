# Stronger ordinary wells do not improve frozen capture transfer

The previous Coma turn produced a limited shear-shape test, not a complete cluster solution. This run returns to actual galaxy rotation and tests a physically motivated modification of the interception model: a stronger ordinary gravitational field increases capture probability. No individual galaxy exposure is fitted.

## Formula and conservation

The additional **postulate**, using familiar power-law mathematics, is

    F_g=[g_b(R_disk)/g_ref]^q,
    k0,g=k0 F_g, C_g=C F_g,
    g_ref=1000 (km/s)^2/kpc.

The same F_g multiplies opacity and the density source normalization. Consequently C_g/k0,g=C/k0 for every galaxy: incoming integrated companion exposure stays common. Changing opacity affects both capture and attenuation. Omitting the factor in C_g would silently change the assumed external exposure from galaxy to galaxy. g_ref is a normalization convention, not a derived fundamental acceleration.

The local radial opacity shape and all-direction transport integration remain the previous model. g_b at the disk scale is constructed from the retained gas and fixed stellar mass/light templates; interpolation and inner/outer continuations are stated in the protocol. This is an ordinary-force proxy, not total binding depth or self-consistent feedback. No microscopic interaction or originality claim is established.

## Fit and transfer

Fitted on 89 training galaxies: q=0.625051, C=1.09272e8 Msun/kpc^3, k0=0.544596 kpc^-1 and a/R_disk=2.208277. All three starts converge away from bounds. Parameters are frozen before prediction on the other groups. Existing partitions have been exposed repeatedly, so this is exploratory transfer rather than blind evidence.

| Group | Previous interception RMS km/s | Well-response RMS km/s | Previous log10 RMS | Well-response log10 RMS |
|---|---:|---:|---:|---:|
| Training, 89 galaxies | 30.709 | 30.832 | 0.141004 | 0.135328 |
| Validation, 29 galaxies | 33.184 | 35.941 | 0.119980 | 0.120438 |
| Test, 31 galaxies | 23.994 | 26.622 | 0.098545 | 0.099329 |

Values use finer evaluations at fixed fitted parameters. The fitting objective is mean squared log10 speed residual with equal weight per galaxy. It improves on training, although absolute km/s error does not. Both evaluation groups worsen on both metrics. Thus the preferred positive q in training does not justify adopting the modification.

This result narrows the model: the particular global capture prefactor based on ordinary acceleration at one disk scale does not give better transfer. It does not rule out response to full potential, tidal structure, stellar populations or actual spatial illumination. Those would require different specified equations and predictions, not retrospective explanations of these residuals.

The previous q=0 candidate remains the comparison model; its predictions were reproduced within 1e-7 km/s by the updated implementation. The new fit's resolution change in group RMS is below 0.00098 km/s. The old result files remain untouched; run.py now accepts `--well-strength` and writes separate outputs. Source supply, retention, lensing, timing/brightness and the other requirements remain unresolved. All six goals stay open.

Reproduce: `python research_work/results/isotropic-galaxy-transfer/run.py --well-strength`. well-strength-results.json records parameters, starts and all group scores; well-strength-predictions.json records every measured/predicted rotation point. The input hashes and earlier source/geometry limitations remain in force.

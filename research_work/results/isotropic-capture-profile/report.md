# Isotropic companion arrival: where the deposits form

## Result

Strong capture can concentrate deposits in outer layers even when the capture coefficient is larger deeper inside. The reason is attenuation: many companions are absorbed before reaching the center. This supports the proposed interception mechanism conditionally, without inserting a desired dark-halo profile.

It also exposes a tradeoff. A spherical outer shell contributes no Newtonian force inside itself, although it contributes to the potential and to projected lensing mass. The strongest edge concentration therefore produces very little extra circular speed deep inside this model. Edge loading alone is not a successful stellar-motion explanation.

The spatial profile is now derived from a specified incoming field and capture rule. The incoming companion supply, retention mechanism, disk/bar geometry and microscopic cross section are still assumptions. No astronomical data were fitted.

## Capture rule and derived profile

**Project postulates, not claimed unique:** a sphere of radius10a receives isotropic companion specific intensity I0 at its boundary. Rays travel straight at c0 and are absorbed at coefficient

\[
\kappa(r)=\frac{k_0/a}{1+\sqrt{1+(r/a)^2}},\quad r\leq10a.
\]

Outside the test sphere kappa is zero. This is the saturating binding-depth rule specialized to a fixed baryonic Plummer potential. It is not the deposited density profile and is not inferred from a dark-matter fit. Its finite cutoff is a model boundary, not an observed galaxy edge. Weak focusing, scattering, re-emission and changing opacity from deposit gravity are omitted.

**Known absorption mathematics, applied conditionally:** at radius r, for a propagation direction with radial cosine mu, the backward distance to the boundary is

\[
s_b=r\mu+\sqrt{R^2-r^2+r^2\mu^2}.
\]

Integrating kappa backward along that ray gives optical depth tau(r,mu). The incoming intensity at the point is I=I0 exp(-tau). The deposited power per volume is therefore

\[
\boxed{q(r)=2\pi\kappa(r)I_0\int_{-1}^{1}e^{-\tau(r,\mu)}d\mu.}
\]

This includes arrival from every direction. **Retention postulate:** energy remains where captured, without escape, giving effective gravitational density rho_D=q*t_age/c0 squared. Orbital support or a stable bound field configuration is not demonstrated by that rule.

## How much reaches the outside versus inside?

| Capture coefficient k0 | Fraction of incident power captured | Fraction of deposits outside radius5a | Local deposition q(9a)/q(a) |
|---|---:|---:|---:|
| 0.01 | 1.65% | 78.78% | 0.243 |
| 0.1 | 15.05% | 79.73% | 0.272 |
| 1 | 73.80% | 88.24% | 0.991 |
| 10 | 99.40% | 99.975% | about1.62 million |

The outer half of the radius contains **87.5% of a sphere's volume**. Consequently, the first two large outer mass fractions do not establish preferential deposition per unit volume: their local deposition is still larger inward. The strongest case does show a pronounced outer layer and a heavily deprived center. The enormous local ratio reflects the nearly empty center, not unlimited energy generation.

For comparison, the known Plummer density used in the earlier diagnostic is much more centrally concentrated. The present result does not reproduce that shape merely because the baryonic well setting kappa has a Plummer potential. Nor does the roughly inverse-radius thin-capture density automatically produce a flat rotation curve.

## Force, lensing and clock consequences from that same profile

**Conditional ordinary weak-field gravity:** enclosed effective mass gives g_r=-G M_enclosed/r squared and v_c squared=G M_enclosed/r. The potential includes both interior and exterior spherical shells. Projected mass gives the equal-potential metric's static Born deflection 4G M_projected(<b)/(c0 squared b). These are established formulas, not new project gravity laws.

To compare shapes, choose the intensity for each case so its initial central potential depth is |Phi(0)|/c0 squared=1e-6 after a retention age40a/c0. This is an explicit normalization, not a prediction of the available cosmic companion intensity. The required intensities and all profiles are archived. Holding a common potential depth does not hold the total mass or incoming intensity fixed.

| k0 | Extra v_c at a | Extra v_c at5a | Extra v_c at9a | Static deflection at impact5a |
|---|---:|---:|---:|---:|
| 0.01 | 42.37 km/s | 147.43 km/s | 212.60 km/s | 0.40689 arcsec |
| 0.1 | 40.85 km/s | 145.27 km/s | 213.38 km/s | 0.40511 arcsec |
| 1 | 25.47 km/s | 119.48 km/s | 218.24 km/s | 0.37830 arcsec |
| 10 | 0.018 km/s | 6.41 km/s | 139.49 km/s | 0.24325 arcsec |

These are deposit-only spherical predictions, not observed stellar velocities or full galaxy models. With other spherical components, squared circular-speed contributions add. Bulge/bar stars require orbit distributions; one cannot compare every star directly to a circular-speed formula. A disk's geometry also does not obey the spherical shell cancellation used here.

Constant retained deposition gives a specified growth history from the same q: rho_D(r,t)=q(r)*(t_age+t)/c0 squared. The previously tested equal-potential diagnostic metric then produces central-path redshifts from9.32e-7 to1.17e-6 for rays from-20a to20a. The depth grows by approximately a factor of two during these illustrative crossings. Clock/event checks include endpoints and pass. The metric remains a prescribed weak-field response, not a full coupled Einstein/transport solution; the temporal receiving sector's own gravity is not included in this spherical diagnostic.

## Conservation and numerical evidence

The incoming power is 4pi squared R squared I0. An independent projected-chord absorption calculation agrees with the volume integral of q to relative2.11e-4 or better on the fine grid. Capture fractions lie between zero and one. This closes incident/transmitted/absorbed companion energy for this transfer problem; it does not demonstrate that photon conversion upstream supplies I0.

Two resolutions use321/641 radial nodes and48/96 quadrature nodes per angle and ray integral. All enclosed-fraction refinement differences pass the declared0.003 gate. The first projected-mass quadrature raised roundoff warnings at interpolation kinks; segmentwise integration with a smooth transformed variable fixes them. Independent8/16-node segment calculations pass relative1e-9. Source and receiver proper-clock intervals agree with the ray spectral factor within5.2e-11.

The steady illuminating field is assumed established before the retained-deposit interval. Maximum companion transit is20a/c0, half the chosen age; startup cannot simply be declared negligible. No full switch-on or self-consistent growth simulation is claimed. Binding/support energy and any gravitational radiation during capture still need an accounting.

[Protocol](protocol.md), [solver](run.py), [profiles](profiles.csv) and [results](results.json) are reproducible with:

```text
python research_work/results/isotropic-capture-profile/run.py
```

## Next discriminating test

Determine whether a single binding-depth capture rule can yield an extended profile with the required radial and vertical forces, without excessive edge shielding or center depletion. Vary the physical depth dependence as a declared common law, rather than fitting a density independently at each radius. Then carry the surviving law into a nonspherical ordinary-matter geometry and compare motions and lensing with uncertainties. The intensity and formation history must ultimately come from the photon/companion budget. No current astronomical fit, prior alpha, or unseen-data status changed; all nine goals remain incomplete.

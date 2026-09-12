# Deposits linked to gravity, clocks and light bending

## Result and scope

A prescribed growing attractive well can redshift a crossing signal even after source and observer clock factors are included. A static well with symmetric equal-depth endpoints gives zero net spectral/event stretch. The growth result is a known type of evolving-potential effect, not a newly discovered redshift mechanism: see [Laguna et al., gravitational-wave Sachs-Wolfe effects](https://arxiv.org/abs/0905.1908). That paper's expanding background is not adopted here.

The same potential now supplies explicit stellar acceleration and lensing predictions. This advances goal5 at a conditional level. It does not derive the deposited spatial profile, its growth rate, or a sufficient radiation supply. Those are prescribed inputs in this experiment, not outputs of the preceding one-dimensional capture simulations.

## Deposit-to-potential rule

**Conditional use of known Newtonian gravity:** assume deposited energy acts as a nonrelativistic effective source rho_D=u_D/c0 squared, with Poisson equation Laplacian Phi_D=4pi G rho_D. This is the ordinary coupling baseline, not an arbitrary enhancement of gravity per unit energy. Pressure, momentum and field energy need a fuller treatment.

Choose the known positive Plummer profile as a diagnostic:

\[
\rho_D(r,t)=\frac{3m(t)a^2}{4\pi(r^2+a^2)^{5/2}},\qquad
\Phi_D(r,t)=-\frac{Gm(t)}{\sqrt{r^2+a^2}}.
\]

**Conditional derivatives, not new force laws:**

\[
\mathbf g_D=-\frac{Gm(t)\mathbf r}{(r^2+a^2)^{3/2}},\qquad
v_c^2(r)=\frac{Gm(t)r^2}{(r^2+a^2)^{3/2}}.
\]

For cylindrical radius R and height Z, substitute r squared=R squared+Z squared. The vertical component is g_Z=-Gm Z/(R squared+Z squared+a squared)^(3/2), and the radial component replaces Z by R. Thus a locally deeper well has directional force consequences through its gradient; no separate choice of a pulling direction is needed. This spherical diagnostic is not a bulge/bar/disk fit or proof that companions form that profile.

## One clock and lensing prescription

**Proposed diagnostic metric using established metric mathematics:**

\[
ds^2=-e^{2\phi}c_0^2dt^2+e^{-2\gamma\phi}d\mathbf x^2,
\qquad \phi=\Phi_D/c_0^2.
\]

Gamma here is written as lowercase gamma and is a spatial-curvature parameter, not the capture rate used in earlier experiments. Gamma=0 retains the prior fixed spatial rods. Gamma=1 gives equal weak-field time and spatial potentials, as in the usual weak-field gravity baseline. The exponential completion is not claimed to solve the full Einstein equations for growing Plummer matter. For gamma=1 local spatial geometry changes as the well grows; there is no homogeneous expansion factor. This is an explicit departure from strictly fixed spatial rods.

Clocks held at fixed positions read d tau=e^phi dt; local physical lengths are e^(-gamma phi)|dx|. Shared light/GW coordinate speed is c0 exp[(1+gamma)phi], and their locally measured speed is c0. The same messenger law is imposed, not observationally inferred in this calculation.

**Known weak-field lensing relation applied conditionally:** a static snapshot gives Born deflection magnitude

\[
\widehat\alpha(b)=\frac{2(1+\gamma)Gm b}{c_0^2(b^2+a^2)}.
\]

At fixed dynamical potential, gamma=0 yields half the gamma=1 bending. This is why a time-only force prescription cannot silently inherit the lensing of a model with equal spatial and temporal potentials. Parametrized metric tests of clock rates and light bending are established methodology; see [Will's review](https://arxiv.org/abs/1403.7377). No observational bound from that review is fitted here.

## Growing-well propagation test

Use units a=c0=G=1. The specified mass history is m(t)=beta[1+f t/20], with beta=1e-6 or1e-5, f=0,0.1 or1, and gamma=0 or1. These are twelve histories. A central radial signal travels from x=-10 to10 starting at time zero. The well approximately doubles its mass during the crossing when f=1; no claim is made that real wells grow this quickly.

Let n_opt=exp[-(1+gamma)phi]. **Conditional ray and clock derivation:**

\[
\frac{dt}{dx}=n_{\rm opt},\quad
\frac{d\ln J}{dx}=\partial_t n_{\rm opt},\quad
S=1+z=J e^{\phi_o-\phi_e}.
\]

In the weak-field limit for this path,

\[
\ln S\simeq\beta f\left[(1+\gamma)\frac{\operatorname{asinh}(10)}{10}-\frac{1}{\sqrt{101}}\right].
\]

The first term is propagation through the evolving potential; the second is the opposing receiver-clock change. Their sum is positive for the tested growing cases. This is a conditional derivation from the imposed mass history, not a universal distance-redshift law or a microscopic photon-conversion rate.

| Initial depth beta | Initial v_c at r=a | z for mass doubling, gamma=0 | z for mass doubling, gamma=1 | Snapshot deflection at b=a, gamma=1 |
|---|---:|---:|---:|---:|
| 1e-6 | 178.26 km/s | 2.00319e-7 | 5.00142e-7 | 0.41253 arcsec |
| 1e-5 | 563.70 km/s | 2.00321e-6 | 5.00150e-6 | 4.12530 arcsec |

These speeds are illustrative predictions from the selected depths, not measured stars. Static controls give z=0, while ten-percent-growth cases give approximately one tenth of the doubling shifts. Lensing is a separate static snapshot at initial mass; it is not a full time-dependent off-axis lens calculation. The circular speeds are instantaneous weak-field diagnostics, not solutions of stellar orbits in the growing potential.

The magnitude is tied to the potential depth and its growth. Increasing it by changing the mass also changes stellar forces and lensing. An independently adjustable amplification of clock change would be an additional postulate requiring its own tests; it has not been inserted to match the old alpha calibration.

## Energy consequence

An independent integral of one half rho_D Phi_D gives the known gravitational potential self-energy

\[
U=-\frac{3\pi Gm^2}{32a},\qquad
\dot U=-\frac{3\pi Gm\dot m}{16a}
\]

at fixed a. The increasingly negative potential energy must be included with deposited rest energy, kinetic/support energy and any outgoing radiation. U alone is not the full binding energy of a virialized system: supporting motions also contribute. The archived numeric field named binding_energy denotes this gravitational potential term only. Nothing here closes the formation-energy ledger, justifies permanent support, or proves the prescribed m(t) can be funded by the earlier transport model.

## Verification and next step

[Protocol](protocol.md), [code](run.py) and [results](results.json) are archived. Run `python research_work/results/deposit-clock-gravity/run.py`. Twelve ray calculations agree with independent finite-event proper-clock integrations within 1.6e-11. Independent mass, deflection and potential-energy quadratures meet their declared 1e-9 relative gates. The weak-field approximation differs from the exact diagnostic-metric log stretch by at most 7.94e-11 in the tested range.

This identifies a sign-consistent deposit-feedback route, with mandatory force/lensing consequences. Next derive the profile and mass-growth history from captured energy, then compare those coupled predictions with stellar and lensing data. Include the existing temporal field's own gravity rather than double-counting deposits or ignoring its forces. Real-data fits, a complete stress-energy theory, cosmological coverage and a genuinely unseen test remain unfinished. All nine goals stay active.

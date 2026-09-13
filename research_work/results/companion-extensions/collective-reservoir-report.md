# Reciprocal transfer from a collective domain to finite protected storage

13 September 2026. Finite-state stationary model with explicit reverse export.

## Result

Collective enhancement can help replenish a leaking protected level, but does not independently set its equilibrium occupancy. With zero protected decay, the occupancy derived below is independent of collective domain size. A returning relaxation field can undo protection. The model therefore separates faster transfer from greater storage capacity rather than identifying them.

## Hypothesis and rates

Use the previous symmetric domain with m=0,...,N excitations, each of energy delta. Add one protected level p=0 or 1 of energy delta/2. The relaxation quantum also has energy delta/2, conserving energy in each export. This changes the earlier instantaneous packet-reset architecture: here an export removes one excitation and fills one finite protected level.

Known collective matrix factors are U_m=(N-m)(m+1) and D_m=m(N-m+1). Adopt pump transitions

\[
(m,p)\to(m+1,p): A U_m,\qquad
(m,p)\to(m-1,p): B D_m.
\]

A reciprocal interaction proportional to S_- b^dagger c^dagger plus its Hermitian conjugate motivates the export rates (b denotes the two-state protected raising operator, c the relaxation mode):

\[
(m,0)\to(m-1,1):\kappa(n+1)D_m,
\]
\[
(m-1,1)\to(m,0):\kappa n U_{m-1}.
\]

Because U_(m-1)=D_m, both directions have the same collective matrix factor across a link. Bosonic stimulation and symmetric-spin algebra are known; this effective coupling and choice of protected state are our hypotheses. The relaxation occupation n is held fixed externally. Spatial propagation, bath heating, recoil and dephasing are not included.

A protected decay transition (m,1) to (m,0) at rate gamma releases delta/2 into a separate outgoing channel, assumed not to return. Gamma is stipulated, not derived. We use A=1, B=r and kappa=0.01, so time is measured in inverse-A units and delta sets energy units.

## Exact no-decay equilibrium

With gamma=0, detailed balance within the symmetric manifold yields

\[
P(m,p)\propto r^{-m}w^p,\quad
w=\frac{n+1}{rn},\quad r=B/A,
\]
\[
\boxed{P(p=1)=\frac{n+1}{n+1+rn}.}
\]

For n=0 the limiting occupancy is one. This is a derived equilibrium of our finite model, not a universal retention law. N and kappa cancel from the occupancy, although they affect kinetics. At r=0.9 and n=10, protection is 55% for either tested domain size. Net stationary transfer is zero when there is no protected decay: a finite level fills rather than storing indefinitely.

## Finite leakage and energy accounting

Let J denote forward exports minus inverse exports. Stationarity implies

\[
J=\gamma P(p=1),\qquad
P_{input}=\delta J,
\]
\[
P_{relax}=\delta J/2,\qquad
P_{decay}=\delta\gamma P(p=1)/2.
\]

Thus net source input equals net relaxation output plus decay output. The nonzero current with leakage is replenishment, not growing stored energy. Stored energy includes both delta times mean m in the active domain and delta P(p=1)/2 in the protected level; the active excitation inventory cannot be ignored.

For r=0.9 and gamma=0.1:

| Domain size | Relaxation occupation n | Protected occupancy | Net transfer rate in A units |
|---:|---:|---:|---:|
| 10 | 0 | 0.671940 | 0.067194 |
| 100 | 0 | 0.984972 | 0.098497 |
| 10 | 10 | 0.534859 | 0.053486 |
| 100 | 10 | 0.549286 | 0.054929 |

The larger domain approaches its no-decay occupancy more closely. However, it contains ten times as many constituents while feeding the same one-level reservoir, so this table does not demonstrate improved storage per unit material. At r=0.9 the active domains contain about 6 versus 91 excitations in the no-decay stationary state, all requiring source energy. A realistic comparison needs reservoir multiplicity and domain abundance fixed consistently.

## What this contributes

We have removed the assumption of one-way export from this finite model and obtained a useful analytic occupancy law. It shows why coherence alone cannot guarantee permanent accumulation: returning waves and reservoir lifetime remain decisive. A zero-decay, empty-returning-bath limit fills the level, but the physical protection and empty-bath assumptions still need justification.

Next, any proposed reservoir should specify its capacity, leakage spectrum and occupation of returning modes. Those quantities must connect to the source/outlet and galaxy inventory calculations. This finite local model does not supply a traveling-companion mechanism, a halo distribution or a new redshift prediction. Its half-energy split is an explicit diagnostic choice and does not replace the exact one-third empirical reference.

## Verification

Fifty-four stationary generators cover N=10 or 100; r=0.9, 1 or 1.030301; n=0, 0.1 or 10; and gamma=0, 0.001 or 0.1. All pass probability, generator residual, net-current and energy-ledger checks with absolute tolerance 1e-8 in the stated units. Eighteen no-decay cases agree with the full analytic product distribution, not just its protected marginal. Tiny signed currents around numerical zero remain visible in the JSON as solver residuals, not physical negative emission.

Files: `collective-reservoir.py`, `collective-reservoir-results.json`. Related: [collective normalization](collective-normalization-report.md), [finite bath](finite-bath-report.md).

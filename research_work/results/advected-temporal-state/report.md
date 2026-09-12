# Traveling temporal energy: propagation and capture audit

## Outcome in ordinary language

We tested a direct way to avoid leaving converted energy behind in voids: let the traveling companion energy itself determine the local time-field state. It balances the stated reference-energy ledger, but its feedback changes how disturbances propagate. In one example, energy flows at the intended light speed while a disturbance in the coupled field propagates at 1.724 times that speed. More importantly, capture drives initially regular examples into a boundary with unbounded high-frequency amplification in the linearized equations.

This particular instantaneous coupling is not ready to serve as the complete physical rule. It is a conditional failure of specified equations, not a rejection of all nonexpanding or companion models. No astronomical observations were fitted in this audit.

## Postulates and meanings

All equations below use nondimensional variables with reference speed c0=1. R is the sum of electromagnetic and gravitational-wave reference-energy densities, C is traveling receiving-sector reference-energy density, and D is deposited reference-energy density. K is a positive energy-density scale. Gamma is a nonnegative capture rate per reference time. These are not yet a covariant stress-energy or momentum completion.

**Project postulates, without a uniqueness claim:**

\[
n=1+C/K,\qquad v=1/n,\qquad d\tau=dt/n,
\]

\[
R_t+(vR)_x=-\frac{n_t}{n}R,
\quad C_t+(\eta vC)_x=\frac{n_t}{n}R-\Gamma C,
\quad D_t=\Gamma C.
\]

The primary eta=1 choice carries receiving energy at local c0. Eta=0.5 is a separate sensitivity, not an adopted change to the user's preferred propagation speed. The shared ray law remains omega=c0|k|/n. Measured redshift must still include source and receiver clock factors as derived in the [clock completion](../inhomogeneous-clock-completion/report.md).

**Conditional derivation using ordinary conservation mathematics:**

\[
\partial_t(R+C+D)+\partial_x[vR+\eta vC]=0.
\]

There is no additional stationary memory energy in this closure. Nevertheless, conserving this reference sum does not alone prove conservation of physically measured energy and momentum. The state-to-energy and clock rules must ultimately be given a consistent operational completion.

## Propagation calculation

**Conditional algebraic consequence of the postulates:** define d=1-R/(Kn). Substitution eliminates n_t and gives

\[
d C_t+\frac{\eta}{n^2}C_x=-\Gamma C.
\]

The source-free local principal matrices for y=(R,C) are

\[
M=\begin{pmatrix}1&R/(Kn)\\0&d\end{pmatrix},\qquad
B=\begin{pmatrix}1/n&-R/(Kn^2)\\0&\eta/n^2\end{pmatrix}.
\]

**Known characteristic analysis applied to this proposed system:** the eigenvalues of M inverse times B are 1/n and eta/(n squared d). In the stipulated local clock units, the receiving-state characteristic speed relative to c0 is

\[
\boxed{v_{\rm state,local}/c_0=\frac{\eta K}{K+C-R}.}
\]

This is the speed of a small coupled disturbance, distinct from energy flux divided by density. In the forward branch d>0 and eta=1, requiring this speed to be at most c0 gives C>=R. At C=R>0 the characteristic speeds coincide and the matrix is defective: it lacks two independent eigenvectors. At R=K+C the time matrix is singular. Neither issue is removed merely by the energy-sum identity.

For K=0.5, R=0.22, C=0.01, eta=1, the local characteristic speed is 1.724137931 c0. Eta=0.5 reduces that particular value to 0.862068966 c0, but does not remove the singular surface or capture-feedback reversal. A negative characteristic speed by itself is not proof of acausality; the relevant failures are the stated propagation requirements and singular/defective behavior.

At the coincident state K=0.5, R=C=0.22, the coordinate matrix has equal diagonal entries 0.694444444 and off-diagonal entry -0.424382716. The standard Fourier solution gives a radiation perturbation proportional to wavenumber times time. At t=1, wavenumbers 1, 10, 100 and 1000 yield normalized radiation gains 0.4244, 4.2438, 42.4383 and 424.3827. These are linear transfer gains, not claims that finite perturbations remain linear at arbitrarily large amplitude. They exhibit the absence of a frequency-independent linear bound in the same norm at this background.

## Can we remain in the good region?

**Conditional homogeneous capture derivation:** with Gamma>0 and no spatial gradients,

\[
\dot C=-\Gamma C/d,\quad R(K+C)=A=R_0(K+C_0).
\]

Starting C0>R0>0, capture lowers C while the radiation reference energy increases. The C=R boundary is reached at

\[
C_*=(\sqrt{K^2+4A}-K)/2.
\]

The time to that boundary follows by integrating

\[
t_* =\int_{C_*}^{C_0}\frac{1-A/(K+C)^2}{\Gamma C}\,dC.
\]

With R0=0.22, C0=1 and Gamma=1:

| K | C at equal-speed boundary | Reference time to boundary |
|---|---:|---:|
| 0.5 | 0.376498204 | 0.709624134 |
| 2 | 0.288409873 | 1.117095682 |

These times are nondimensional, not years. Direct integration stops at the boundary. Thus the strict C>R region is not preserved even by homogeneous capture. For the first case, a subsequent formal singular surface occurs at C=sqrt(A)-K=0.074456265; we did not integrate through it. In the second case there is no positive singular C along the formal homogeneous invariant, but leaving the desired speed domain still occurs.

## Evidence and limits

[run.py](run.py) regenerates [states.csv](states.csv) and [results.json](results.json). The [protocol](protocol.md) separates the original grid from the subsequent domain-preservation test. Run from the repository root:

```text
python research_work/results/advected-temporal-state/run.py
```

All 72 local states were retained. Numerical eigenvalues agree with the analytic values to 2.23e-16. For eta=1, 18/36 grid states have forward speeds at most c0, including nine coincident states; these are not 18 viable models. Eta=0.5 has 28/36 such grid states and no coincidences in this particular grid, not proof that coincidences are impossible. Both grids have eight states with capture-feedback reversal. Matrix-exponential checks of the Fourier solution agree to 1.12e-16. Independent quadrature and direct capture integration agree in crossing time within 4.97e-13; reference-budget error is at most 6.67e-16.

No redshift calibration, survey selection, observation exposure or holdout status changed. No full spatial stability, physical energy/momentum conservation, deposit force or lensing claim is established.

## Next derivation

Replace instantaneous algebraic feedback with explicitly defined finite-response dynamics, or another energy-carrying constitutive rule. Any extra response variable must have a stated energy budget and transport law; it must not silently restore free stationary memory. Derive its characteristic speeds and capture-domain behavior before fitting redshift. Then repeat observable clock, frequency and finite-event timing calculations together. An ultimate explanation of why time exists is not needed for this step.

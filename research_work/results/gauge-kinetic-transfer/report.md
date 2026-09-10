# A photon interaction exists, but it is not the required propagation law

The simplest scalar-photon coupling in our existing candidate action does not make an isolated ideal traveling light pulse steadily lose energy into companions. Its scalar source vanishes for that null electromagnetic wave. A direct backreacting field calculation confirms this cancellation, while colliding pulses do produce scalar energy and conserve total energy and momentum.

This is concrete progress on the mechanism question, but it does **not** derive astronomical redshift, event stretching, lossless companion travel or permanent galactic storage. The collision example also allows energy to return to the electromagnetic field, so its exact dynamics are not our proposed one-way reservoir equations.

## The existing candidate and what is tested

Candidate B in the [action specification](../microphysics/action-specification.md) contains a chi F^2 interaction. That is a **known scalar-photon interaction**, not a novel graviton coupling. Scalars have different spin from ordinary tensor gravitons; no identification between them is made here.

For this dimensionless classical calculation, use signature (-,+,+,+), c=hbar=1 and

\[
\mathcal L=-\frac14 Z(\phi)F_{\mu\nu}F^{\mu\nu}
-\frac12\partial_\mu\phi\partial^\mu\phi-\frac12m^2\phi^2,
\qquad Z(\phi)=e^{g\phi}>0.
\]

The exponential is a specified positive nonlinear completion with the same linear phi F^2 vertex as Z=1+g phi near the vacuum. It changes higher-order interactions and is not silently equated to every completion of the earlier effective action. Such gauge-kinetic scalar couplings are established literature; see, for example, [the scalar-radiation interaction discussed in EPJC](https://doi.org/10.1140/epjc/s10052-025-13814-w). No priority claim is made for this action or the numerical method.

There is no imposed time-dependent refractive index, gravitational well, absorber or external energy source. Gravity is neglected in this local interaction test. Initial scalar field and scalar velocity are zero. The scalar mass m=0.2 and coupling g=0.3 are dimensionless demonstration parameters, **not the fitted lens particle mass or a physical conversion efficiency**.

## Why a traveling wave has no scalar source

Variation of the known action gives

\[
\partial_\mu(ZF^{\mu\nu})=0,
\qquad
\ddot\phi-\nabla^2\phi+m^2\phi=
\frac{Z_{,\phi}}2(\mathbf E^2-\mathbf B^2).
\]

In natural units a null electromagnetic wave has E^2=B^2 and E dot B=0. Thus it can carry nonzero energy while giving no classical source to this scalar. This distinction is known in scalar generation studies, including [low-frequency dilaton generation](https://doi.org/10.1140/epjc/s10052-022-10193-4).

For one transverse potential A(x,t)=F(x-t), we have A_t=-F' and A_x=F'. With phi=0 and V'(0)=0, both Maxwell and scalar equations are satisfied exactly for any smooth waveform F. Its length of travel does not change that cancellation. The scalar source depends on an electromagnetic invariant, not simply on the positive light intensity.

For two counterpropagating waves A=F(x-t)+G(x+t),

\[
A_t^2-A_x^2=-4F'G'.
\]

The cross term need not vanish during overlap. The combined electromagnetic field can therefore generate the scalar. This is an algebraic consequence of the known action, not a new first principle. The scalar must then backreact on the electromagnetic field for conservation to hold.

The exact null-wave statement is restricted to the ideal plane-wave/classical vacuum configuration. It is not a claim that every focused beam, near field, curved background, quantum process or externally driven scalar configuration has zero interaction.

## A closed Hamiltonian calculation

In one spatial dimension, retaining one transverse electromagnetic polarization, define canonical momenta Pi=Z A_t and pi=phi_t. Per transverse area, the Hamiltonian is

\[
H=\int dx\,\frac12\left[
\frac{\Pi^2}{Z}+Z A_x^2+\pi^2+\phi_x^2+m^2\phi^2\right].
\]

It includes the interaction through Z in the electromagnetic energy. The evolution equations used in the code are

\[
A_t=\Pi/Z,\quad \Pi_t=\partial_x(ZA_x),\quad \phi_t=\pi,
\]

\[
\pi_t=\phi_{xx}-m^2\phi+
\frac g2\left(\frac{\Pi^2}{Z}-ZA_x^2\right).
\]

The known total momentum for this reduction is

\[
P=-\int dx\,(\Pi A_x+\pi\phi_x).
\]

No separately prescribed loss term is added. Energy leaves one sector only through the same interaction that changes the other. These are conditional Hamilton equations for the stated candidate, not an astrophysical rate derivation.

## Numerical experiment and results

The domain is periodic with length 100. Gaussian potentials have width 2 and initial centers -15 and +15. One run contains only the right-moving pulse; the other includes the left-moving pulse as well. Evolution stops at time 45, before a repeated collision through the periodic boundary. This is finite energy per transverse area, not a three-dimensional isolated galaxy or a quantum one-photon calculation.

Spatial derivatives use Fourier differentiation and time evolution uses fourth-order Runge-Kutta. Both a 512-point, dt=0.02 run and a 1024-point, dt=0.01 run are retained for each setup.

| Setup | Scalar fraction of initial energy at final time, refined run | Interpretation |
|---|---:|---|
| One traveling pulse | 7.4e-29 | Numerical roundoff, consistent with exact zero source |
| Two colliding pulses | 4.04297e-4, about 0.0404% | Nonzero field generation in this chosen collision |

In the collision run, the largest **recorded** scalar fraction is about 1.093% during the interaction, but the final fraction is much smaller. Some scalar energy returns to the electromagnetic sector. The collision fraction changes by less than 7.5e-10 relative under the paired resolution/time-step refinement.

The largest total-energy drift over all runs is below 4.2e-10 of initial energy. Momentum drift normalized to initial energy in natural units is also below 4.2e-10. The refined runs have total-energy drift below 1.3e-11. These checks support the Hamiltonian evolution; they do not establish a real astronomical conversion rate or long-time irreversible behavior.

The pulses have broad spectra rather than a specified optical carrier. We have not measured a redshift or event-duration change in this experiment. A reduction in total electromagnetic energy alone cannot establish redshift: energy could leave through conversion of whole photons, altered amplitudes or spectral redistribution. The surviving light's frequencies and event timing must be calculated separately.

## Consequences for the accumulation problem

This action contains a way for suitable electromagnetic configurations to exchange energy with a scalar, but an isolated ideal traveling wave in the vacuum does not continuously feed that scalar. Extending propagation distance alone does not supply the missing source. A model relying on ambient fields, other radiation or a dynamical background must specify those inputs and derive the resulting rate, spectrum and detector/clock response.

The scalar produced in the collision is not captured: there is no gravitational well or protected storage channel in this calculation. We have not shown that it becomes the ultra-light stationary source used in the galaxy fits, nor that its energy remains lossless while traveling. Those are separate missing connections.

The reversible microscopic exchange also does not automatically reduce to the proposed monotonic equations for photon loss and permanent deposits. A justified coarse-graining, transport/capture mechanism and account of reverse processes would be needed. This result does not forbid an effective one-way approximation in a specified environment; it shows why that approximation must be derived rather than imposed.

The full research task therefore remains open: derive the source law and its environmental dependence, check spectra and supernova timing, connect transport to supported storage, and test the resulting gravity across galaxies with fresh data after model development. Neither this local experiment nor the earlier disk mass bound evaluates the deferred total photon-supply budget.

Reproduce with `python research_work/results/gauge-kinetic-transfer/run.py`. Requires NumPy. `energy-histories.json` retains both sectors' sampled energies at both resolutions; `results.json` records conservation, scalar generation and refinement. No astronomical catalog or fitted redshift coefficient is used as an input here.

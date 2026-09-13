# Discrete transfer: a joint spectral-width and survival calculation

13 September 2026. Exact solution of a stipulated Markov energy ladder, not microscopic quantum dynamics.

## Result

The earlier smooth energy-loss calculation remains accurate for its small-gap optical examples. Solving the discrete stochastic version additionally couples photon survival to spectral broadening. For an illustrative relative energy width 1e-5 and 90% survival at observed mean redshift z=1, a 2 eV initial photon needs a gap no larger than 2e-10 eV and a removal/useful-transfer rate ratio below 2.11e-11.

That is about one removal per 47 billion useful transfers. This does not exclude the hypothetical theory: the width and transmission are chosen design allowances, and the interaction has not predicted a removal ratio. It shows that the brightness and narrow-spectrum requirements must be met together.

## Model and provenance

**Proposed event process:** allowed photon energies are k*Delta. Initially k=n, E_0=n*Delta. In the dimensionless path coordinate tau=integral alpha ds, a surviving state k makes a useful transition k->k-1 at rate k and is removed at rate r*k. Each useful event sends Delta into the companion sector; removal sends the full remaining photon energy into a separately counted sector.

The energy units k are bookkeeping levels, not a claim that a photon contains smaller photons. The process is a known pure-death Markov chain with a killing channel. Applying it to companion transfer is a hypothesis. It assumes energy-proportional event rates, constant r, no reverse transfer, no frequency-dependent gap, and no correlation between events. Those rates have not been derived from the resonance population.

The k=0 state is a mathematical exhaustion boundary. Formal survival below includes that state; it is not a detectable zero-energy photon. Its conditional probability is negligible in the reported narrow-line design cases. Small-n verification fixtures retain it explicitly to check the mathematics.

## Exact solution

Let P_k(tau) be the probability of an unremoved state with energy k*Delta. **Known master-equation accounting:**

\[
\frac{dP_k}{d\tau}=(k+1)P_{k+1}-(1+r)kP_k.
\]

**Closed-form consequence derived for this application; no claim of mathematical novelty:** define

\[
q=e^{-(1+r)\tau},\quad t=\frac{1+rq}{1+r},\quad
T=t^n,\quad p=\frac qt.
\]

Then P_k=T*binomial(n,k)*p^k*(1-p)^(n-k). Conditional on not being removed, the energy distribution is binomial:

\[
\langle E\rangle_{\rm survive}=E_0p,\qquad
w^2\equiv\frac{\operatorname{Var}(E\mid\mathrm{survive})}{\langle E\rangle_{\rm survive}^2}
=\frac{1-p}{np}.
\]

If z is defined by E_0/<E>_survive-1, then p=1/(1+z), so w^2=z/n=Delta*z/E_0 exactly within this model. This is an energy-centroid redshift and energy RMS width; nonlinear averages in wavelength have their own transformations. In the narrow-line limit the distinction is small, but it must not be ignored for broad distributions.

## Solve both requirements together

At fixed target z and transmission T_min, choose the smallest integer n>=z/w_max^2. The maximum constant removal ratio is

\[
r_{\max}=\frac{T_{\min}^{-1/n}-1}{1-p},\qquad
\tau=\frac{\ln[(1+r(1-p))/p]}{1+r}.
\]

For large n this becomes

\[
r_{\max}\simeq[-\ln T_{\min}]\,w_{\max}^2\frac{1+z}{z^2}.
\]

This approximation is a derived joint design requirement, not a new interaction law. It accounts for the small survivor-selection change in the mean redshift, rather than assuming all survivors follow the deterministic drift exactly.

| Target z | Largest optical gap, approximately (eV) | Maximum removal/useful ratio |
|---|---:|---:|
| 0.00766048 | 2.61e-8 | 1.81e-7 |
| 0.1 | 2.00e-9 | 1.16e-9 |
| 1 | 2.00e-10 | 2.11e-11 |

All rows use initial E_0=2 eV, width allowance 1e-5 and T_min=0.9. These rows separately choose the largest allowed gap; they are not three independently tunable gaps in one universal model. A shared constant gap satisfying all three would need the smallest listed value, which would tighten the other rows' removal allowances. The numerical output also records width allowance 1e-3.

The earlier fixed gap 1e-8 eV at z=1 has width about 7.07e-5 in this model, so satisfying its earlier survival constraint did not satisfy the narrower design width. This is a refinement of the interpretation, not a reversal of its energy accounting.

## Energy ledger and verification

Let C be accumulated companion energy and H accumulated removed-sector energy per initial photon. Exact ensemble rates give dC/dtau=Delta*sum kP_k and dH/dtau=r*Delta*sum k^2 P_k. Thus C=(Delta/r)(1-T) for r>0, with the continuous limit at r=0, and

\[
E_0=T E_0p+C+H.
\]

At the z=1, narrow-width, 90%-survival design point, the outgoing radiation carries 0.900000 eV per initial photon, companions receive 0.949122 eV, and the removed sector receives 0.150878 eV. Its physical destination is still unspecified; the ledger does not assume it is heat or gravitating retained mass.

Nine independent finite-state master-equation integrations, using initial n=5,20,50 and three removal ratios, agree with the exact probabilities and energy sectors to better than 1e-9 normalized error. Probability including the removed state is conserved. Six joint-design cases meet their stipulated mean redshift, width and transmission, and three fixed-gap cases compare against the previous drift approximation. Tests validate this specified process, not its physical realization.

## Research consequence

Keep the near-achromatic resonance rate as an incomplete candidate, now with coupled spectral and survival requirements. A successful microscopic version must derive the required small gap and selective branching from shared couplings, rather than independently choosing them to pass each observable. Its damping and elastic scattering still need their own energy and image calculations.

The model is stationary and its events do not supply a systematic stretch of successive emission intervals. Matching the mean redshift and line width therefore would still not complete the supernova timing requirement. A changing coherent field or a derived time-dependent propagation law remains a separate branch that must eventually be connected consistently. Photon supply, stable deposits, exact-third retention and joint lensing/motions remain open.

Reproduce with `python research_work/results/companion-extensions/discrete-transfer.py`. [Source](discrete-transfer.py), [numerical record](discrete-transfer-results.json).

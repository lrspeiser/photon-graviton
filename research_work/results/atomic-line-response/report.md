# Atomic transition response: one common clock correction is insufficient

## Established inputs and model assumptions

The Dirac-Coulomb spectrum is established physics, not a formula unique to this project. See [Littlejohn, Solutions of the Dirac Equation](https://bohr.physics.berkeley.edu/classes/221/1617/notes/spdirac.pdf) and [NIST's hydrogen spectroscopy compilation](https://www.nist.gov/publications/critical-compilation-experimental-data-spectral-lines-and-energy-levels-hydrogen). The proposed application retains fixed electron mass, charge and matter limiting speed while the optical medium has epsilon=mu=n. Electrostatic Coulomb strength is then divided by n. This explicit preferred-frame matter completion is optional; it is not a covariant completed theory or an inferred property of actual atoms.

Established Dirac energy formula, applied with a proposed effective Coulomb coupling a_eff=a0/n:

E(N,j)/(m_e c_m^2) = [1 + (a_eff/[N-j-1/2+sqrt((j+1/2)^2-a_eff^2)])^2]^(-1/2).

Here N is the atomic principal quantum number, distinct from the optical factor n. a0=1/137 is an illustrative numerical input, not a new measurement. Energies include rest energy; transition differences remove it. This approximation omits radiative corrections, nuclear size and recoil. The field is treated as locally quasistatic for the atomic calculation; rapid nonadiabatic atomic response is not tested.

Established leading spectroscopic scalings applied to this candidate give electronic frequencies proportional to n^-2 and fine-structure splittings proportional to n^-4. Internally computed exact Dirac ratios verify those limiting scalings within the predeclared tolerance for the tested range. Their application here has unverified originality; the spectroscopic structure itself is known.

## Result

In an illustrative homogeneous change from n_e=1 to n_o=1.1, reference propagation has a stretch factor of 1.1. Once measured against the corresponding local atomic transition, the calculated factors are:

| Transition or splitting | Detector/source atomic-frequency ratio | Measured one-plus-redshift |
|---|---:|---:|
| 1s to 2p(3/2), Lyman-alpha component | 0.8264439 | 0.9090883 |
| 1s to 3p(3/2), Lyman-beta component | 0.8264442 | 0.9090886 |
| 2p(3/2) minus 2p(1/2) fine splitting | 0.6830095 | 0.7513105 |

All factors below one indicate blueshift relative to that detector transition, despite reduced photon reference frequency. This reproduces the earlier electronic clock-sign failure and adds a differential spectral diagnostic: the fine splitting and gross electronic transitions do not share one response factor. The fine splitting is a level-frequency comparison, not a claim that a common observed astronomical line directly realizes every tabulated transition.

Established normalization, applied separately to each transition j:

1+z_j = S_ref [nu_j(n_o)/nu_j(n_e)].

Taking ratios between two such inferred factors cancels the common propagation stretch. Therefore adjusting a single propagation coefficient cannot remove the differential atomic effect. Arbitrary mass/charge compensations would be new couplings and would require their own energy and spectral checks. No such compensation is adopted here.

## Verification and limits

Fifteen transition/field combinations were evaluated at 50-digit precision, including the unchanged-field control and four finite changes. Exact and leading-order ratios agree within 0.1 percent; all unchanged-field redshifts vanish, and each nonzero field change produces unequal fine/gross factors. This is a theoretical diagnostic, not a measured constraint or a full atomic numerical solution including QED. Reproduce with `python research_work/results/atomic-line-response/run.py` (requires mpmath).

Matched endpoint fields could leave both source and detector spectra unchanged even while intervening propagation evolves. That protection must be derived with physical environmental coupling and finite boundaries; it cannot be obtained by declaring q_atom fixed. This result does not reject all cumulative-time mechanisms. It narrows the present candidate: an unscreened minimal fixed-parameter atomic completion is incompatible with the intended common positive redshift.

Next join screened emitter/detector regions to a dynamically generated, nonuniform propagation field, retaining matter forces and clock response. The measured redshift test must use local transition frequencies and wavelengths at those endpoints, not just coordinate carrier periods. Fresh observational validation and a derived astronomical rate remain open.

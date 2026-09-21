# CM-1: derive memory from the existing finite companion Hamiltonian

Date: 2026-09-21. Scientific input is the five-sector `common_hamiltonian` in `check_companion_bridge.py` as inspected at `abd026f324c6bf90f71a68f0b0d6a4f8b889c6c3`, blob `13df72233738167c8a56e39ea0a835c675ee5143`. This is a derived finite-model verification, not an astrophysical fit or proof of a novel projection formalism. The matrix, analytic kernel and approximate long-time behavior were explored during design; the checks below are not represented as blinded discoveries.

## Fixed construction

Preserve all five sectors, component energies and couplings. In the rotating frame H'=H-5I, retain P={bound chi plus recoil, receiver plus bound chi plus recoil}, indices 3 and 4, and eliminate Q={source, high photon, traveling chi}, indices 0,1,2. The source term from the eliminated sector's INITIAL state must be retained. No damping, exponential retention time, new bath, independent halo or adjusted matrix element is introduced.

Use the exact projected equation:

```
i p_dot = A p + V exp(-i B t) q0 - i integral_0^t K(t-s) p(s) ds
A = P H' P, B = Q H' Q, V = P H' Q
K(t) = V exp(-i B t) V_dagger
```

For the existing couplings a=.19, b=.17, c=.15, d=.18, only K_00 is nonzero:

```
Omega = sqrt(a*a+b*b)
K_00(t) = c*c*(a*a+b*b*cos(Omega*t))/(a*a+b*b)
```

This is a finite oscillatory kernel, not a decaying bath response. Standard projection/Schur-complement mathematics has prior art (Feshbach 1958/1962); deriving generalized memory equations by eliminating degrees of freedom also predates this project (Zwanzig 1973). The project-specific question is the consequence of its unchanged source/conversion/capture/receiver matrix.

## Verification independent of merely printing the formula

1. Extract the original common Hamiltonian's construction without rerunning or altering its historical best-event search. Require H to equal the sum of its component matrices and verify its stated tridiagonal rotating form to 1e-12.
2. Compare the analytic kernel against matrix exponentiation at 41 points on [0,100], absolute error <1e-12.
3. Compare the projected full resolvent to the Schur-complement resolvent at z=.07+.11i, .33+.05i and -.21+.17i, relative error <1e-10.
4. Use independent 128- and 256-node Gauss-Legendre convolution at times .5,2,20,40,100 to reconstruct eliminated amplitudes and verify the reduced equation against full five-state eigen-evolution, errors <1e-10. The omitted-initial-source control must produce a residual >1e-4, otherwise the test cannot certify the initial fuel term.
5. Record unitarity and component-energy closure on the full trajectories, errors <1e-10.
6. Evolve initial states 0,3,4 on [0,2200] at spacing .05. Record extrema with their sampled times, the exact diagonal-ensemble time-average bound probability, and return probabilities. Extrema are only grid-sampled extrema, not continuous-time certified global bounds. Preserve the original t<=220 peak and explicitly compare peak probability with persistence.

The source/recoil sectors are abstract finite states, not physical positions or an escaping 3D bath. Norm or high peak bound occupancy is not a derived astrophysical capture rate, irreversible retention time, mass density profile or galaxy force. This calculation does not identify the SM-1 AQUAL scalar with chi, nor derive its a0 or constitutive law. It instead supplies an explicit required bridge: the response must be derived from the microscopic Hamiltonian, and high conversion probability cannot stand in for sustained storage.

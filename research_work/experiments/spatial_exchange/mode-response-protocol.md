# SE-1M: gravitational response per unit converted-mode energy

Declared 20 September 2026 before this local calculation; the spatial SE-1 and resolution SR-2 campaigns continue unchanged. This tests whether the new interaction actually changes the field's response per unit energy, using its own Hamiltonian rather than an older gravity formula.

Freeze a locally homogeneous U and set the vector shift to zero. Write kappa=k0 tanh(U/U0), nu=mu^2 exp(2 chi U)(1+kappa^2), a=exp(4U), C=exp(2U). The two-channel mixing matrix has normalized eigenvectors (1,kappa)/sqrt(1+kappa^2) with eigenvalue0, and (-kappa,1)/sqrt(1+kappa^2) with eigenvalue nu. For spatial wave number K, the dispersion is Omega^2=a(K^2+nu_j). Its group velocity is C K/sqrt(K^2+nu_j), bounded by the common front speed C.

For a fixed canonical wave state, define the averaged source-to-energy coefficient S=<partial H_wave/partial U>/<H_wave>. It is a derivative of the declared energy, not a fitted mass multiplier. The scalar-field source receives g S E. Averaging a linear mode over phase gives

S_massless=2,

S_massive=2+[chi+kappa kappa_U/(1+kappa^2)] nu/(K^2+nu).

The derivation uses <Pi^2>=Omega^2 A^2/(2 a^2), <F^2>=A^2/2 and the derivative of the eigenvalue at fixed eigenvector (the matrix Hellmann-Feynman identity). Half the derivative of nu contributes divided by K^2+nu. Crucially, do not recompute the mode state at perturbed U when numerically taking the Hamiltonian derivative: hold its field, canonical momentum and gradient fixed.

Evaluate all288 fixtures: U={-.005,-.003,-.001,-.0003,0,.0003}, chi={-50,0,50,200}, k0={0,.5}, K={.25,.75,2}, branch={massless,massive}. Use mu=.75,U0=.001 and phase samples theta=2 pi j/64. Compute the phase-averaged Hamiltonian derivative by centered U differences of1e-7 and compare S to the expression above, scaled error <1e-5. Independently verify matrix eigenvalues/eigenvectors, positive wave energy and group velocity <=C. Report enhanced, suppressed and negative S values; do not reject or accept a gravity law solely from those signs. Preserve every result.

This is a local frozen-background necessary diagnostic, not a self-generated background, full nonlinear stability result or observation fit. A large S cannot be transplanted into a galaxy as a multiplier: U, mixing, emission energy, mode composition and reaction must evolve in the full model. A negative S is a source-sign prediction of this candidate, not negative total mode energy. Bare X/Y amplitudes are not generally these local propagation eigenmodes.

Attribution: coupled-wave eigenmodes, phase averaging and the Hellmann-Feynman matrix identity are established mathematics. The mass/mixing functions and their complete Hamiltonian derivative are the specific candidate assumptions declared in SE-1. No novelty claim is made for the mathematical tools. All twelve original requirements remain in scope.

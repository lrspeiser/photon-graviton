# SE-ST1: frozen shared-drift compatibility passes

20 September 2026. Protocol167c720; executable7b24c79. All648 declared
frozen principal fixtures pass, including their two particle-mass checks.
Source hash, vectors, eigenvalue diagnostics and particle derivatives are
preserved in shifted-transfer-v1/results.json.

The candidate extension adds bounded common advection to the local transfer
Hamiltonian using kinetic K in the internal-field current:

    H = C^2(Pi_A^2+K^2)/2 + (|grad A|^2+c_Q^2|grad Q|^2)/2
        - beta_j(Pi_A dot partial_j A + K dot partial_j Q) + V,
    K=P-lambda A-epsilon(curl A) cross Q,
    |beta|<c_Q C, 0<c_Q<=1.

This is an explicit candidate change, not the already evolved LR3 model.
Using canonical P instead of K in that current would give a different
Hamiltonian and would not be covered by this positivity argument. Completing
squares bounds the cross terms when |beta|<c_Q C. The frozen directional
energy Hessian S is positive and SM symmetric for the derived principal M.

| Diagnostic | Worst tested value |
|---|---:|
|Minimum symmetrizer eigenvalue|0.0118862|
|Scaled symmetry residual|2.09e-16|
|Longitudinal characteristic residual|1.56e-15|
|Extremal-speed coverage residual|8.89e-16|
|Scaled particle derivative error|1.75e-9|

The provisional particle Hamiltonian is
h_m=sqrt(C)*sqrt(m^2+C|p|^2)+beta dot p. For zero mass its velocity is
beta+C p/|p|; massive particles lie inside that velocity ball. Because the
curl operator vanishes on longitudinal perturbations, the local field has
reference modes beta dot n +/- C in every tested direction, with additional
coupled branches. This resolves the frozen directional coverage question for
this candidate. It does not make all field modes share one speed, exclude
radiation into slower modes or prove the full nonlinear initial-value problem.

The numerical tests use the analytically proposed matrices, not the principal
symbol extracted from an implemented full coupled PDE. The next necessary
step is to implement the full Hamiltonian derivatives, including all derivatives
of C and beta if they depend on scalar/vector fields, and verify the actual
joint equations. Matter must receive its force and supply reciprocal field
sources from the same Hamiltonian. Ordinary-matter production and the fuel
budget of Q/P remain additional unsolved requirements. Locality of numerical
source treatment needs its own control.

Hamiltonian squares, drift dispersion, symmetrizers and particle dispersion
ansatzes are established mathematical structures; no historical novelty is
claimed. No dark matter, expanding background or agreement with an older
gravity law enters these checks. This is a candidate route toward integration,
not a completed matter/light theory or an observational fit.

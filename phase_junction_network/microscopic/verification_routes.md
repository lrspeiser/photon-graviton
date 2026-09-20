# Independent Verification Matrix

**Status:** active protocol  
**Date:** 2026-09-20

## Is a finite Hamiltonian the only way to verify the model?

No. A finite Hamiltonian is one way to establish a microscopic quantum realization, but it is not the only verification method and it cannot, by itself, establish the entire theory.

Different questions require different checks:

| Layer | Independent method | What it establishes | What would fail |
|---|---|---|---|
| Algebraic identity | Exact rational coefficient comparison | The first-order connection eliminates to the Fierz–Pauli kernel for every momentum | Any nonzero coefficient mismatch |
| Numerical algebra | Random floating-point Schur and gauge tests | Implementation correctness and conditioning | Residuals that do not converge to roundoff |
| Finite kinematics | Exact Weyl/CSS rank and commutator calculation | Correct finite constraint count and two local logical modes | Noncommuting constraints or wrong rank |
| Physical spectrum | Projected eigenvalue calculation | Two positive tensor modes and no physical scalar | Extra gapless modes, negative residue, or instability |
| Continuum scaling | Multi-L finite-size scaling | Gap closes as 1/L, with cutoff anisotropy vanishing | Nonlinear gap closing or persistent anisotropy |
| Static response | Lattice Green function and conserved sources | Long-range 1/r potential and source consistency | Yukawa gap, wrong sign, or source dependence |
| Canonical realization | Exact diagonalization, tensor networks, or reduced Hamiltonian | A stable finite quantum phase and cutoff control | Boundary leakage, ghosts, or a gapped phase |
| Covariant realization | Euclidean/Lorentzian path integral and propagators | Correct poles, residues, and correlation functions | Complex uncontrolled measure or wrong propagator |
| Universality | Renormalization/coarse-graining | Insensitivity to regulator details and common infrared fixed point | Fine tuning of every irrelevant lattice choice |
| Empirical test | Frozen predictions against data | Whether nature uses the mechanism | Failure of held-out observations |

No single row substitutes for all the others.

## Checks now completed

### 1. Exact rational first-order identity

`check_first_order_exact.py` uses unnormalized rational bases and no floating-point arithmetic. The frame kernel and the Schur-complement kernel are homogeneous quadratic 6-by-6 matrices in three momentum components.

A homogeneous quadratic kernel is determined by evaluations at

\[
e_x,\ e_y,\ e_z,\ e_x+e_y,\ e_x+e_z,\ e_y+e_z.
\]

The script verifies exact equality at all six points over \(\mathbb Q\). It also checks exact gauge-null identities at four independent integer momentum vectors.

This is stronger than a random numerical test: it proves equality of every quadratic momentum coefficient in the implemented linear theory.

### 2. Random numerical cross-check

The existing 500-momentum floating-point test independently verifies:

- the gamma–gamma identity;
- the Schur complement;
- gauge null directions;
- two equal positive TT eigenvalues.

The exact and numerical implementations use different bases and arithmetic, reducing the chance of one shared coding error.

### 3. Finite-size tensor scaling

`check_tensor_scaling.py` uses the lattice momentum

\[
\widehat k_i=2\sin\left(\frac{\pi n_i}{L}\right).
\]

For the lowest axial mode, the measured fit is

\[
\omega\propto L^{-0.995891},
\]

consistent with a massless linear branch. The maximum relative split between the two tensor polarizations is

\[
6.11\times10^{-16}.
\]

For equal continuum-norm modes \((3,0,0)\) and \((2,2,1)\), the cubic directional anisotropy scales as

\[
L^{-2.01111},
\]

consistent with an \(O(\ell^2k^2)\) cutoff effect that vanishes in the continuum limit.

### 4. Static Green function

The existing 96-cubed lattice calculation fits the point-source kernel to

\[
\frac{A}{r}+B
\]

over radii 4–20 with normalized RMS residual

\[
0.0024096.
\]

This tests a different observable from the wave spectrum.

### 5. Auxiliary stability/no-go test

`check_auxiliary_stability.py` proves and samples the fact that a positive real auxiliary quadratic block with no bare frame term induces a nonpositive Schur complement.

This rejected a naïve microscopic Hamiltonian before expensive exact diagonalization. It also separates three legitimate continuation branches:

- constrained first-order reduction;
- complex Euclidean auxiliary representation;
- bounded Hamiltonian with explicit bare frame stiffness.

## Verification methods still to add

### Reduced constrained transfer matrix

Carry out the Dirac reduction of the first-order frame–connection system and build a finite transfer matrix directly on the physical constrained Hilbert space. Check:

- Hermiticity/unitarity;
- positivity of physical energies;
- 1/L tensor gap scaling;
- absence of propagating connection states;
- finite-state boundary leakage.

### Euclidean two-point function

Compute the gauge-fixed frame correlator

\[
\langle h_{ij}(k,\tau)h_{kl}(-k,0)\rangle
\]

after integrating the connection. The physical tensor poles must be at

\[
\omega^2=c_g^2k^2
\]

with positive residues. Scalar and vector poles must be absent or constrained. Reflection positivity should be tested before interpreting the Euclidean model as a Hilbert-space theory.

### Nonperturbative finite methods

Depending on the final sign structure, use one or more of:

- exact diagonalization on small periodic cells;
- sparse Lanczos methods;
- tensor-network ground states and excitations;
- quantum Monte Carlo if the weights are nonnegative;
- complex-Langevin or contour deformation only with independent convergence checks;
- real-time variational simulation;
- finite-size scaling of gauge-invariant correlators.

### Regulator universality

Repeat the linear and nonlinear calculation on more than one spatial regulator:

- cubic links;
- simplicial/Regge cells;
- irregular isotropic graphs.

A genuine continuum phase should agree after finite parameter matching, rather than requiring a new theory for every graph.

### Empirical verification

Only after microscopic parameters and calibration inputs are frozen should observations be loaded. The first tests remain:

- common photon/gravity cone;
- equivalence principle;
- weak-field light propagation;
- gravitational-wave polarizations and dispersion;
- precision bounds on electromagnetic cutoff operators.

## Current acceptance standard

The project should not declare the Phase Junction theory verified merely because one implementation reproduces the target quadratic kernel.

The next phase passes only when at least three independent lines agree:

1. exact constraint and algebraic reduction;
2. stable microscopic or transfer-matrix spectrum with the correct continuum scaling;
3. source and correlation functions with the correct poles, residues, and long-range behavior.

Empirical data then tests whether the resulting fixed theory describes nature.

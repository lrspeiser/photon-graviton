# CM-1E: a conditional collective route, without demanding irreversible capture

**21 September 2026. Post-CM-1 analytic exploration, executed locally, not a preregistered physical result.** The unchanged five-state Hamiltonian and all original couplings are retained. [Executable](check_companion_ensemble.py); scientific input and exact projection are documented in [CM-1](protocol_companion_memory.md).

## Why this additional calculation matters

CM-1 shows that one source-initialized finite system can reach a high bound occupation and later return almost completely to its source. It does not follow that a macroscopic stationary population requires each microscopic event to be irreversible. We should test a many-source limit rather than build irreversibility into the requirements for every possible companion mechanism.

This calculation supplies one conditional limit, using established diagonal-ensemble and independent-copy statistics. It is not claimed as new statistical mechanics, a microscopic derivation of gravity or evidence that real sources have the required independence.

## Exact stationary ensemble

Write the unchanged Hamiltonian eigenvectors as |E_j> and the initial source state as |S>. Since its five energies are nondegenerate, the infinite-time dephased state is

    rho_bar = sum_j |<E_j|S>|^2 |E_j><E_j|.

It obeys

    [H,rho_bar] = 0,
    Tr(rho_bar) = 1,
    Tr(H rho_bar) = <S|H|S> = 5.

Let P_b project onto the two sectors containing bound chi. The calculated stationary mean is

    p = Tr(rho_bar P_b) = 0.4087849273304729.

This is a stationary ensemble or time average, not a claim that one isolated pure state irreversibly becomes rho_bar. Source fuel and interaction energies are included. The mean is not a physical mass fraction or already a derived gravitational source.

An explicit finite preparation is a population whose source-excitation ages are uniformly distributed on [0,T]. Its state follows by integrating the unchanged unitary trajectory. In the energy basis its off-diagonal elements acquire

    exp[-i(E_j-E_k)T/2] sinc[(E_j-E_k)T/2].

Here sinc(z)=sin(z)/z. No damping or retention time is fitted. For maximum ages T=220, 2200 and 22000 in the original dimensionless model time, the predicted mean bound occupations are respectively **0.41781601925, 0.40798055663 and 0.40860755454**. All three ensembles preserve energy 5. The uniform age distribution is an explicit preparation assumption; it is not derived from stellar histories.

## Many-source fluctuation condition

For N identical independent copies in rho_bar, the total bound occupation N_b has

    mean(N_b) = N p,
    Var(N_b) = N p(1-p),
    RMS(N_b)/mean(N_b) = sqrt[(1-p)/(N p)]
                       = 1.2026113871/sqrt(N).

Thus the conditional occupation RMS is **1.2026 percent for 10,000 copies** and **0.12026 percent for one million copies**. A merely illustrative one-percent occupation target requires at least **14,463 independent copies**. This target is not an astronomical observational bound, and these are occupation statistics, not a demonstrated amplitude or spectrum of gravitational fluctuations.

Independence is a substantive assumption. With a common pairwise correlation coefficient r for the occupation variables,

    [RMS/mean]^2 = [(1-p)/p] [1+(N-1)r]/N.

Positive correlations leave a large-N floor. The illustrative one-percent target would require r below about **6.9143e-5** in that large-N, equal-correlation approximation. A shared frame, finite capacity or interacting sources could violate this; none is included in the independent-copy calculation.

## Executed verification

The local run uses the exact original Hamiltonian blob `13df72233738167c8a56e39ea0a835c675ee5143`. Source SHA-256 for the new script is `9d6054501cb6e7df005aae8c220c9eb7e1d67d1d1a72c2d6824b680ed43942a5`.

| Check | Result |
|---|---:|
| Trace error | 0 |
| Hermiticity error | 3.45e-17 |
| Minimum ensemble eigenvalue | 0.12540838425 |
| Commutator norm with H | 9.14e-17 |
| Initial-to-ensemble energy difference | 0 |
| Component-energy ledger discrepancy | 0 |
| Independent direct-matrix-exponential age quadrature versus analytic sinc state | 4.94e-15 |
| Explicit two-copy mean discrepancy | 0 |
| Explicit two-copy variance discrepancy | 0 |

The quadrature independently integrates 256 direct matrix-exponential pure states over ages [0,220]. The two-copy check uses the full 25-dimensional tensor-product density matrix and occupation operator. All numerical checks pass. This addendum was executed in the working container; unlike the original CM-1 and SM-2 runs, it is not described as a completed GitHub Actions run.

Reproduce from the repository root with a fresh output directory:

    python phase_junction_network/microscopic/check_companion_ensemble.py --output-dir cm1e-replay-001

## Consequence for the next physical construction

A stable macroscopic mean may come from either persistent individual capture or a controlled collective population of reversible exchanges. The latter is now an explicit calculable possibility for the finite model, conditional on its preparation and independence. The previous completion plan's capture requirement should be read as allowing this alternative, not insisting on a one-way microscopic process.

The next spatial model must calculate the source-age distribution, interactions, shared-frame correlations, transport, capacity and actual gravitational response. It must distinguish occupation fluctuations from fluctuations of the physical source operator. The dimensionless stationary fraction does not determine a halo profile, conversion budget, acceleration scale, light deflection or SM-1's constitutive law. Those links remain to be derived, not supplied by relabeling this ensemble.

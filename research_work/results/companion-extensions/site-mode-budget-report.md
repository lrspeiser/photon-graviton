# Connecting protected sites to the Milky Way energy inventory

13 September 2026. Conditional stationary-storage bound, not a fitted lifetime.

## Result and scope

The adjustable site/mode ratios in the previous feedback experiment are extremely small compared with the ratios implied if the entire fitted Milky Way companion mass is stored as tiny finite excitations. With the assumptions below, maintaining 90% protected occupation requires a protected-state mean lifetime above roughly 10^32 to 10^41 years for the 120 kpc examples. This is not a universe-age limit and does not prohibit a perfectly stable state. It quantifies how stable this particular mechanism must be.

We assume the deposited gravitational inventory corresponds to stored energy Mc^2. All that energy is in finite sites of energies 0, E_a and (1-epsilon)E_a, with epsilon=0.5. We do not count pre-existing seed rest mass as energy made by the transfer. If seeds supply most of the gravitating mass, the site count below no longer applies, and their origin and mass inventory must be supplied separately.

## Derivation and formula provenance

**Known mode-counting mathematics:** for g massless polarizations with speed c in volume V,

\[
N_{mode}=\frac{gV(E_{hi}^3-E_{lo}^3)}{6\pi^2\hbar^3c^3}.
\]

This follows by integrating the standard density of states; see [Tong, Quantum Gases, sections 3.1–3.2](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S3.html). Our choice g=2 and application to companion relaxation waves are hypotheses, not observations of gravitons.

**Our conditional finite-site bound:** require P_b at least p>1/2. Because P_a<=1-P_b, energy per site cannot exceed E_a(1-epsilon p). Thus

\[
N_{site}\ge\frac{Mc^2}{E_a(1-\epsilon p)},\qquad
\mu\ge\mu_{min}=N_{site,min}/N_{mode}.
\]

This deliberately allows the remaining sites to be fully bright, maximizing energy per site and minimizing the required site count. It is favorable to the model.

In the preceding one-zone feedback hypothesis,

\[
J=\kappa[(n+1)P_a-nP_b],\quad
\dot P_b=J-\gamma_bP_b,\quad
\dot n=\mu J-n/\tau.
\]

A stationary state with positive protected decay needs J=gamma_b P_b>0 and n=mu tau gamma_b P_b. Even giving the pump all otherwise available sites, P_a<=1-P_b implies

\[
J/\kappa\le (1-P_b)+n(1-2P_b).
\]

For this upper bound to be positive,

\[
n<\frac{1-P_b}{2P_b-1}\le\frac{1-p}{2p-1},
\]
\[
\boxed{\gamma_b^{-1}>\mu_{min}\tau\frac{p(2p-1)}{1-p}.}
\]

At p=0.9 the coefficient is 7.2. This is a necessary condition, not sufficient: finite pumping, reverse optical transitions and finite kappa can tighten it. Kappa cancels only because the sign of the largest possible relaxation current is already a necessary test. For gamma_b=0 the lifetime is infinite and this stationary leakage argument does not constrain the filling time.

## Geometry and numerical results

We use a homogeneous sphere, V=4 pi R^3/3, and the earlier mean free-streaming escape time tau=3R/(4c), represented by exponential escape in the population closure. Modes span E_hi/lo=(E_a/2) exp(+/-b/2). Uniform occupation of every counted mode and full coupling across that band are favorable assumptions. A broad band is a sensitivity case requiring a broadened transition or a site ensemble; a truly sharp identical transition cannot use arbitrary off-resonance modes. Nonuniform occupations and spatial transport require a new calculation.

The two inventories are imported from the coupled-torque calculation with receiver region 30–60 kpc; they remain fitted gravity inventories, not independently established companion masses. For baseline I, M=2.09261e11 solar masses, R=120 kpc, p=0.9:

| Bright energy (eV) | Log bandwidth | Minimum sites/mode | Necessary lifetime above (years) |
|---:|---:|---:|---:|
| 1e-08 | 0.01 | 1.2107e+28 | 2.5587e+34 |
| 1e-08 | 1 | 8.529e+25 | 1.8026e+32 |
| 2e-10 | 0.01 | 7.5667e+34 | 1.5992e+41 |
| 2e-10 | 1 | 5.3306e+32 | 1.1266e+39 |

The narrow-band 1e-8 eV case requires at least 4.244e85 sites, while only 3.506e57 radiation modes are counted. The earlier illustrative ratios of 0.001–1000 therefore do not model this full-inventory branch.

At fixed fractional bandwidth, p and epsilon, this lower bound scales as M/(R^2 E_a^4): larger escape volumes and higher transition energies ease the requirement. Increasing R also changes the transport geometry and cannot automatically represent the observed halo. Increasing energy per storage unit may require collective accumulation because the photon line-width track favors small individual energy transfers. Those are different quantities and must be connected by an explicit accumulation mechanism.

## Consequences and alternatives

1. A truly stable protected state remains possible as a postulate; derive the selection rule or conserved charge that protects it, alongside a permitted formation channel. It does not by itself explain filling, capture or supply.
2. Collective reservoirs may hold many small photon contributions per state. Their mode count, saturation and release law must replace this single-excitation calculation consistently.
3. A separate outgoing channel can reduce stimulated return only if the interaction explains why its waves escape without efficiently reversing the transition. Simply renaming them does not change the coupling.
4. Pre-existing massive seeds can reduce the excitation-energy inventory, but their rest mass is then an additional gravitational source requiring its own accounting.

These are conditional alternatives, not successful fits. The exact one-third empirical reference is unchanged. No radius or lifetime here is a claimed observation or imposed cosmic age.

## Verification

The executable evaluates 72 cases: two inventories, two energies, three radii, two bandwidths and three target fractions. Mode counts agree with dimensionless quadrature within 1e-12 relative tolerance. Radius scaling is checked independently across each matching set. In every case, raising the decay rate 1% above the necessary ceiling makes even the maximal allowed relaxation current negative, so it cannot replenish positive leakage. These tests establish the stated algebraic consequence of the closure, not the physical existence of its states.

Files: `site-mode-budget.py`, `site-mode-budget-results.json`. Predecessor: [escape feedback](escape-feedback-report.md).

# SR-1: reciprocal spatial response

Declared 20 September 2026 before calculation. This extends CC-2; its archived sources and failures remain untouched. No dark component or expansion is introduced. No observation fitting occurs.

## Hamiltonian and canonical equations

Let F=(phi,A), U=g phi-lambda A^2/2, u=dU/dF=(g,-lambda A), alpha=exp(U), z=exp(2sU), a=exp((1+3s)U), d=exp((1-s)U), C=exp((1+s)U). Set b=kappa eta A/sqrt(1+eta^2 A^2), beta=C b. Here s is a fixed constitutive parameter, not an object-specific lens multiplier. Use s=0 control and s=1 candidate.

Hf = integral [a Pi^2/2 + d |grad F|^2/2 - beta dot J + omega^2 F^2/2], J_j=sum Pi_a partial_j F_a.

Hp = sum_i sum_cells W_i [alpha E_i + beta dot p_i], E_i=sqrt(m_i^2+z p_i^2).

With B_ja=partial beta_j/partial F_a = (1+s) beta_j u_a + C partial b_j/partial F_a:

Fdot = a Pi - beta dot grad F.

Pidot = div(d grad F-beta Pi) - omega^2 F - u[(1+3s)a Pi^2/2+(1-s)d |grad F|^2/2] + B^T J - particle_source.

Particle_source_a = sum_i W_i/dV [u_a alpha (E_i+s z p_i^2/E_i)+sum_j B_ja p_ij].

qdot_i = sum W_i [alpha z p_i/E_i+beta]; pdot_i=-sum (partial_q W_i)[alpha E_i+beta dot p_i]. Differentiate normalization of the same spherical W used in deposition. Retain the original nearest-neighbor symmetric gradient energy and centered shift coupling. All field coefficient derivatives must be included.

The frozen longitudinal principal block has characteristic speeds beta dot n +/- sqrt(a d)=beta dot n +/- C. Its symmetrizer diag(a,d) is positive and symmetrizes the evolution block. Photon qdot=beta+C p/|p|; the spherical averaged photon cone is centered on average beta with radius average C. This proves frozen principal hyperbolicity, not nonlinear or orbital stability. The field energy block is positive because |beta|/sqrt(a d)=|b|<kappa<1. The potential coefficient is constitutive, not a derived covariant gravity action. This effective Hamiltonian does not establish microscopic gravitons.

## Preregistered checks and runs

1. Seed 20260920. 240 local samples, s in {0,1}, lambda in {0,.5}: random F,Pi,gradient,p and mass 0 or 1. Test all field, momentum and gradient derivatives against central differences (scaled error <2e-6); principal eigenvalues and symmetrizer residual <1e-10; local positive field block; photon cone <1e-12. Check particle momentum Hessian positivity up to photon longitudinal zero mode.
2. Six coupled grid directional-gradient checks per (s,lambda), n=12,L=8,radius=1.2, seeded compact random fields and nonzero momenta. Include source position and momentum variations. Scaled Hamiltonian directional error <2e-6. Compare s=0,lambda=0 energy/RHS against unchanged CC-2 (<1e-11 scaled). Negative control omits field coefficient derivatives and must have a larger derivative error on a selected field-rich fixture.
3. Six empty-field 3D photon evolutions, lambda=0, n=24,L=12,radius=.9,T=2.5,dt=.02: s=0/1 crossed with eta=0 and source p=.2; eta=.08 and source p=.2; eta=.08 and source p corresponding to 200 km/s with c=299792.458 km/s. Repeat s=1 fast eta=.08 with dt=.01, and n=32,dt=.01 (eight trajectories total). No claim of long-time isolated evolution: stop before boundary contact and record edge amplitudes and conservative travel clearance.
4. Gates: scaled H+Q drift <1e-5, averaged cone residual <1e-10, positive geometric clearance using maximum characteristic speed and maximum source extent; edge amplitude <1e-5. Refined probe-position difference <.01 and field-energy relative difference <.1. Record momentum/angular drift, but do not assert exact grid rotation invariance. Record deflection relative to initial heading, matter energy lost, field energy and circulation, including weakening as well as strengthening.
5. Independently reconstruct final energy from raw arrays without calling model energy/ingredients; verify archives, source hashes, gates and refinements. Preserve first run on failures, no threshold adjustment after results. The weak static limit gives delta=(1+s)2C0/b at fixed U=-C0/r; check this analytically, without describing it as cluster-data agreement.

## Attribution and limits

Metric null geodesics, lapse/shift optical geometry and Hamiltonian evolution are established mathematics: Gibbons, Herdeiro, Warnick and Werner, https://arxiv.org/abs/0811.2877 . s=1 equal weak temporal/spatial response is established, not our invention. Our specific proposal is this effective scalar/vector constitutive family and reciprocal numerical test. Optional amplitude coupling resembles spontaneous scalarization/vectorization; see https://arxiv.org/abs/gr-qc/9602056 and https://arxiv.org/abs/1706.01056 . Lambda is exercised only for derivative controls here; no amplitude-formation claim. This campaign cannot resolve outer rotation, threshold tuning, long-time boundary failures, physical emission, or observational validation. Those remain required by the twelve-item goal.

# Round 24: Casimir-like boundary reservoir and a candidate companion-field EFT

**26 September 2026. Status: registered exploration, not adopted.**

This round records the attempt to replace the hot-companion law's "quiet internal store" with an explicit
Casimir-like boundary reservoir and then close the remaining gap with a field-theory architecture. It does
**not** change the adopted round-12 law, constants, regression baseline, locked forecasts, or collision
prescriptions. No astronomical regression suite was run for this candidate.

The useful outcome is that the problem has narrowed. A microscopic boundary model can reproduce the
required source power scale surprisingly closely, and the heat/collision rule has a simple response-model
derivation. But the companion's ordinary stress-energy is far too small to produce the observed extra
gravity through Einstein gravity alone. The missing physics is therefore a nonlinear gravitational response,
not more companion energy.

## 1. Fixed targets from the adopted law

Use the round-12 constants:

- companion transport speed: u = 169.4 km/s;
- galaxy acceleration scale: a_fit = 6.298e-11 m/s^2;
- strong-field hold scale: g_d = 2.027e-10 m/s^2.

The source power per unit mass implied by the law is

    ell_fit = a_fit u / 2 = 5.3344e-6 W/kg.

These are inputs to this round, not refitted here.

## 2. A spherical boundary reservoir instead of two literal plates

The two-plate Casimir analogy is useful conceptually but its effective spacing was arbitrary. Replace it by a
spherical bag-like boundary with a proton-scale radius. As a benchmark, take the lowest massless-quark bag
eigenvalue x_q = 2.04. Minimizing the simplest three-quark bag energy gives

    R_b = 4 x_q hbar / (m_p c) = 1.7161 fm.

For a vector-like companion boundary mode use two imported benchmark numbers that must ultimately be
re-derived for the actual companion field:

- chi = k R = 2.744 for the lowest vector spherical mode;
- C_8 = 0.51 for the finite eight-vector-mode Casimir coefficient.

Let the dimensionless leakage probability per mode time be the proton gravitational fine-structure strength

    alpha_G = G m_p^2 / (hbar c) = 5.90615e-39.

Then

    E_X       = C_8 hbar u / R_b
    Gamma_X   = alpha_G chi u / R_b
    P_p       = Gamma_X E_X
    ell_micro = P_p / m_p
    a_micro   = 2 ell_micro / u.

Numerically:

- R_b = 1.7161 fm;
- E_X = 33.14 keV;
- Gamma_X = 1.600e-18 s^-1;
- ell_micro = 5.078e-6 W/kg;
- a_micro = 5.995e-11 m/s^2;
- a_micro / a_fit = 0.9519.

Thus this particular spherical eight-vector benchmark lands 4.8% below the adopted astronomical acceleration
scale without inserting a_fit into the microscopic calculation.

This is **not a derivation from QCD**. The values x_q, chi and C_8 are imported benchmarks; using eight vector
channels and alpha_G-strength leakage is the candidate mechanism. In particular, the finite Casimir
coefficient is renormalization-dependent in bag calculations and must not be treated as a precision
prediction.

The boundary mode is best interpreted as a **leak valve**, not the ultimate energy reservoir. Continuous
companion emission must be replenished by the matter sector; static Casimir energy cannot supply free
continuous power.

## 3. Heat and collision suppression from a driven boundary response

Uniform translation of an ordinary passive Casimir cavity does not continuously radiate. The mechanism here
therefore requires the already-existing companion stream to modulate a matter-boundary degree of freedom.

Let the two time-reversed internal channels respond to relative velocities u+v and u-v. If weak parametric
emission is quadratic in the boundary response, then

    P(v) / P(0)
      = ( |u+v|^2 + |u-v|^2 ) / (2 u^2)
      = 1 + v^2/u^2.

For an isotropic population, <v^2> = 3 sigma^2, giving

    P / P0 = 1 + 3 sigma^2/u^2,

which is the adopted heat factor.

Give the internal response q a finite relaxation rate gamma,

    dq/dt = -gamma (q - v),

and let collisions randomize v at rate nu. The stationary filtered variance is

    <q_i^2> = sigma^2 gamma/(gamma + nu),

so the candidate microscopic heat rule is

    k_eff = (3 sigma^2/u^2) gamma/(gamma + nu).

Limits:

- free stars / galaxies: nu << gamma, so k_eff -> 3 sigma^2/u^2;
- rapidly colliding gas: nu >> gamma, so k_eff -> 0.

For an unsuppressed k = 12, nu/gamma = 1, 3, 10, 30 gives 6, 3, 1.091, 0.387,
matching the collision-suppression numbers already seen in the earlier toy model.

What is derived here is the response-model algebra. What is **not** yet derived is the microscopic identity of
q or why the companion stream drives it with exactly this quadratic law.

## 4. The ordinary gravitational effect of the companion energy is far too small

This is the key correction.

For a mass M emitting P = ell M into an isotropic companion stream at speed u,

    U_X(r) = ell M / (4 pi u r^2).

If this energy gravitated only through ordinary Einstein gravity, its equivalent enclosed mass would be

    M_X(<r) = ell M r / (u c^2),

and its acceleration would be

    g_X,GR = G ell M / (u c^2 r).

For a reference M = 1e11 solar masses at r = 10 kpc, using ell_fit:

- ordinary-gravity response: g_X,GR = 1.507e-17 m/s^2;
- required companion pull: sqrt(G M a_fit)/r = 9.369e-11 m/s^2.

The required effect is larger by

    6.217e6.

Therefore the companion's ordinary stress-energy **cannot** be the explanation by itself. The earlier relation
between companion intensity and gravitational acceleration was implicitly a modified-gravity constitutive
law. A complete theory needs to make that nonlinear response explicit.

## 5. Candidate slow transport sector

Introduce a scalar companion variable X and a unit timelike field U^mu defining the local companion rest frame,

    U^mu U_mu = -1,
    h^munu = g^munu + U^mu U^nu.

A minimal quadratic transport action is

    L_X =
      (Z_t/2) (U^mu nabla_mu X)^2
      - (Z_s/2) h^munu nabla_mu X nabla_nu X.

Linear X disturbances then have, in the U-frame,

    u^2/c^2 = Z_s/Z_t = 3.193e-7.

At quadratic order, Z_t > 0 and Z_s > 0 avoid the obvious ghost and gradient instabilities. The hierarchy
Z_s/Z_t is not explained here.

The slow speed u belongs to the **companion transport sector**, not to tensor gravitational waves. A
relativistic completion must keep the physical tensor mode luminal.

## 6. Candidate matter coupling: boundary modulation without a static fifth-force charge

A useful interaction is a quadratic disformal-type operator,

    L_int = (1/(2 M_D^4)) T^munu nabla_mu X nabla_nu X.

Because it is even under X -> -X, there is no linear single-X vertex in this minimal form. The intent is that
matter changes X's local boundary/propagation conditions rather than carrying an ordinary scalar charge.

As a normalization estimate, take the QCD bag scale B^(1/4) = 145 MeV and define

    epsilon = B/M_D^4.

If parametric production has

    Gamma_X / omega_X ~ epsilon^2

and the target leakage is alpha_G, then epsilon ~ sqrt(alpha_G), giving

    M_D = B^(1/4) / alpha_G^(1/8) = 8.71 TeV.

This 8.71 TeV value is a **dimensional target**, not a bound or a completed matching calculation. It depends on
canonical normalization of X, the actual microscopic stress tensor, and the parametric-production matrix
element.

The X -> -X symmetry makes the leading high-energy emission channel pair production,

    p -> p + X + X,

rather than single-X Cherenkov emission. Whether the pair-emission rate is safe for ultra-high-energy cosmic
rays is one of the two decisive calculations below.

## 7. Candidate nonlinear gravity sector

The six-million-fold normalization gap can be closed by making the companion state source a nonlinear
gravitational scalar rather than treating companion energy as ordinary gravitating mass.

In the quasistatic limit consider

    L_phi = - |grad phi|^3 / (12 pi G a) - rho_phi phi.

Variation gives

    div( (|grad phi|/a) grad phi ) = 4 pi G rho_phi.

For a spherical source of scalar charge M_phi,

    g_phi = sqrt(G M_phi a) / r.

Thus the required square-root force law and 1/r falloff follow from the field equation rather than being
declared as a force conversion rule.

A local candidate source is

    rho_phi = rho_cold + (1 + k_eff) rho_free.

If this works, the field equation itself supplies direction and no separate hand-built direction rule is
needed. But this source is **only a candidate completion**. It is not algebraically identical to the adopted
nonlocal S and g_hot construction in nonspherical systems.

The required astronomy test is therefore to rerun SPARC, X-COP, KiDS/SLACS, the Milky Way, the Bullet Cluster,
MACS J0025, Abell 520 and El Gordo using this local equation. If the collision geometry or the galaxy/cluster
normalization is lost, this local EFT completion fails even if its spherical limit is attractive.

## 8. One physical metric for matter and light

The receiver side should not give each material species its own X charge. Instead all matter and light should
couple to one physical metric built from the tensor/vector/scalar gravitational sector. Then composition
independence and lensing are properties of the metric theory rather than separate receiver mechanisms.

Requirements on any relativistic completion:

1. tensor gravitational waves propagate at c;
2. the slow u = 169.4 km/s mode remains a hidden transport mode, not the tensor graviton;
3. photons and nonrelativistic matter respond to the same physical metric;
4. the quasistatic limit reduces to the nonlinear equation above;
5. the kinetic matrix has no ghosts or gradient instabilities over the background used by galaxies and the
   Solar System.

This round does not supply that full covariant metric action.

## 9. A possible origin of the strong-field scale

The spherical leakage linewidth gives a suggestive detuning scale. Across a cavity of radius R_b, ordinary
gravity produces a frequency mismatch of order

    delta omega_g ~ m_p g R_b / hbar.

Equating this to Gamma_X gives

    g_* = hbar Gamma_X / (m_p R_b)
        = 5.88e-11 m/s^2
        = 0.290 g_d.

This is only an order-of-magnitude lead. It does not derive the adopted exponential hold exp(-g/g_d), and it
does not derive the separate 0.15 pc release length. It does suggest that an extraordinarily narrow
matter-boundary transition could naturally become gravity-sensitive near the same acceleration regime.

## 10. What is now actually left

This construction reduces the open problem to a few sharp statements rather than more phenomenological
switches.

### Pass/fail test A: astronomy

Replace the adopted S + g_hot magnitude/direction prescription by the local nonlinear equation

    div( (|grad phi|/a) grad phi ) = 4 pi G [rho_cold + (1+k_eff) rho_free]

and rerun the full frozen suite without refitting object-by-object parameters.

**Fail condition:** the local theory cannot simultaneously preserve the galaxy rotation relation, cluster
normalization, early/late lensing difference, and colliding-cluster lensing geometry.

### Pass/fail test B: particle physics

Canonically normalize X and calculate the energy-loss rate for

    p -> p + X + X

with u = 169.4 km/s and the coupling scale implied by the microscopic matching (nominal target M_D ~ 8.7 TeV),
including the EFT cutoff and the preferred-frame dispersion.

**Fail condition:** observed ultra-high-energy cosmic rays would lose order-unity energy over distances shorter
than their plausible propagation distances.

### Foundation test C: microscopic source

Derive, rather than assume, the effective number of vector channels and the leakage coefficient from a
specific confined matter calculation. The current eight-vector / alpha_G choice is the strongest numerical
lead but remains an input.

## 11. Reproducible arithmetic

The following is sufficient to reproduce the numerical claims in this note.

~~~python
import math

G = 6.67430e-11
hbar = 1.054571817e-34
c = 299792458.0
m_p = 1.67262192369e-27
M_sun = 1.98847e30
kpc = 3.085677581491367e19

u = 169.4e3
a_fit = 6.298e-11
g_d = 2.027e-10

x_q = 2.04
chi = 2.744
C8 = 0.51

alpha_G = G*m_p*m_p/(hbar*c)
R_b = 4*x_q*hbar/(m_p*c)
E_X = C8*hbar*u/R_b
Gamma_X = alpha_G*chi*u/R_b
ell_micro = Gamma_X*E_X/m_p
a_micro = 2*ell_micro/u

ell_fit = a_fit*u/2

M = 1e11*M_sun
r = 10*kpc
g_gr = G*ell_fit*M/(u*c*c*r)
g_need = math.sqrt(G*M*a_fit)/r

g_star = hbar*Gamma_X/(m_p*R_b)
M_D_TeV = 145.0/(math.sqrt(alpha_G)**0.25)/1e6

print("alpha_G", alpha_G)
print("R_b_fm", R_b/1e-15)
print("E_X_keV", E_X/1.602176634e-16)
print("Gamma_X_s^-1", Gamma_X)
print("ell_micro_W_per_kg", ell_micro)
print("a_micro", a_micro)
print("a_micro/a_fit", a_micro/a_fit)
print("ordinary_GR_gap", g_need/g_gr)
print("g_star", g_star, "g_star/g_d", g_star/g_d)
print("M_D_TeV", M_D_TeV)
~~~

Expected central outputs:

    R_b_fm             1.7161
    E_X_keV            33.14
    Gamma_X_s^-1       1.600e-18
    ell_micro_W_per_kg 5.078e-6
    a_micro            5.995e-11
    a_micro/a_fit      0.9519
    ordinary_GR_gap    6.217e6
    g_star/g_d         0.290
    M_D_TeV            8.71

## 12. Status

**Registered, not adopted.**

What this round adds:

- a spherical, proton-scale boundary benchmark that reproduces 95.2% of the adopted acceleration normalization;
- a response-model derivation of the heat and collision-suppression forms;
- a proof that ordinary gravity of the companion energy is about 6.2 million times too weak in a reference galaxy;
- a concrete slow-field kinetic sector;
- a candidate quadratic matter coupling with a nominal matching scale near 8.7 TeV;
- a nonlinear gravitational action whose spherical limit gives sqrt(G M a)/r;
- two decisive pass/fail calculations.

What it does **not** establish:

- that QCD actually generates the proposed companion field;
- that eight vector channels leak with alpha_G strength;
- that the 8.7 TeV matching survives a canonical microscopic calculation;
- that the local nonlinear source reproduces the existing nonspherical astronomy;
- that pair Cherenkov emission is safe;
- a complete covariant, stable relativistic metric action;
- the exact exponential strong-field hold or the 0.15 pc release length.

The adopted hot-companion law and every frozen forecast remain unchanged until those tests are passed.

# A production channel must create occupation, balance recoil and predict color

## Why this follows the bound-cloud exchange test

The [fixed-occupation cloud](../bound-cloud-exchange/derivation.md) can receive photon energy as internal motion but cannot create its own population. Here we construct an optional pair-production channel and calculate its limitations. The receiver is a seeded companion store, so recoil remains in the companion sector. Its initial rest energy is explicitly present; the origin of that seed is not solved or assigned to independent dark matter.

This is a scalar effective model, not a derivation for ordinary gravitons. Bound levels, their spatial overlap and the store's electromagnetic response are hypotheses. We do not claim that a favorable finite-mode transition establishes a fundamental interaction, capture rate or halo.

## 1. Charge conservation can allow pair creation

A relativistic complex scalar has both particle and antiparticle modes. In schematic notation,

```
chi = sum_j [u_j b_j + v_j* d_j^dagger].
Q = sum_j (b_j^dagger b_j - d_j^dagger d_j).
```

The neutral bilinear chi^dagger chi contains b^dagger d^dagger and its adjoint, in addition to occupation terms. It can create a particle and antiparticle together while preserving Q. The sum of their occupations can grow. This is the distinction between conserved charge and total particle count, described in [Tong's complex-scalar quantization](https://davidtong.org/pdfs/teaching/quantum-field-theory/qft2.pdf) and [interacting-field example](https://davidtong.org/pdfs/teaching/quantum-field-theory/qft3.pdf). Those sources support the field-theory distinction, not our bound states or proposed photon coupling.

The previous nonrelativistic single-species density model has no antiparticle production sector. Its occupation theorem remains correct. Keeping pair terms requires a relativistic completion and enough supplied energy, rather than changing that theorem by hand. Opposite scalar charge here is a theoretical label; neither constituent is assumed electrically charged.

## 2. The store supplies recoil, not free energy

An isolated standard photon cannot leave a lower-energy photon plus a massive free pair in otherwise empty flat space: the transferred four-momentum has (k-k')^2 <= 0, whereas a massive pair has timelike total four-momentum. This extends the already derived [vacuum kinematic constraint](../interaction-rate/derivation.md); it is conditional on the stated dispersion and environment.

Instead consider photon + store -> photon' + store containing an additional bound pair. Work in c=hbar=1. Let initial store rest energy be M, additional internal rest-frame energy be delta>0, initial photon energy E, final photon energy E', and scattering angle theta. The final store's invariant mass is M+delta. Four-momentum conservation gives

```
E_threshold = delta + delta^2/(2M).
E' = M (E-E_threshold) / [M+E(1-cos(theta))].
P_store = E n - E' n'.
K_recoil = sqrt((M+delta)^2+|P_store|^2) - (M+delta).
E-E' = delta + K_recoil.
```

For exactly forward propagation,

```
Delta = E-E' = delta + delta^2/(2M).
```

There is no compulsory deflection in this kinematic example, but its likelihood requires an amplitude. The loss for this fixed transition is independent of incoming photon color. Adding binding lowers delta relative to the free pair rest energy; it does not provide extra positive mass. For a bound constituent rest energy m and binding magnitude b, delta=2(m-b). Direct creation into a bound mode already accounts for the negative binding energy; do not also add a second binding-release reservoir. Creation of free particles followed by capture would be a different reaction with its own outgoing channels.

The diagnostic tests 180 combinations of receiving mass, internal gap, photon energy and angle. Energies below threshold are forbidden by the expression, not assigned a negative final photon energy.

## 3. An explicit optional interaction

Extend the earlier [worldline polarizability operator](../matter-assisted/derivation.md) to the companion store:

```
S_int = -(beta/2) integral d tau [chi^dagger chi](X) E_mu E^mu,
E_mu = F_mu_nu u^nu.
```

Here X and u describe the moving store. The operator is electromagnetic-gauge invariant and neutral under the scalar's global phase symmetry. In the store rest frame it contains the electric-field squared; a forward two-photon matrix element can be nonzero for parallel transverse polarizations. This does not borrow the vanishing vacuum chi F_mu_nu F^mu_nu amplitude as a nonzero rate.

For a spatially extended store the worldline expression stands for a mode-overlap integral with a spatial response kernel. That kernel, beta and the bound wavefunctions have not been derived for a self-gravitating companion store. They cannot be taken from the polarizability of ordinary atoms. Projection on two photon modes and one bound particle/antiparticle mode yields the schematic resonant term

```
H_int = g a_low^dagger a_high b^dagger d^dagger exp(i q X) + adjoint,
q = k_high-k_low.
```

The recoil operator transfers q to the store. The adjoint is essential: bound pairs can disappear and return energy to light. This projection retains selected near-resonant modes; it is not a controlled error estimate for omitted modes or a derived continuum scattering probability. Pair levels are presumed available, so their production matrix element is an assumption to test further, not a demonstrated capture mechanism.

## 4. A finite closed quantum calculation

Start with N high-energy photons, no low-energy photons and empty pair modes. A seed store already exists in other modes and has mass M. The reachable basis is

```
|k> = |N-k high photons, k low photons, k particles, k antiparticles;
       store momentum k Delta>,  k=0,...,N.
```

There is no numerical occupation cutoff within the retained Hamiltonian: its photon-number conservation closes this specified mode sector exactly. The projection from the parent operator is itself a physical approximation, not an exact restriction of that full operator. Q=0 and total momentum=N E_high in every basis state. Each pair adds the fixed internal gap delta to the store. That rigid-level approximation omits changes in binding and mode overlap as the store grows; it is not the earlier self-consistent ground-state family with binding proportional to N^3.

Tune E_high-E_low=Delta=delta+delta^2/(2M) so that the first forward transition is exactly resonant. Subtract the constant initial energy M+N E_high. Then

```
H_kk = -k Delta + k delta + K_k,
K_k = sqrt((M+k delta)^2+(k Delta)^2) - (M+k delta),
H_(k+1,k) = g (k+1) sqrt((N-k)(k+1)).
```

The two factors for particle and antiparticle creation multiply to k+1; the photon factors are sqrt(N-k) and sqrt(k+1). Successive transitions are generally detuned because recoil grows and the store mass changes. Tuning the first transition is not proof that an entire radiation population sees the same conversion rate.

We evolve the full finite Hermitian matrix by eigendecomposition and compare its states with independent matrix exponentials. The ledger includes

```
photon loss = <k> Delta,
pair constituent rest-energy gain = 2m <k>,
binding-energy change = -2b <k>,
store kinetic-energy gain = <K_k>,
photon loss = rest gain + binding change + kinetic gain + <H_int>.
```

During interaction, shared interaction energy must be included. No energy enters an unrelated receiver or a hidden pump. The seed rest energy M is reported separately and remains part of the full total. The initial bare product state has zero interaction-energy expectation but need not be an eigenstate of the interacting Hamiltonian; no interaction switching protocol is simulated.

For N=1 exact resonance gives pair probability sin^2(g t). A pair can be produced from empty pair modes, reach unit probability at g t=pi/2, and be reabsorbed by g t=pi. This is production, not irreversible deposition. In the discrete closed model the photon modes remain available to interact. Real escaping radiation can change reversibility, but spatial escape and the associated mode dynamics need to be calculated rather than deleting the adjoint.

### Omitted channels matter for storage lifetime

The parent electric-field-squared operator also contains two-photon creation and annihilation terms. Combined with b d, these allow a stored particle/antiparticle pair to annihilate into two photons. This need not wait for an incoming low-frequency photon: two oppositely directed photons of energy delta/2 each can carry away the internal pair gap while leaving a resting seed store with no recoil. The point operator has nonzero polarization choices for this channel. Actual mode overlap, phase space and the resulting lifetime are not calculated by the finite-mode run.

Consequently, photon number is conserved in the retained two-frequency conversion Hamiltonian, not in the full proposed interaction. Letting the original light escape does not by itself prove permanence. The next rate calculation must compare illuminated pair production and radiative pair loss using the same coupling. The existence of a loss channel alone says nothing about whether the lifetime is short or cosmologically long.

## 5. A simple rate prediction has a color problem

The fixed-gap forward loss Delta already gives a larger fractional shift to lower-energy photons. A color-dependent event rate could in principle compensate. The simplest point-electric-dipole limit of the proposed operator instead has the following leading heavy-store prediction, for a fixed bound transition matrix element:

```
two photon electric-field factors: |matrix element|^2 proportional to E E',
final photon density of states: proportional to E'^2,
Gamma(E) proportional to E (E-delta)^3,  E>delta,
alpha(E) = (delta/E) Gamma(E) proportional to delta (E-delta)^3.
```

The polarization angular integral is energy-independent in this point-dipole, negligible-recoil limit. Photon speed is c=1. Thus Gamma can equally be expressed as a rate per unit path up to a constant target density. The coefficient depends on beta, occupations and bound-state overlap, and is not supplied here. The formula concerns the pair-creation channel with initially empty receiving modes and includes the final photon directions; inverse processes and populated modes need the full kinetic balance. It predicts a strong color dependence, not the required common fractional loss.

This rate approximation is separate from the exact finite-recoil kinematic test and the discrete coherent-mode dynamics. A galaxy-sized store is not pointlike to optical light. Its actual spatial form factor and receiving spectrum may change the dependence and must be derived before applying any rate to observations. The point limit is an explicit candidate failure, not a theorem excluding all bound-pair mechanisms.

## 6. Connection to the proposed time effect remains missing

The effective operator makes a charge-neutral production channel possible and supplies explicit recoil. It does not yet derive n(s,t), explain enhanced coupling in voids, or make event-arrival intervals stretch with photon wavelength. A discrete shift between two frequency modes cannot be asserted to satisfy the earlier whole-signal time map. Next calculate a spatially resolved production spectrum and its effect on the propagation field, using the same interaction for mean energy shift, spectral width, angular change and carrier/event timing. A successful continuum limit must also fund the seed and deposits and reproduce motion/lensing without a gravity multiplier. Preserve all 20 tasks and 32 requirements.

The earlier optical n-field Hamiltonian and this pair-conversion model are distinct candidates at this stage. Do not add their losses for the same photon or count one transferred energy unit twice. A joined theory must derive both propagation and production from one common energy-conserving interaction.

The color constraint above applies when the selected inelastic events are used directly as the redshift law. It does not exclude pair production as a source for a separately calculated evolving mean field. The same parent operator also has density/forward-response terms; any coherent propagation effect and its energy exchange must be derived together with pair production, rather than fitted independently or double counted.

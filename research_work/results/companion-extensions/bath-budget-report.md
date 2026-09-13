# Can thermal companion waves provide the required cold heat capacity?

13 September 2026. Conditional physical bath budget, with no fixed universe size or age.

## Result

A local equilibrium bath of two massless bosonic polarizations traveling at c cannot hold the tested settling heat while remaining at the very low temperatures required by the small-gap protected-state examples. The shortfall is enormous even for the prior orbital cooling energy, which is much smaller than the total deposited rest-energy scale.

This rejects the proposed local thermal-wave heat sink under these assumptions. It does not reject nonthermal escaping radiation, slow collective modes, massive baths, or protection mechanisms with different excitation energies. Each alternative changes the physics and requires its own energy and transport account.

## Known bath physics and proposed identification

For g equilibrium bosonic polarizations with zero chemical potential, linear dispersion E=hbar*v*k, and a sufficiently large continuum volume, **known statistical mechanics** gives

\[
u(T)=\frac{g\pi^2}{30}\frac{(k_BT)^4}{\hbar^3v^3},\qquad
C_V=4Vu(T)/T.
\]

This follows from the density of states and Bose occupation integral; see [Tong's quantum-gas notes](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S3.html). Assigning this bath to companions is a hypothesis. The main calculation sets g=2 and v=c; this does not identify the particles as detected gravitons.

The prior thermal-backflow requirement k_up/k_down<=1e-3 gives T_max=epsilon/(k_B ln1000). Write epsilon=f*E_a, where f is the released fraction of the bright excitation. Keep the previous starting temperature T_i=0.01 E_a/k_B. A closed bath able to accept heat Q requires

\[
V_{\min}=\frac{Q}{u(T_{\max})-u(T_i)},\qquad
R_{\min}=(3V_{\min}/4\pi)^{1/3}.
\]

These are derived budget requirements, not new gravity laws or assumptions about the universe's boundary. Raising an arbitrary constant heat capacity is no longer available: here its value follows from the bath spectrum and volume.

## Two energy accounts, kept distinct

1. **Previously calculated orbital cooling:** use the coupled 30–60 kpc receiver cases, Q=4.471e50 J and 4.432e50 J for ordinary-matter baselines I/II. This asks what would happen if that heat were assigned to this bath. The orbital model has not independently identified this as its cooling channel.
2. **All deposit energy supplied by protected excitations:** only if the modeled inventory's energy M_dep*c^2 consists of those excitations, retaining fraction 1-f requires release Q=[f/(1-f)] M_dep*c^2. This is a deliberately explicit additional identification, not silently adopted. It does not apply to transitions that merely perturb a separately existing massive seed. Such a seed's original gravity could not be counted as created by the small transitions.

The reference inventory is approximately 2.09e11/2.06e11 solar masses from the preceding model, not a new mass measurement. Container radii 30,120,1000 kpc are diagnostic local volumes. The total inventory is a normalization benchmark; this calculation does not claim its entire original profile lies inside each container.

## Representative numbers

The following table uses baseline I, f=0.5 and a 120-kpc bath. Both tiny excitation energies come from earlier illustrative transfer calculations, not from a measured companion spectrum.

| E_a (eV) | Maximum temperature (K) | Heat capacity over allowed temperature interval (J) | Orbital cooling / capacity | Required bath radius for orbital cooling (kpc) |
|---|---:|---:|---:|---:|
| 1e-8 | 8.40e-6 | 8.01e29 | 5.58e20 | 9.88e8 |
| 2e-10 | 1.68e-7 | 1.28e23 | 3.49e27 | 1.82e11 |

If all deposit energy instead came from the protected excitations, the f=0.5 formation heat would be 3.74e58 J. The corresponding 120-kpc shortfalls are 4.67e28 and 2.92e35. These are conditional formation costs, not additions to the orbital heat already counted or a demonstrated source history.

An arbitrarily large universe does not increase the capacity of this local bath. Spreading heat into a larger volume requires a spatial transport and causal history. The enormous required radii are not excluded merely by comparison to a conventional cosmological horizon; no such horizon is imposed here.

## Can it radiate the heat instead?

For an isotropic equilibrium bath with a freely emitting surface, the **known radiation flux** is v*u/4. At v=c, a spherical surface constrained to remain no hotter than T_max has the conditional net thermal power

\[
L_{\max}=\pi R^2c[u(T_{\max})-u(T_i)].
\]

For the 120-kpc, E_a=1e-8 eV example, this is 4.86e16 W. Exporting the prior orbital heat at that maximum rate takes at least 2.91e26 years; the smaller-gap example gives 1.82e33 years. These are required durations under this surface-temperature constraint, not an assumed universe age or an observational age exclusion. They show how much history this proposed thermal outlet would need. Nonthermal or hotter emission does not obey this particular cold-thermal-outlet estimate and must be calculated separately.

## Alternatives quantified without adopting them

Because u scales as g/v^3, the script records the degree count or propagation speed required to fit the heat into each benchmark volume. For the 120-kpc, larger-gap orbital case, holding g=2 would require v approximately 36.4 m/s, or holding v=c would require g approximately 1.12e21. These are inverse targets, not discovered properties. Slow material-like modes would need a supporting medium, a valid low-energy dispersion range, actual thermalization and their own gravitational source. Such a bath would not be identical to the freely traveling c-speed companions.

Direct escape of nonthermal relaxation quanta could avoid storing heat locally. That does not automatically guarantee protection: their local occupation, reabsorption, stimulated return and angular escape must be computed. A thermal temperature cannot be assigned to those waves merely to reuse the present formula.

## Verification and decision

The code evaluates 24 energy/gap/release-fraction/baryon cases, each with three volumes. The g=2,v=c coefficient agrees with 4*sigma_SB/c, providing an independent physical-constant check; capacity and required-volume identities are checked. No fit parameters are chosen to match observed speeds or redshifts in this budget.

Do not use a small local equilibrium massless bath as the unexplained large heat capacity in the protected-state model. Retain explicit nonthermal escape and physically specified slow or massive baths as separate options. The next useful calculation is the escape and reabsorption of emitted relaxation waves, with their population determined by the emission rate and volume. Source supply, lifetime, selection rules, gravitational support, timing and joint lensing/motions remain unresolved.

Reproduce with `python research_work/results/companion-extensions/bath-budget.py`. [Source](bath-budget.py), [complete results](bath-budget-results.json).

# A reciprocal interaction that makes assembly bias explicit

13 September 2026. Candidate mechanism and conditional rate derivation.

## Outcome

A proposed two-mode energy-exchange interaction gives a specific condition for the favorable forward bias used in the assembly calculation: the higher-energy input mode must be more occupied than the lower-energy output mode. Simply supplying a smooth thermal spectrum in the same direction does not provide that condition. Distinct source and escape channels may do so, but their angular, spectral and polarization geometry must be specified.

This is a local loading mechanism. It does not yet implement the separate traveling-companion leg of the full theory. The carriers could instead be companions, but their input populations and spectrum would then need calculation. Neither reading establishes cosmological redshift, a capture cross-section, lensing or event-duration stretching.

## Hamiltonian and provenance

**Proposed effective interaction**, using familiar ladder-operator mathematics:

\[
H_{int}=\hbar g S_+ a_L^\dagger a_H+\hbar g^*S_-a_H^\dagger a_L,
\qquad E_H-E_L=\delta.
\]

S_+ raises the store by delta; a_H removes a high-energy quantum and a_L dagger creates a lower-energy quantum. The conjugate term permits reverse transfer. This is an effective resonant interaction of the same general mathematical type used for frequency conversion, not a claim to have invented ladder operators or derived a fundamental graviton coupling.

Each forward event conserves energy: E_H=E_L+delta. A translationally complete interaction must also include store recoil or a spatial field carrying the momentum difference. The present internal-state Hamiltonian omits it; resonance would require the corresponding recoil correction. Calling the operators photons does not resolve that requirement.

**Conditional weak-coupling rate derivation:** assume independent incoherent modes, equal paired phase-space weights and flat storage matrix elements. Bosonic number-state amplitudes give

\[
a=\Gamma n_H(1+n_L),\qquad b=\Gamma n_L(1+n_H),
\]
\[
\boxed{r=b/a=\frac{n_L(1+n_H)}{n_H(1+n_L)}},\qquad
r<1\ \Longleftrightarrow\ n_H>n_L.
\]

The factors n and n+1 are known bosonic mathematics. The choice of interaction, constant Gamma, selected paired modes and storage ladder are our hypotheses. Integrating over unequal mode densities or transition overlaps can change the aggregate rates, so this equation cannot be applied directly to total astronomical fluxes. Correlated or coherent fields also require a different treatment.

For a desired r, the low-mode population must be

\[
n_L=\frac{r n_H}{1+(1-r)n_H}.
\]

In dilute modes r is approximately n_L/n_H. A 10% population contrast can thus produce the earlier r=0.9 example. At high occupation the stimulated factors weaken the bias unless the contrast is larger.

## Smooth thermal spectra do not supply that bias

For an undiluted equilibrium occupation n(E)=1/[exp(E/kT)-1],

\[
r=\exp(\delta/kT)>1.
\]

This follows from the standard Bose occupation, described in [Tong's quantum-gas notes](https://www.damtp.cam.ac.uk/user/tong/statphys/statmechhtml/S3.html). Equal dilution of both selected modes preserves n_H<n_L and hence r>1, though it changes the magnitude. This conclusion concerns occupation per mode, not the slope of plotted energy intensity.

The executable checks illustrative temperatures of 3000, 6000 and 10000 K, two small energy transfers, and dilution factors 1 and 1e-9. At 6000 K and delta=1e-8 eV, r-1 is about 1.93e-8 without dilution and 1.98e-8 with dilution. Thus the spectrum does not generate the desired forward bias, but its small reverse bias should not be confused with the much larger 3% reverse bias in the previous sensitivity example. These are Planck-shaped examples, not fitted stellar spectra.

A cold receiving ladder can still absorb transiently even for r>1; that inequality does not mean every net transfer is forbidden. It means the constant link rates favor downward steps once both are available. Repeated first-passage filling and endpoint protection remain separate questions.

## Finite reservoirs consume their own bias

Let M_H and M_L count modes in closed, internally uniform reservoirs. After x net forward events,

\[
n_H=n_{H0}-x/M_H,\qquad n_L=n_{L0}+x/M_L.
\]

Equality occurs at

\[
x_{equal}=\frac{n_{H0}-n_{L0}}{1/M_H+1/M_L}.
\]

This is the point where paired-rate bias changes sign, not a universal maximum stored energy. Stochastic transfer, a cold initial ladder or an absorbing endpoint can continue changing storage after this point.

For equal mode counts, n_H0=1e-12 and the initially r=0.9 low population, equality occurs after only 5% of the initial high-mode photon number transfers. Those events store only delta per photon, not its full optical energy. A sustained bias therefore requires replenishing high modes or exporting low modes. The source and outlet must remain in the energy ledger. Net x events remove x E_H from the high reservoir, add x E_L to the low reservoir and add x delta to the store.

## Next discriminating calculation

Specify physically distinct source and escape modes and follow their occupation during loading. Test whether the required outlet preserves observed image directions and spectra; an empty output channel cannot be assumed for every ray without checking what channel receives that ray. Derive Gamma and its frequency dependence, include momentum, and connect the field populations to the available radiation supply. Only then can the favorable assembly times be treated as a prediction rather than a freely chosen rate ratio.

This is a possible non-equilibrium loading route, not a completed theory. It leaves the exact one-third empirical retention reference unchanged.

## Verification

Twelve occupation designs verify inversion for target rate ratios and finite-reservoir equalization. Twelve thermal cases check the sign and the undiluted detailed-balance identity. Nine finite number-state matrix elements independently verify the forward stimulated factors and the Hermiticity of the interaction. These are checks of the proposed model's algebra, not measured evidence for its existence.

Executable: `mode-driven-bias.py`; output: `mode-driven-bias-results.json`. Predecessor: [assembly and recycling](assembly-recycling-report.md).

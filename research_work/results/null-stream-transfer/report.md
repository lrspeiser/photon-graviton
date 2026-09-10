# Massless companion transport repairs a stress identity, but not event timing

An alternative kinetic model can send lost photon energy into a co-propagating massless stream with consistent energy and momentum transport. It passes the specific stress/center-of-energy test that the earlier canonical reservoir failed. It does not stretch the time between emitted events: all signals still traverse the same fixed path at the same speed. The transfer coefficient cannot fix that timing identity, even if it changes with time.

This is a joined transport-and-timing test of a different mathematical model. It is not a microscopic photon interaction, a graviton discovery, a complete relativistic field theory, or an adopted explanation of gravity. The old canonical model and its failure remain unchanged.

## Formula provenance

**Known relativistic radiation accounting.** In units c=1, take a forward null direction k=(1,1,0,0). A unidirectional massless stream has stress tensor T^{mu nu}=u k^mu k^nu. Its energy flux and momentum density agree; in ordinary units J_energy=c^2 p_density. Local energy and momentum accounting are discussed in [Feynman's field-energy lecture](https://www.feynmanlectures.caltech.edu/II_27.html). This standard stress structure is not unique to our project and does not determine particle spin.

**Hypothetical kinetic postulate, using standard transport mathematics.** Let eta(t,x) >= 0 transfer photon energy to companions while both streams propagate at c=1:

    (partial_t + partial_x) u_gamma = -eta u_gamma
    (partial_t + partial_x) u_c     = +eta u_gamma
    T_gamma = u_gamma k k
    T_c     = u_c k k

The opposite source four-vectors are Q_gamma=-eta u_gamma k and Q_c=+eta u_gamma k. Their sum vanishes, so the total displayed radiation stress is symmetric and locally conserved. There is no deposited reservoir in this transport test. Already-created companion energy has no loss term and continues onward.

The coefficient is prescribed in one reference frame. Under a boost along the beam, with d=gamma(1-beta), u'=d^2 u and eta'=eta/d give the correctly transformed source four-vector. This tensor consistency does not derive the physical origin or dynamics of the rate-setting environment. Nor does it include gravitational backreaction or an action for the interacting fields.

**Additional postulate required to call the energy loss redshift.** Photon number density obeys (partial_t+partial_x)N_gamma=0. Then per-photon energy E_gamma=u_gamma/N_gamma falls, rather than merely removing photons. Applying the known relation E_gamma=h_P nu assigns a lower frequency. Photon-number retention and this kinetic frequency interpretation are assumptions here, not a derived coherent electromagnetic wave solution.

**Known fractional-loss solution, not a new formula.** Along a ray,

    E_out/E_in = exp(-integral eta dt)
    1+z_assigned = exp(integral eta dt)

In ordinary units a distance coefficient would be alpha=eta/c along that ray. The exponential is standard integration. No new first-principles alpha or astronomical normalization has been obtained.

## Executed transport and moving-observer checks

A finite pulse initially occupies -2 <= x <= -1. Conversion occurs only within 1 < x < 3, with W(x)=sin^4[pi(x-1)/2] there and zero elsewhere. The tested rate is eta=a W(x)(1+b t). These arbitrary bounded profiles let the pulse leave the interaction region; they are not a derived gravitational screen or a fit to void observations. Units are normalized, not megaparsecs or cosmic years.

| Case | Final photon energy | Final companion energy |
|---|---:|---:|
| No conversion | 1.000000 | 0.000000 |
| a=0.2, b=0 | 0.860708 | 0.139292 |
| a=0.2, b=0.04 | 0.842822 | 0.157178 |

Both 129- and 257-label integrations reproduce the analytic exit fraction. Maximum local energy-transfer error is 8.9e-16; maximum exit-fraction error is about 2.1e-12. The total pulse center follows x=-1.5+t with maximum error 5.4e-15.

Five boosts, beta=-0.8,-0.3,0,0.3,0.8, are checked for each of the six pulse runs. Integrals use surfaces of equal boosted time, respecting relativity of simultaneity. Total boosted energy and momentum are both d times the original energy; the center-of-energy speed remains one, with maximum error 1.5e-14. These 30 checks verify the displayed null-stream accounting.

This differs materially from the old canonical reservoir. That model already transferred missing photon momentum into its profile support, while its outgoing reservoir energy lacked the corresponding canonical momentum. Adding massless momentum to that reservoir without changing the support dynamics would double-count momentum. Here the transfer equations are replaced: both energy and momentum move between radiation species, and the old support impulse is not retained. A microscopic implementation might introduce additional medium forces; those would need explicit accounting.

## Why timing still fails

**Conditional kinematic identity, using known constant-speed propagation.** For fixed source/detector separation D and ordinary equal endpoint clock standards,

    t_arrival = t_emission + D/c
    d t_arrival / d t_emission = 1

The rate eta does not appear. Twenty-four independently integrated rays and eighteen event-interval comparisons confirm this for zero, stationary and evolving conversion rates; maximum event-stretch error relative to one is 5.4e-15. A changing eta can change the received energy and distort a brightness profile, but it does not change the arrival mapping of identifiable source stages in these equations. Such distortion is not automatically the measured stretching of spectral evolution.

There is also an unresolved coherent-wave requirement. **Known free-wave identity:** a smooth forward wave in fixed-speed one-dimensional vacuum has phase phi=F(t-x/c), so its frequency is constant along a forward characteristic. The imposed per-photon law d nu/dt=-eta nu therefore cannot simultaneously be the unchanged free Maxwell phase law in the conversion region. An interaction can modify electromagnetic propagation, but its phase/coherence dynamics must be supplied; the kinetic energy ledger does not derive them. This observation does not rule out interacting waves or more general nonexpanding models.

## Connection to existing observations

The 35 previously exposed spectral-aging rows are reused without refitting or reopening holdouts. Their provenance is the existing Table 3 transcription from [Blondin et al. (2008)](https://arxiv.org/abs/0804.3595), checked in the earlier electromagnetic audit.

Using the recorded individual uncertainties, the unchanged diagonal chi-square is 150.569 for no event stretching and 26.949 for aging rate 1/(1+z). These exactly reproduce earlier scores; they are not new independent evidence. The comparison is conditional on the spectral-age templates, does not include full shared covariance, and uses measured z rather than predicting distance/redshift from eta. It shows why the repaired transport must still confront the existing timing requirement.

## Consequence for the full objective

The specific missing-momentum problem is repairable at the level of this alternative radiation stress model. That pass is necessary but insufficient. The next physical candidate must derive an interaction that changes the arrival map as well as the photon spectrum while accounting for physical energy, momentum and actual clock standards. Assigning eta a more flexible fitted function cannot change t_arrival=t_emission+D/c while these trajectories and endpoint clocks stay fixed.

A particle-transition rate, coherent electromagnetic response, companion identity, capture, retention in wells, spatial support, and the joint motion/lensing response remain underived. The galaxy field is still an empirical response. The total photon-supply budget remains deferred, not passed. No stellar weight, gravitational coefficient, expansion assumption, or reserved observational outcome was changed.

Run `python research_work/results/null-stream-transfer/run.py` from the repository to reproduce the pulse, boost, ray and exposed-data checks. Full numerical outputs and source hashes are in results.json.

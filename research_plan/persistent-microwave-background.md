# Persistent microwave background and collective spatial modes

Status: user-requested extension of the hypothetical nonexpanding photon–companion theory. Candidate mechanisms are adopted for investigation, not asserted as measured causes. No CMB fit has been performed for this extension. No fixed cosmic age, curvature radius, or boundary is adopted. The six existing research objectives remain intact; this supplies an additional joint constraint on propagation and the mechanism-specific prediction.

## Physical idea

The universe may contain a persistent thermal radiation bath whose angular structure comes from collective modes or spatial organization. The radiation, the disturbances that affect it, and the statistical pattern on the sky are distinct. A persistent glow need not imply that the original disturbances remain active, and persistent statistics need not imply fixed individual hot and cold patches.

The observed acoustic peaks concern angular patch sizes, not a succession of eight arriving waves or eight microwave frequencies. A 2016 analysis reported eight TT peaks, five EE peaks and twelve TE extrema; these counts are measured features of that analysis, not a fundamental limit on modes ([Pan et al.](https://arxiv.org/abs/1603.03091)). The observations are targets; the paper's expansion-based interpretation is not a premise of this branch.

An arbitrarily old or eternal universe is permitted, not required. A finite-age transition in an eternal universe is also permitted. Infinite available time alone neither synchronizes modes nor establishes a thermal spectrum.

## Three alternatives retained separately

| Branch | Proposed origin of the pattern | Physics to specify |
|---|---|---|
| A: recording transition | A medium supported correlated disturbances; a finite-age transparency transition left radiation carrying their pattern. | Correlation/phase mechanism, transition, thermalization, and preservation during later conversion. The transition is not identified with creation of the universe. |
| B: driven resonances | Collective modes remain statistically sustained and influence microwave emission, scattering or frequency shifts. | A physical medium/field, driving spectrum, damping, energy source, and temperature/polarization coupling. |
| C: spatial organization | A persistent thermal bath acquires angular structure from a correlated medium. | A dynamical origin for spatial correlations and a transfer law, rather than freely inserting each observed peak. |

These alternatives are not silently combined into a single established mechanism. Each needs a forward calculation and frozen predictions.

## Optional spherical geometry

A physical ball with a boundary requires explicit absorbing, emitting, reflecting or transmitting boundary conditions and observer-location predictions. No such wall is adopted here.

The boundary-free alternative is a round three-sphere S^3 of curvature radius R. Its geometry has no preferred center. It is an optional toy background, not a derived static gravitational solution. Static radius and stability still need a supporting stress-energy/dynamical law. Existing static-universe stability results depend on perturbation type and matter assumptions; they do not license a blanket stability claim ([Barrow et al.](https://arxiv.org/abs/gr-qc/0302094)).

**Known mathematical identity:** for radius R,

    Laplacian Q_nA = -n(n+2)/R^2 * Q_nA,
    n=0,1,2,...; degeneracy (n+1)^2.

This is established S^3 spectral geometry ([Lachieze-Rey](https://arxiv.org/abs/math/0401153)). For an ideal sound-like scalar with speed c_s, the derived toy frequency is

    omega_n^2 = c_s^2*n(n+2)/R^2.

The sequence is not truncated at eight. n is a three-dimensional spatial label, not the observed angular multipole ell. A universe-scale collective mode is not the microwave photon's individual wavelength. Circulation around compact space does not itself cool photons. Previous flat/specified-distance calculations cannot be relabeled S^3 results without recomputing ray paths, angular distances and capture geometry.

## Equations for a calculable branch

**Known driven-oscillator framework, with hypothetical coefficients/couplings:** for a real, unit-inertia mode amplitude a_nA,

    d²a_nA/dt² + 2 gamma_n da_nA/dt + omega_n² a_nA = F_nA(t),
    E_nA = [(da_nA/dt)² + omega_n² a_nA²]/2,
    dE_nA/dt = F_nA da_nA/dt - 2 gamma_n (da_nA/dt)².

Physical energy normalization must be supplied by the field model. The equation allows ideal lossless persistence, decaying transients or balanced driving. Random driving can excite a resonator; it does not automatically generate the observed TT/TE/EE angular structure.

**Known angular-statistics definition:** expand the observed temperature or E polarization as a_lm^X, and form

    C_l^{XY} = (1/(2l+1)) sum_m <a_lm^X (a_lm^Y)*>, X,Y in {T,E}.

For independent spatial modes with covariance P_n and transfer coefficients T_lm,nA^X, a schematic discrete forward model is

    a_lm^X = sum_nA T_lm,nA^X a_nA,
    C_l^{XY} = (1/(2l+1)) sum_m,n,A P_n T_lm,nA^X (T_lm,nA^Y)*.

The transfer coefficients must include geometry, visibility, propagation and polarization generation; normalization follows the chosen mode basis. Continuously driven branch B generally requires a double integral over source times and the unequal-time covariance <a_nA(t) a_nA*(t')>, not merely a snapshot P_n. These are general linear-response identities, not a computed prediction. A compact shared law must determine P_n/driving; assigning a separate strength to every peak is not a mechanism test.

**Known thermal target, not a derivation of its origin:**

    u_nu(T) = (8*pi*h*nu^3/c^3) / [exp(h*nu/(k_B*T))-1].

The model must produce the observed near-thermal microwave spectrum around 2.725 K and its distortion limits ([COBE/FIRAS archive](https://lambda.gsfc.nasa.gov/product/cobe/firas_overview.html)). Redshifting arbitrary starlight alone does not generally turn it into this distribution. Thermalization and pattern generation must be compatible: interactions that equilibrate the spectrum may erase directional structure. A lossless equilibrium bath can persist without continuous driving, but that statement does not apply unchanged if our conversion/capture continually removes energy from it.

## Connection to photons, companions and deposits

The established project postulate remains photon energy transfer to traveling companions, followed by capture in wells. Collective pattern energy must be accounted for without counting it again as independent companion energy if it is an excitation of that same field.

For illustration, define disjoint photon, freely traveling companion, deposited and material-reservoir energies. In a closed-volume bookkeeping model,

    dot E_gamma = P_star - P_conversion + P_return,
    dot E_companion = P_conversion - P_capture,
    dot E_deposit = P_capture - P_release,
    dot E_material = -P_star + P_release - P_return.

Their sum is constant. These equations are **conservation bookkeeping**, not supplied rates or a new mechanism. P_release/P_return are optional new channels; set them to zero to retain permanent deposits and no recycling. They are not adopted merely by adding this extension. If organized modes exchange energy with any sector, add equal-and-opposite transfer terms or include the modes inside that sector's energy.

A closed finite-volume system with permanently positive capture and no release cannot have a stationary finite deposited reservoir forever: dot E_deposit=P_capture>0. Eternal steady recycling therefore requires an explicit departure from permanent capture, or a different nonstationary/open description. Infinite age and repeated circulation do not create usable energy or resolve entropy/free-energy requirements.

The existing [homogeneous optical-growth supply mismatch](../research_work/results/brightness-distance-consistency/optical-growth-report.md) remains on record. A persistent background does not remove it by assertion. This extension may replace that restricted source interpretation, but must calculate the replacement's energy flow and its effect on the existing redshift, timing, brightness, rotation and lensing formulas.

## Formula-first work sequence

1. Specify one of A/B/C, geometry, field, a small shared parameter set, source/visibility law and energy transfers. Preserve alternatives as separate branches.
2. Calculate the frequency spectrum and joint TT, TE and EE spectra. Include peak positions/heights, damping and signed cross-correlation; do not use peak count alone as success.
3. Declare calibration and prediction partitions before fitting. Fit an explicitly selected TT range, then freeze shared parameters and predict the remaining TT range and polarization. Parameters invisible to TT need independent constraints rather than a claim of a unique TT-derived polarization prediction.
4. Compare with published measured spectra/covariance while recording foreground, mask, calibration and inference assumptions. Do not use cosmological parameter chains as independent observations or insert expansion-derived distances as facts.
5. Carry the same propagation and energy law back into the six existing goals. Check distortion, angular smearing, time variability, topology and observer-location consequences when predicted. Never label an already viewed spectrum untouched.

No new data-resolution campaign or claimed CMB success is created by this addition. Its deliverable is a common formula that fits selected observations and predicts others with parameters frozen.

## Provenance

Integrated from the user's supplied discussion of spherical/eternal backgrounds, resonances and the three candidate constructions. The supplied text is a proposal and explanatory source, not an empirical result. The primary sources linked above were checked for the mathematical and observational claims used here. The original explanatory discussion also contrasted standard expanding-universe histories; those histories are not adopted as physical premises of this extension.

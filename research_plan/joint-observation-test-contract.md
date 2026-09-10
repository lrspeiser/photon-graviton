# Forward contract for distinguishing nonexpanding mechanisms

Status: 2026-09-10 working protocol. This is not a claim of completed observational validation or a finalized preregistration. The goal is one predictive mechanism, not separately adjustable equations for each phenomenon.

## 1. Require predictions before ranking

Every candidate must specify a local propagation or interaction law, its clock/rod convention, companion energy/momentum ledger, source conditions, and parameter vector. Derive wavelength shift, event-duration ratio and relative messenger delays from that law. Record "unspecified" where derivation is missing. A Fourier receiver map does not by itself specify a causal flight-time map.

Initial findings: [joint audit](../research_work/results/joint-propagation-audit/report.md). The uniform positive photon-only inverse-affine rule is incompatible with the associated GW170817/Gamma event under its stated assumptions. A shared messenger rule remains open, not validated. Fixed-speed energy-only transfer has no event stretching. Power-law frequency dependence and aging evidence have already been exposed and are training/diagnostic material.

## 2. Test distance against environment without circularity

**Exploratory empirical ansatz; known regression mathematics, not first-principles physics:**

\[
Y_i=aD_i+bV_i,\qquad V_i=\int_{\mathrm{path}\ i}W(\mathbf{x})\,ds.
\]

Y is ln(1+z_transfer). Separate source/observer Doppler and gravitational factors before interpreting total observed redshift as transfer, or marginalize their independently specified distributions. Retain negative observed redshifts rather than clipping them. a and b have inverse-length units; the physical domain must keep total local loss nonnegative.

- Distance-only model: b=0. Void model: a and b shared across targets. If the density definition or threshold has adjustable parameters, count and freeze those too.
- Choose W from independently mapped ordinary matter and survey completeness; do not define a void by a dark-halo inference, by the redshift residual being fitted, or by a vanishing acceleration that could merely reflect force cancellation.
- Begin with independent distance indicators (parallax, geometric masers, calibrated stellar/SBF indicators). Published distances may be stipulated, but preserve the brightness and calibration assumptions. Do not silently identify luminosity distance with static path length outside that explicit contract.
- D and V must have different proportions among the selected sightlines. If V=fD for every target, only a+bf is identifiable; more data on the same relation cannot measure the void effect. Check the weighted design-matrix singular values and profile uncertainty before interpreting any fitted b.
- Use target groups matched in distance but different void exposure, and groups matched in exposure but different distance. Include selection, shared distance zero points, correlated source motions and environment-map errors. A redshift-space map requires an explicit distortion/model sensitivity; it is not automatically an independent spatial void map.
- Do not add a separate velocity offset for every galaxy to force agreement. Fix its population model on training information and propagate uncertainty.

Deliverable: a provenance table of target ID, sky direction, independent distance and error, redshift and frame, ordinary-matter path estimate and uncertainty, calibration assumptions and prior exposure role; then distance-only/void model training and validation predictions. Presently no such independently mapped common path sample has been certified. The previously studied 164 distance groups, six masers and eight ELVES/ALFALFA objects are exposed; none can become a fresh holdout after retuning.

## 3. Predict spectra, clocks and brightness with the same parameters

- Propagate an emitted spectral distribution and an event envelope through the same interaction. Predict line centroids, widths, continuum shape, polarization and any scattering. Distinguish deterministic coherent shifts from stochastic jumps: an upper bound for one is not an upper bound for all.
- Derive what an atomic clock, ruler and detector measure. A coordinate time change applied equally to the signal and its measuring apparatus can cancel the proposed observable. Do not append clock immunity without a physical coupling prescription.
- **Conditional bookkeeping identity, not a new law:** in static Euclidean geometry, with isotropic source luminosity L, surviving photon fraction p, spectral stretch S_nu and event stretch S_t,

\[
F=\frac{L\,p}{4\pi D^2 S_\nu S_t}.
\]

  This is bolometric flux with no lensing or angular redistribution. Any scattering, absorption, beaming or geometric correction must come from the model. If both stretches equal S and photons survive, the factor is S^-2. Band flux requires the shifted spectrum and detector response. This brightness prediction can affect distance calibration; compare observer fluxes with explicit source/calibration nuisances instead of importing an incompatible cosmological distance modulus.
- Compute thermal-spectrum evolution with the candidate's phase-space transport and source history. Prior failed FIRAS histories remain failed but do not exclude every imaginable thermal origin. Do not retune an independent emissivity curve to erase each residual without counting it as a new model.
- Derive companion propagation/capture separately and close momentum as well as energy. A common light/GW delay does not alone prove ordinary gravitons can store the requested deposits. Galaxy force fits do not determine photon conversion while an independent gravity-amplitude degeneracy remains.

## 4. Freeze before genuinely withheld evaluation

Freeze candidate equations, all interaction and nuisance settings, source preparation, environment definition, distance calibration, quality cuts, exclusions, error/selection model, metrics and failure thresholds in a versioned commit. Specify untouched target IDs or a sealed selection/query and record a hash before retrieving outcomes. Do not use earlier validation outcomes again as new validation after repair.

Run the frozen pipeline once on the reserved observations. Preserve all predicted and observed rows, failures, outliers and sensitivity results. If a mechanism is revised, give it a new version and acquire a new reserved test; previous failures remain visible. Training fit quality is not the final success criterion. Success requires predictive uncertainty and joint agreement, not exact reproduction of noisy individual measurements.

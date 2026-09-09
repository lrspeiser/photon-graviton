from pathlib import Path
import os
import json,shutil
ROOT=Path(__file__).resolve().parents[3]
CONTRACT=json.loads((ROOT/'research_plan/universe-contract.json').read_text(encoding='utf-8'))
BASE=ROOT/'research_work'
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'data-audit'
registry=json.loads((ROOT/'research_work/results/specification/observable-registry.json').read_text(encoding='utf-8'))
# Measurement/proof targets, units, audited local source IDs, nuisance/calibration
# dependencies, and missing prerequisites for every roadmap requirement.
rows={
'R01':('Dimensionless source separation in material rods; clock ratios','dimensionless; s','', 'Matter metric, material rod response, atomic transition convention','Final operational contract'),
'R02':('Transition amplitudes and differential conversion rates','s^-1; m^-1; m^2','', 'Field normalization, backgrounds, occupation numbers','Derived interaction kernel and reverse channels'),
'R03':('Total energy and momentum balance including boundary flow','J; kg m s^-1; W','', 'Control volume, frame, gravitational work, interaction partition','Full on-shell stress-energy derivation'),
'R04':('Constraint count and kinetic/gradient eigenvalues','degrees of freedom; normalized eigenvalues','', 'Background solution, gauge, cutoff','Background-specific stability and radiative analysis'),
'R05':('Companion mass, spin, dispersion and polarizations','eV/c^2; omega(k) in s^-1','', 'Matter frame and field identity','Selected field identity and observational transfer function'),
'R06':('Observed multiband fluxes and inferred luminosity history','W m^-2; W; J','dustpedia_themis,dustpedia_dl14,dustpedia_schema', 'Distance, SED model, dust, band completeness, initial mass function','Underlying photometry and historical fuel/source model'),
'R07':('Spectral redshifts versus independent distance information','dimensionless; Mpc','cf4_catalog,cf4_schema,supernova_standardized', 'Reference frame, peculiar velocities, distance indicators','Raw spectral calibration and candidate distance reconstruction'),
'R08':('Transient widths and spectral-age progression','days; dimensionless stretch','transient_widths,spectral_aging', 'Rest-band mapping, templates, selection, source population','Raw multiepoch observations and joint covariance'),
'R09':('Flux/counts and standardized apparent brightness','W m^-2; counts s^-1; mag','supernova_standardized,supernova_covariance,supernova_schema', 'Absolute calibration, light-curve standardization, dust, selection','Candidate-consistent source calibration and count model'),
'R10':('Angular sizes and surface brightness','rad; W m^-2 sr^-1','sparc_catalog,sparc_curves', 'Physical size, inclination, distance reciprocity','Independent angular-size targets beyond processed galaxy tables'),
'R11':('Line width, image profile and polarization transfer','km s^-1; rad; dimensionless Stokes fractions','spectral_aging', 'Instrument line spread, intrinsic source width, scattering','Resolved spectra, PSFs, polarimetry and their covariances'),
'R12':('Optical/hyperfine/cavity frequency ratios and drifts','dimensionless; yr^-1','clock_paper_text', 'Spacer aging, temperature, atomic sensitivities, detrending','Actual frequency time series and secular-drift likelihood'),
'R13':('Local orbit, timing and gravitational deflection residuals','m; s; rad','', 'Ephemeris, source masses, reference-frame conventions','Primary local gravity/binary data and forward models'),
'R14':('Capture fraction, residence time and stored-state pressure','dimensionless; yr; Pa','', 'Cross-section, target distribution, recoil, reverse processes','Physical capture and support mechanism'),
'R15':('Halo radial density and predicted velocity field','kg m^-3; km s^-1','sparc_catalog,sparc_curves', 'Inclination, distance, baryon mass-to-light ratio','Spatial transport/self-gravity solution'),
'R16':('Galaxy rotation profiles and population scaling residuals','km s^-1; m s^-2','sparc_catalog,sparc_curves,sparc_split', 'Shared distance and inclination errors, stellar mass calibration','Unused validation sample and covariance-aware likelihood'),
'R17':('Dispersions, satellite/stream positions and disk vertical motion','km s^-1; rad','', 'Anisotropy, membership, geometry, selection','Primary non-rotation dynamical samples'),
'R18':('Shear, image positions and lens time delays','dimensionless; rad; days','', 'Source redshifts, shear calibration, mass-sheet degeneracy','Independent lensing measurements and same-field lens equations'),
'R19':('Cluster gas profiles, spectra, velocities and lensing','keV; counts; km s^-1; dimensionless shear','cluster_source_text,cluster_gas_source_text,cluster_derived', 'Hydrostatic bias, distance, R500 definition, gas calibration','Profile-level observations and lensing; published masses are derived'),
'R20':('Offsets and evolution in merging systems','rad; km s^-1','', 'Merger orientation/time, gas response, lens reconstruction','Merger images/spectra and dynamical simulations'),
'R21':('Dynamics conditioned on illumination and environment','km s^-1; W m^-2','sparc_catalog,dustpedia_themis,matched_energy_derived', 'Catalog completeness, selection, distance, environment definition','Independently selected environmental comparison sample'),
'R22':('Thermal spectra, stellar evolution and gas heating/cooling','K; W; abundance ratios','', 'Opacity, composition, transport, stellar calibration','Primary thermal/stellar constraints and evolution calculation'),
'R23':('GW strain, phase, polarization and EM arrival difference','dimensionless strain; s','', 'Source emission lag, detector calibration, binary model','Strain data and consistent emission/propagation model'),
'R24':('Compact-object spectra, orbits and ringdown','Hz; s; angular size','', 'Mass/spin inference, accretion, detector calibration','Compact-object observations and strong-field solution'),
'R25':('CMB monopole spectrum and distortion residuals','cm^-1; MJy sr^-1; kJy sr^-1','cmb_spectrum', 'Absolute calibration, Galactic subtraction, supplied blackbody reconstruction','Full error covariance and radiation-production history'),
'R26':('CMB angular TT/EE/TE spectra and correlations','multipole ell; microkelvin^2','cmb_TT,cmb_EE', 'Foregrounds, masks, beam, bin windows, calibration','Full likelihood/TE and candidate perturbation transfer functions'),
'R27':('Angular/redshift clustering and correlation features','dimensionless redshift; rad; correlation amplitude','', 'Survey selection, fiducial coordinate conversion, tracer bias','Galaxy catalogs and candidate-consistent clustering pipeline'),
'R28':('Abundances, stellar populations and chronology indicators','abundance ratios; flux; model-dependent years','', 'Nuclear rates, stellar tracks, dust, distance','Independent abundance/population inputs and new-history model'),
'R29':('Static background solution, redshift drift, topology signals','yr^-1; curvature length; angle','', 'Clock/rod convention, boundary conditions, observer motion','Solved background, stability and primary drift/topology tests'),
'R30':('Integrated radiation backgrounds, fuel and entropy budget','W m^-2 sr^-1; J; J K^-1','cmb_spectrum,dustpedia_themis', 'Obscured populations, source history, absorption, completeness','Multi-band background compilation and closed global ledger'),
'R31':('Parameter identifiability and selection/calibration effects','dimensionless likelihood; parameter covariance','sparc_split,supernova_covariance', 'Shared nuisance structure and reused samples','Candidate-specific joint likelihood and simulation calibration'),
'R32':('Previously frozen discriminating predictions','units of chosen observable','', 'Multiplicity, exposure history, model flexibility','Unexposed validation targets and frozen predictive distribution'),
}
assert set(rows)=={r['id'] for r in registry}
audit=json.loads((OUT/'input-audit.json').read_text(encoding='utf-8'))
ids={s['id'] for s in audit['inventory']}
for r in registry:
    target,units,sources,nuisance,missing=rows[r['id']]
    source_ids=sources.split(',') if sources else []
    assert set(source_ids)<=ids
    r.update(target=target,units=units,local_input_ids=source_ids,nuisance_and_model_dependencies=nuisance,missing_for_completion=missing,availability='partial local products; not complete measurement pipeline' if source_ids else 'no primary observational input assigned; theoretical proof or acquisition required',status='mapped; completion unproven',validation_exposure='all assigned recovered products exploratory; no new blind data')
for r in registry:
    r['universe_contract_file']='research_plan/universe-contract.json'
    r['distance_policy']='Adopted published galaxy distances fixed; provenance retained, not fitted or re-inferred'
    r['excluded_premises_policy']='No dark-matter source, expansion or Big-Bang assumption; retain observations and independently derive permitted interpretation'
    if r['id']=='R07':r['missing_for_completion']='Source/detector spectral calibration and a redshift law at fixed published galaxy distances'
    if r['id']=='R28':r['missing_for_completion']='Abundance/population inputs and chronology without Big-Bang or expansion-derived initial conditions'
(OUT/'observable-registry.json').write_text(json.dumps(registry,indent=2),encoding='utf-8')
protocol=dict(schema_version=1,status='Frozen process rules; candidate-specific likelihood and validation sample remain unregistered',scope='Applies prospectively; does not relabel earlier fits as preregistered',observational_contract='Published galaxy distances are fixed fictional facts; no assumed dark matter, expansion or Big Bang premises.',existing_data='All previously recovered/used data exploratory, including old train/validation/test partitions',prohibited_claims=['Do not count Planck BestFit columns as measured data','Do not treat hydrostatic or rotation-inferred mass as a gravity-independent observation','Do not interpret a new split of exposed galaxies as blind validation','Do not transfer a calibration, luminosity or inferred radius unchanged after changing its distance/matter law'],before_new_model_selection=['Freeze action, finite parameter set, priors and backgrounds','Freeze source history and capture/response parameter sharing','Freeze sample inclusion and exclusions by observable quality, not residual outcome','Record data/code hashes and all transformations','Define primary likelihood, nuisance covariance, comparator and held-out predictive metric','Declare optimization/evaluation budget and stopping rule before fitting','Choose genuinely unexposed data with overlap checks or label the next run exploratory'],comparison_rules=['Report physical absolute-normalization predictions separately from fitted amplitudes','Report galaxy-level and point-level summaries with shared nuisance parameters','Use full relevant covariance when available; disclose diagonal approximations','Record invalid or unstable parameter regions; do not silently discard their failures','A candidate revision gets a new version and does not retain blind status on exposed targets'],numerical_vs_observational='Numerical tolerances validate computations only. No observational acceptance cutoff is invented before the candidate and likelihood exist.',holdout=dict(selected=False,reason='No genuinely unused validation sample established; data acquisition and overlap audit required'),completion=False)
protocol['universe_contract_file']='research_plan/universe-contract.json'
protocol['excluded_active_premises']=CONTRACT['excluded_active_premises']
protocol['prohibited_claims'] += ['Do not refit adopted galaxy distances or reject them because of their original Hubble-flow derivation.','Do not use assumed dark-matter halos, expansion-derived chronology or Big-Bang standard rulers as physical inputs.']
(OUT/'validation-protocol.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')
checks=dict(requirements=len(registry),unique_requirement_ids=len({r['id'] for r in registry}),all_requirements_have_units_target_owner_gap=all(r['units'] and r['target'] and r['tasks'] and r['missing_for_completion'] for r in registry),all_source_references_resolve=True,local_input_products=len(ids),all_assigned_inputs_exist=audit['all_listed_paths_present'],prospective_blind_sample_established=False)
(OUT/'registry-checks.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
print(json.dumps(checks,indent=2))

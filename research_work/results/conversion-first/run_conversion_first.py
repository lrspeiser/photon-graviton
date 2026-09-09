"""Explore the conversion-first benchmark without modifying recovered evidence.

Run with PHOTON_GRAVITON_RESULTS pointing to a fresh output parent.
This is an empirical transport benchmark, not a microscopic theory.
"""
from pathlib import Path
import csv
import hashlib
import json
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
OUT = Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT / 'research_work/generated')) / 'conversion-first'
C = 299792.458  # km/s
MPC_PER_MLY = 0.3066013938
ARCHIVE_ALPHA = 7.7315e-5 / MPC_PER_MLY


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8', newline='\n')


def predicted(k, distance, linear=False):
    x = k * distance / C
    return C * (x if linear else np.expm1(x))


def fit(distance, observed, linear=False):
    result = minimize_scalar(lambda k: float(np.sum((predicted(k, distance, linear) - observed)**2)),
                             bounds=(0., 150.), method='bounded', options={'xatol': 1e-9})
    if not result.success:
        raise RuntimeError(result.message)
    # Include the physical boundaries, which bounded optimization may not visit.
    candidates = [0., float(result.x), 150.]
    return min(candidates, key=lambda k: np.sum((predicted(k, distance, linear) - observed)**2))


def metrics(residual):
    return dict(count=len(residual), rmse_km_s=float(np.sqrt(np.mean(residual**2))),
                mae_km_s=float(np.mean(abs(residual))), bias_km_s=float(np.mean(residual)),
                nominal_95pct_residual_coverage=float(np.mean(abs(residual) <= 1.96*300)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'redshift_paper/all_164_groups.csv'
    with source.open(encoding='utf-8', newline='') as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 164
    features_dir = ROOT / 'redshift_paper_sources/time-redshift-expanded'
    features = {str(r['pgc']): r for r in json.loads((features_dir / 'features.json').read_text())}
    labels = {}
    for partition in ['train', 'validation', 'test']:
        labels.update(json.loads((features_dir / (partition + '_labels.json')).read_text()))
    for row in rows:
        original = features[row['pgc']]
        assert float(row['catalog_distance_mpc']) == original['distance_mpc']
        assert row['split'] == original['split']
        assert float(row['observed_cmb_cz_kms']) == labels[row['pgc']]
        assert abs(float(row['observed_cmb_z'])*C - labels[row['pgc']]) < 1e-9
    assert len({r['group_pgc'] for r in rows}) == len(rows)
    distance = np.array([float(r['catalog_distance_mpc']) for r in rows])
    observed = np.array([float(r['observed_cmb_cz_kms']) for r in rows])
    tiles = np.array([int(r['sky_tile']) for r in rows])
    splits = np.array([r['split'] for r in rows])
    for tile in np.unique(tiles):
        assert len(set(splits[tiles == tile])) == 1
    train = splits == 'train'
    k_exp = fit(distance[train], observed[train])
    k_linear = fit(distance[train], observed[train], True)
    assert 0 < k_exp < 150 and 0 < k_linear < 150
    models = {
        'conversion_exponential': predicted(k_exp, distance),
        'linear_small_depth_control': predicted(k_linear, distance, True),
        'archived_conversion_rate': predicted(C*ARCHIVE_ALPHA, distance),
        'zero_conversion_diagnostic': np.zeros_like(distance),
    }
    rng = np.random.default_rng(2026090901)
    training_tiles = np.unique(tiles[train])
    boot = []
    for _ in range(500):
        ids = np.concatenate([np.flatnonzero(train & (tiles == tile))
                              for tile in rng.choice(training_tiles, len(training_tiles), replace=True)])
        boot.append(fit(distance[ids], observed[ids]))
    stats = {part: {name: metrics((value-observed)[splits == part]) for name, value in models.items()}
             for part in ['train', 'validation', 'test']}
    rng = np.random.default_rng(2026090902)
    test = splits == 'test'
    test_tiles = np.unique(tiles[test])
    differences = []
    for _ in range(1000):
        ids = np.concatenate([np.flatnonzero(test & (tiles == tile))
                              for tile in rng.choice(test_tiles, len(test_tiles), replace=True)])
        differences.append(metrics((models['conversion_exponential']-observed)[ids])['rmse_km_s'] -
                           metrics((models['linear_small_depth_control']-observed)[ids])['rmse_km_s'])
    alpha = k_exp/C
    fitted_z = models['conversion_exponential']/C
    transfer = -np.expm1(-alpha*distance)
    assert np.max(abs(transfer - fitted_z/(1+fitted_z))) < 1e-14
    assert np.max(abs(transfer + np.exp(-alpha*distance)-1)) < 1e-14
    result = dict(scope='Exploratory reused galaxy data, fixed published distances, stipulated conversion law',
                  protocol_sha256=digest(HERE/'protocol.json'), source_sha256=digest(source),
                  source_files={str(p.relative_to(ROOT)): digest(p) for p in
                                [source, features_dir/'features.json', features_dir/'train_labels.json',
                                 features_dir/'validation_labels.json', features_dir/'test_labels.json']},
                  fitted=dict(alpha_per_mpc=alpha, alpha_per_mly=alpha*MPC_PER_MLY,
                              c_alpha_km_s_per_mpc=k_exp, linear_c_alpha=k_linear,
                              training_tile_bootstrap_95_c_alpha=np.percentile(boot, [2.5,97.5]).tolist(),
                              boundary_hit=False),
                  distance_range_mpc=[float(distance.min()), float(distance.max())],
                  fitted_conversion_fraction_range=[float(transfer.min()), float(transfer.max())],
                  statistics=stats,
                  paired_test_rmse_exponential_minus_linear_95=np.percentile(differences,[2.5,97.5]).tolist(),
                  distance_refitted=False, new_independent_validation=False,
                  imported_expansion_or_dark_matter_model=False,
                  physical_cause_identified=False)
    diagnostic = []
    for z in [.01, .1, 1., 2.]:
        shift = 1+z
        width_limit = 10/C  # illustrative added 10 km/s energy-width scale
        epsilon_max = np.log1p(width_limit**2)/np.log(shift)
        diagnostic.append(dict(conversion_z=z, supplied_energy_fraction=z/shift,
                               illustrative_distance_mpc=np.log(shift)/alpha,
                               event_duration_ratio_fixed_conditions=1.,
                               illustrative_comparison_duration_ratio=shift,
                               flux_relative_to_no_loss_static_source=1/shift,
                               flux_if_additional_duration_stretch_S=1/shift**2,
                               maximum_event_fraction_for_illustrative_10kms_width=float(epsilon_max),
                               minimum_mean_event_count=float(np.log(shift)/epsilon_max)))
    # Independent ODE vs matrix exponential for a closed four-reservoir model.
    # Capture-rate and lifetime values are illustrative, in arbitrary common time units.
    ledger_checks = []
    for h in [.01, 1., 10.]:
        for capture in [.01, 1., 10.]:
            for decay in [0., .1]:
                matrix = np.array([[-h,0,0,0], [h,-capture,0,0],
                                   [0,capture,-decay,0], [0,0,decay,0]])
                initial = np.array([1.,0,0,0])
                sol = solve_ivp(lambda t,y: matrix@y, (0.,5.), initial, rtol=1e-10, atol=1e-12)
                exact = expm(5*matrix)@initial
                error = float(np.max(abs(sol.y[:,-1]-exact)))
                residual = float(np.max(abs(sol.y.sum(axis=0)-1)))
                assert sol.success and error < 1e-9 and residual < 1e-12 and np.min(sol.y)>-1e-11
                ledger_checks.append(dict(h=h,capture=capture,decay=decay,
                                          ode_matrix_error=error,energy_residual=residual))
    # Finite storage with constant incoming captured power: exact integral.
    retention = [dict(age_over_lifetime=x, fraction_of_captured_energy_still_stored=float(-np.expm1(-x)/x))
                 for x in [.01, .1, 1., 10., 100.]]
    capture = [dict(path_optical_depth=x,capture_probability=float(-np.expm1(-x))) for x in [.01,.1,1.,3.,10.]]
    # A time-stationary delay distribution broadens but preserves separation of centroids.
    delay = np.linspace(0,100,10001)
    weights = np.exp(-delay/10)
    weights /= weights.sum()
    separations = []
    for separation in [1.,10.,100.]:
        first = float(weights@delay)
        second = float(weights@(delay+separation))
        ratio = (second-first)/separation
        assert abs(ratio-1) < 1e-12
        separations.append(dict(emission_separation=separation, centroid_separation_ratio=ratio))
    # Dimensionless massless three-body energy/momentum kinematics, c=E_initial=1.
    kinematics = []
    for remaining in [.1,.5,.9]:
        for theta in [0., .01, 1.]:
            required_p2 = (1-remaining)**2 + 4*remaining*np.sin(theta/2)**2
            one_companion_p2 = (1-remaining)**2
            mismatch = float(required_p2-one_companion_p2)
            assert (abs(mismatch)<1e-15) if theta==0 else mismatch>0
            kinematics.append(dict(remaining_energy=remaining, photon_deflection_radians=theta,
                                   momentum_squared_excess_over_single_massless_companion=mismatch))
    result.update(illustrative_predictions=diagnostic,capture_scenarios=capture,retention_scenarios=retention)
    # Closed resonant three-mode toy restricted to |one high photon> and
    # |one lower photon, one companion>. This does not assign a graviton spin/action.
    toy_cases = []
    for omega_companion in [.01, .1, .5]:
        omega_initial = 1.
        omega_final = omega_initial-omega_companion
        coupling = .2
        hamiltonian = np.array([[omega_initial,coupling],[coupling,omega_final+omega_companion]])
        for time in [0.,1.,4.,8.]:
            state = expm(-1j*hamiltonian*time)@np.array([1.,0.],dtype=complex)
            probability = float(abs(state[1])**2)
            assert abs(probability-np.sin(coupling*time)**2) < 1e-13
            photon_energy = abs(state[0])**2*omega_initial+probability*omega_final
            companion_energy = probability*omega_companion
            interaction_energy = float(np.real(np.vdot(state,(hamiltonian-np.eye(2))@state)))
            assert abs(photon_energy+companion_energy+interaction_energy-1) < 1e-13
            assert abs(np.vdot(state,state)-1) < 1e-13
            toy_cases.append(dict(companion_frequency=omega_companion,time=time,
                                  conversion_probability=probability,photon_energy=float(photon_energy),
                                  companion_energy=companion_energy,interaction_energy=interaction_energy))
    result['resonant_mode_toy'] = dict(cases=toy_cases,
        status='Closed energy-conserving inelastic mode toy; no microscopic graviton vertex, physical rate or irreversible drift derived')
    save('results.json', result)
    checks = dict(source_rows_match_original_features_and_labels=True, group_tile_integrity=True,
                  no_distance_changes=True, energy_redshift_identity=True, ledger_cases=ledger_checks,
                  stationary_delay_centroid_cases=separations, massless_kinematic_cases=kinematics,
                  resonant_three_mode_energy_cases=len(toy_cases),
                  scope='Mathematical and input-integrity checks, not observational verification of microscopic physics')
    save('checks.json', checks)
    with (OUT/'predictions.csv').open('w',newline='',encoding='utf-8') as stream:
        fieldnames=['pgc','group_pgc','split','sky_tile','catalog_distance_mpc','observed_cmb_z',
                    'predicted_conversion_z','residual_km_s','companion_energy_fraction_if_full_transfer']
        writer=csv.DictWriter(stream,fieldnames=fieldnames,lineterminator='\n'); writer.writeheader()
        for i,row in enumerate(rows):
            item={key:row[key] for key in fieldnames[:6]}
            item.update(predicted_conversion_z=fitted_z[i],residual_km_s=models['conversion_exponential'][i]-observed[i],
                        companion_energy_fraction_if_full_transfer=transfer[i])
            writer.writerow(item)
    registry=json.loads((ROOT/'research_work/results/data-audit/observable-registry.json').read_text())
    preliminary={
        'R01':'Conversion-only baseline uses fixed material standards; complete matter law still missing.',
        'R02':'Massless momentum diagnostic and resonant three-mode Hamiltonian toy; no physical rate or field action.',
        'R03':'18 closed four-reservoir ODE checks against independent matrix exponential.',
        'R06':'Per-photon production fractions calculated; no integrated stellar source history derived.',
        'R07':'New constant-loss fit to all 164 previously exposed groups with inherited train/test split.',
        'R08':'Fixed conditions and stationary delay model preserve event separation; duration explanation missing.',
        'R09':'Conditional static flux scenarios; no brightness data fit.',
        'R11':'Existing Poisson width requirement evaluated at four illustrative shifts; no measured line-profile fit.',
        'R14':'Capture optical-depth and finite-retention scenarios; no microscopic cross section or supported state.',
        'R15':'Prior point-source radial shape calculation exists; realistic source and normalization unresolved.',
        'R16':'Prior rotation diagnostic exists; this pass does not supply an absolute gravity prediction.',
        'R18':'Prior motion/lensing mismatch remains; no newly derived field response.',
        'R21':'Recovered sample has no path environment fields; seeking/focusing cannot be inferred from this fit.',
        'R30':'Explicit decay receiver included; no global steady-state radiation/entropy solution.',
        'R31':'Tile bootstrap and exponential/linear comparison; shared calibration/motion degeneracies remain.',
        'R32':'New frozen protocol on reused data; genuinely new validation remains missing.',
    }
    coverage=[dict(id=r['id'],area=r['area'],tasks=r['tasks'],complete=False,
                   this_pass=preliminary.get(r['id'],'No new completed prediction in this pass.'),
                   next_requirement=r['required_response'], original_registry_gap=r['missing_for_completion'])
              for r in registry]
    assert len(coverage)==32 and len({r['id'] for r in coverage})==32
    save('coverage.json',coverage)
    lines=['# Conversion-first exploratory research pass','',
           'The user permits conversion without special void/time stretching as its cause. The original time idea was intended to cause light stretching, including for nearby galaxies; it did not remove local redshift. This pass does not establish a microscopic mechanism or complete the research goal.','',
           '## Galaxy redshift result','',
           f'Fitted alpha = {alpha:.10g} per Mpc = {alpha*MPC_PER_MLY:.10g} per million light-years. The equivalent c*alpha is {k_exp:.4f} km/s/Mpc (a conversion coefficient, not an assumed expansion rate).',
           f'Training sky-tile bootstrap 95% interval for c*alpha: {np.percentile(boot,2.5):.3f} to {np.percentile(boot,97.5):.3f} km/s/Mpc.',
           f'All 164 fixed published distances are unchanged, spanning {distance.min():.2f} to {distance.max():.2f} Mpc. The primary error model uses a fixed illustrative 300 km/s residual scale. Distance uncertainties are retained as provenance, not fitted.',
           '', '| Partition | Count | Exponential RMSE km/s | Linear control RMSE km/s | Archived-rate RMSE km/s | Exponential nominal 95% coverage |',
           '|---|---:|---:|---:|---:|---:|']
    for part,values in stats.items():
        m=values['conversion_exponential']
        lines.append(f"| {part} (previously exposed) | {m['count']} | {m['rmse_km_s']:.3f} | {values['linear_small_depth_control']['rmse_km_s']:.3f} | {values['archived_conversion_rate']['rmse_km_s']:.3f} | {m['nominal_95pct_residual_coverage']:.3f} |")
    ci=result['paired_test_rmse_exponential_minus_linear_95']
    lines += ['',f'Paired test-tile bootstrap interval for exponential minus linear RMSE: {ci[0]:.3f} to {ci[1]:.3f} km/s. Negative favors the exponential benchmark. This interval is descriptive on reused data, not a new blind test.',
              'The exponential law approaches the linear control at these short distances. The test difference interval includes zero, so this diagnostic does not resolve their curvature. Nominal 95% coverage is only 80% on the 25 reused test objects; the fixed 300 km/s error model is not a demonstrated adequate description. A distance-redshift trend cannot identify companions, their identity, or a first-principles cause. Catalog frame corrections, shared distance calibrations, motions and selection are not fully modeled. No rotation, lensing, spectrum or timing observation was fitted here.',
              '', '## Predictions from the same rate','',
              '| Conversion z | Energy transferred | Fixed-condition duration ratio | Static flux / unshifted-source flux |', '|---|---:|---:|---:|']
    for d in diagnostic:
        lines.append(f"| {d['conversion_z']} | {100*d['supplied_energy_fraction']:.3f}% | 1 | {d['flux_relative_to_no_loss_static_source']:.6f} |")
    lines += ['', 'These are scenarios, including high-z extrapolations beyond the fitted nearby sample. The flux formula assumes unchanged photon count, isotropic emission, Euclidean geometry, unchanged clocks and fixed speed. A separate duration stretch by S would reduce flux by another factor S. It cannot be inserted without deriving a timing mechanism.',
              '', '## Capture, retention and energy','',
              'All 18 ODE/matrix-exponential ledger comparisons passed, including decay into an explicit receiving reservoir. They establish consistent bookkeeping for chosen rates, not a microscopic interaction.',
              'Capture probability is 1-exp(-optical_depth). An optical depth of 3 gives about 95% capture, but neither the optical depth nor a capture cross section is derived. If captured power is constant for age T and storage lifetime is t_d, the stored fraction of all captured energy is (1-exp(-T/t_d))/(T/t_d). At T=t_d it is 63.2%; at T=10 t_d it is about 10%.',
              'Conservative gravitational focusing can bend paths toward a well. It is not an absorption law and does not by itself establish a permanent deposited reservoir. A claim that companions seek larger wells needs a ray/transport law and a defined absorbing or bound state. Ordinary freely propagating radiation and an added attractive interaction must not be silently interchanged.',
              '', '## Microscopic cause and timing requirements','',
              'A time-independent resonant three-mode Hamiltonian can exchange one higher-frequency photon for a lower-frequency photon plus a companion without special clock behavior. Twelve explicit matrix-exponential checks conserve energy. This is a limited quantum-mode toy: its conversion is reversible/oscillatory, gives a line mixture, and does not establish a nonzero physical photon-graviton coupling or the observed rate. See inelastic-mode-derivation.md for the derivation and missing steps.',
              'For one initial massless photon, one lower-energy massless photon and one massless companion in otherwise empty space, energy/momentum balance gives |p_initial-p_final|^2 = (Delta E/c)^2 + 2 E_initial E_final (1-cos(theta))/c^2. A single companion with energy Delta E has squared momentum (Delta E/c)^2. Thus this simple channel requires collinear photon emission. This does not establish a nonzero transition rate; a background, recoil target or different dispersion needs separate accounting.',
              'A time-stationary delay kernel has I_out(t)=integral K_R(u) I_in(t-u) du. Shifting an input pulse shifts its output by the same amount. Two otherwise identical input pulses therefore retain their centroid separation, even if each is broadened. The numerical delay checks confirm that statement for an illustrative exponential delay distribution. Constant-speed drift is its zero-width limit.',
              'Consequently, adding the same random travel delays cannot generally stretch every feature of an arbitrary light curve by 1+z. A completed conversion-only explanation needs another derived time-transfer effect, source-population explanation, or explicitly time-dependent/nonlinear interaction. This is not a requirement to restore special void clocks, nor a proof against every conversion mechanism.',
              'Transient duration and spectral-aging observations remain targets: see [Blondin et al.](https://arxiv.org/abs/0804.3595) and [DES time-dilation analysis](https://arxiv.org/abs/2406.05050). Their expansion interpretation is not adopted here; their measurement reductions must be audited before a candidate-specific likelihood. No time-dilation data fit was performed in this pass.',
              '', '## Work-through of the eight checklist stages','',
              '1. Redshift forward model: specified in protocol.json, with fixed distances and constant-rate benchmark.',
              '2. Spectra/counts/energy/timing: conditional predictions and mathematical checks computed; physical interaction missing.',
              '3. Training fit: executed on recovered training rows with frozen protocol and descriptive bootstrap.',
              '4. New validation: not available; historical validation/test rows evaluated and explicitly labeled reused.',
              '5. Source map: per-photon conversion computed; joint luminosity histories and corrected environments missing.',
              '6. Capture/storage: rate/lifetime scenarios checked; spatial cross sections and supported state missing.',
              '7. Motion/lensing: prior conditional diagnostic retained; no unified response derived in this pass.',
              '8. Full register: all 32 requirements mapped in coverage.json; none newly claimed complete at full scope.',
              '', '## Next actionable work','',
              'Construct a candidate inelastic interaction with its receiving/recoil sector, and derive its joint energy, direction, polarization and arrival-time kernel. Use the fitted loss coefficient as a target rather than a derived constant. In parallel in the research sequence, obtain unexposed spectra and transient time-series with measurement provenance. A supported deposit law and joint gravity/light-propagation action are needed before fitting extra gravity. No assumption of expansion, Big Bang or an independent dark-matter source is required by this work plan.',
              '', 'Inputs, predictions, metrics and checks are in results.json, predictions.csv, checks.json and coverage.json alongside this report. The runnable source and protocol are preserved in the repository.']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(fitted=result['fitted'],statistics=stats,
                         paired_difference_95=ci,ledger_cases=len(ledger_checks),coverage_areas=len(coverage)),indent=2))


if __name__ == '__main__':
    main()

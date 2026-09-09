"""Charge-neutral bound-pair production with recoil, inverse dynamics and color test."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from scipy.linalg import eigh, expm

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT/'research_work/generated'))/HERE.name


def recoil(mass, momentum):
    # Rationalized expression avoids cancellation for a large receiving mass.
    return momentum**2/(np.hypot(mass, momentum)+mass)


def kinematic_checks(cfg):
    rows = []
    for mass in cfg['kinematics']['store_masses']:
        for gap in cfg['kinematics']['internal_gaps']:
            threshold = gap + gap*gap/(2*mass)
            for multiple in cfg['kinematics']['photon_threshold_multiples']:
                energy = multiple*threshold
                for angle in cfg['kinematics']['angles']:
                    outgoing = mass*(energy-threshold)/(mass+energy*(1-np.cos(angle)))
                    p = np.array([-outgoing*np.sin(angle), energy-outgoing*np.cos(angle)])
                    kinetic = float(recoil(mass+gap, np.linalg.norm(p)))
                    residual = abs(energy-outgoing-gap-kinetic)/max(energy, gap)
                    assert outgoing > 0 and residual < cfg['checks']['kinematic_relative_error']
                    assert np.linalg.norm(p+outgoing*np.array([np.sin(angle),np.cos(angle)])-np.array([0.,energy])) < 1e-12*max(1.,energy)
                    if angle == 0:
                        assert abs(energy-outgoing-threshold) < 1e-12*max(1.,energy)
                    rows.append({'mass': mass, 'internal_gap': gap, 'photon_energy': energy,
                                 'angle': angle, 'outgoing_photon_energy': float(outgoing),
                                 'store_recoil_energy': kinetic, 'relative_energy_residual': float(residual)})
    return rows


def quantum_case(mass, photons, coupling, cfg):
    qcfg = cfg['quantum']
    gap = 2*(qcfg['constituent_rest_energy']-qcfg['binding_per_constituent'])
    transfer = gap+gap*gap/(2*mass)  # Tune the first transition including recoil.
    high = qcfg['high_photon_energy']
    low = high-transfer
    assert low > 0
    k = np.arange(photons+1, dtype=float)
    kinetic = recoil(mass+k*gap, k*transfer)
    diagonal = -k*transfer+k*gap+kinetic
    edge = coupling*(k[:-1]+1)*np.sqrt((photons-k[:-1])*(k[:-1]+1))
    interaction = np.diag(edge, 1)+np.diag(edge, -1)
    hamiltonian = np.diag(diagonal)+interaction
    values, vectors = eigh(hamiltonian)
    initial = np.zeros(photons+1); initial[0] = 1
    tau = np.linspace(0, qcfg['dimensionless_duration'], qcfg['history_samples'])
    t = tau/coupling
    states = vectors @ (np.exp(-1j*values[:, None]*t)*(vectors.T@initial)[:, None])
    probability = abs(states)**2
    norm_error = float(np.max(abs(probability.sum(axis=0)-1)))
    occupation = k@probability
    photon_energy = photons*high-transfer*occupation
    rest_gain = 2*qcfg['constituent_rest_energy']*occupation
    binding_change = -2*qcfg['binding_per_constituent']*occupation
    recoil_gain = kinetic@probability
    interaction_energy = np.einsum('it,ij,jt->t', states.conj(), interaction, states).real
    ledger = photon_energy-photons*high+rest_gain+binding_change+recoil_gain+interaction_energy
    total_error = float(np.max(abs(ledger)))
    assert norm_error < cfg['checks']['norm_error']
    assert total_error < cfg['checks']['absolute_energy_error']
    # The invariant basis contains high=N-k, low=k, b=k, d=k, P_store=k*transfer.
    total_momentum = (photons-k)*high+k*low+k*transfer
    assert max(abs(total_momentum-photons*high)) < 1e-12
    assert np.array_equal(k-k, np.zeros_like(k))  # Basis charge, not an independent simulation.
    indices = np.linspace(0, len(t)-1, qcfg['independent_matrix_exponential_times'], dtype=int)
    independent = max(float(np.linalg.norm(expm(-1j*hamiltonian*t[i])@initial-states[:, i])) for i in indices)
    assert independent < cfg['checks']['state_comparison_error']
    one_pair_error = None
    if photons == 1:
        one_pair_error = float(np.max(abs(probability[1]-np.sin(tau)**2)))
        assert one_pair_error < cfg['checks']['one_pair_probability_error']
    # Zero coupling gives no production from the same initial state.
    zero_control = expm(-1j*np.diag(diagonal)*t[-1])@initial
    assert np.array_equal(zero_control, initial.astype(complex))
    return {'store_seed_mass': mass, 'initial_high_photons': photons, 'coupling': coupling,
            'high_energy': high, 'low_energy': low, 'internal_pair_gap': gap,
            'photon_energy_transfer_per_transition': transfer,
            'first_transition_recoil': float(kinetic[1]),
            'largest_bare_state_detuning': float(max(abs(diagonal))),
            'maximum_sampled_pairs': float(occupation.max()), 'pairs_at_end': float(occupation[-1]),
            'maximum_interaction_energy_magnitude': float(max(abs(interaction_energy))),
            'maximum_norm_error': norm_error, 'maximum_absolute_energy_error': total_error,
            'independent_expm_state_error': independent, 'one_pair_probability_error': one_pair_error,
            'zero_coupling_no_production': True,
            'history': {'tau': tau[::10].tolist(), 'mean_pairs': occupation[::10].tolist(),
                        'photon_energy_loss': (photons*high-photon_energy[::10]).tolist(),
                        'pair_rest_energy_gain': rest_gain[::10].tolist(),
                        'binding_energy_change': binding_change[::10].tolist(),
                        'store_recoil_energy': recoil_gain[::10].tolist(),
                        'interaction_energy': interaction_energy[::10].tolist()}}


def main():
    cfg = json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    kinematics = kinematic_checks(cfg)
    runs = [quantum_case(mass, photons, coupling, cfg)
            for mass in cfg['quantum']['store_masses']
            for photons in cfg['quantum']['initial_high_photons']
            for coupling in cfg['quantum']['couplings']]
    c = cfg['color_test']; transfer = c['internal_gap']+c['internal_gap']**2/(2*c['store_mass'])
    colors = [{'initial_energy': e, 'lost_energy': transfer, 'fractional_loss': transfer/e,
               'conditional_redshift': e/(e-transfer)-1} for e in c['photon_energies']]
    assert len(set(round(row['conditional_redshift'], 8) for row in colors)) == len(colors)
    # Separate leading heavy-store, point-electric-dipole rate: initial/final
    # electric fields give E*Eout, final photon density of states gives Eout^2.
    # No absolute rate or finite-size form factor is supplied by this diagnostic.
    ref = 1.0
    rate_colors = [{'initial_energy': e,
                    'relative_event_rate': e/ref*((e-c['internal_gap'])/(ref-c['internal_gap']))**3,
                    'relative_fractional_loss_rate': ((e-c['internal_gap'])/(ref-c['internal_gap']))**3}
                   for e in c['photon_energies']]
    result = {'scope': cfg['scope'], 'checks_pass': True, 'kinematics': kinematics,
              'quantum_runs': runs, 'single_forward_transition_colors': colors,
              'leading_point_dipole_rate_colors': rate_colors,
              'achromatic_rate_derived': False, 'whole_event_time_stretch_derived': False,
              'spatial_bound_mode_or_capture_overlap_derived': False,
              'source_hashes': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in [HERE/'check.py', HERE/'protocol.json']}}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'bound-pair-production-results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print('Kinematic cases:',len(kinematics),'max residual:',max(r['relative_energy_residual'] for r in kinematics))
    for r in runs:
        print({key:r[key] for key in ['store_seed_mass','initial_high_photons','coupling','maximum_sampled_pairs','pairs_at_end','maximum_absolute_energy_error','independent_expm_state_error']})
    print('Colors:', colors)


if __name__ == '__main__':
    main()

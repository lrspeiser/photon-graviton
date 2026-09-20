"""TF-1 stage 1: the exact-property gates G1, G2, G5, G6, G7, G8, G9 with their controls."""
import numpy as np
from common import save, read
import transverse as T


def run(out):
    field = T.load_field()
    axis, phi = field['axis'], field['phi']
    arch = read(T.ARCHIVED_2D)['cases']['primary-256']
    # live index rays reproduce the archived angles (the field is the archived one)
    index = T.index_rays_2d(axis, phi)
    arch_angles = np.array([r['angle'] for r in arch['rays']])
    live_angles = np.array([r['angle'] for r in index])
    repro = float(np.max(np.abs(live_angles - arch_angles))/np.max(np.abs(arch_angles)))
    steer = T.steering_rays_2d(axis, phi)
    steer_ref = T.steering_rays_2d(axis, phi, refined=True)
    straight = T.steering_rays_2d(axis, phi, zero=True)
    cmp = T.compare_rays(steer, index)
    scale = cmp['max_bend_steering']
    g5 = dict(**cmp, born_relative=float(np.max([abs(r['angle'] - r['born_angle']) for r in steer])/scale),
              step_refinement_relative=float(max(abs(a['angle'] - b['angle']) for a, b in zip(steer, steer_ref))/scale),
              zero_field_max_bend=float(max(abs(r['angle']) for r in straight)),
              archived_mirror_relative=arch['bend']['mirror_relative'], weak_regime=float(T.G_INDEX*field['max_abs_phi']),
              rays=steer, index_rays_live=index)
    probes = T.probes_on_field(axis, phi)
    g1 = dict(steering_speed_deviation=0., steering_speed_note='direction integrated as an angle: unit speed by construction',
              index_control_speed_deviation=T.index_speed_control(axis, phi))
    g6 = {str(a): T.lagrangian_family(a) for a in (0., .5, 1., 2.)}
    g7 = T.velocity_averages()
    g8 = T.focusing_geometry()
    g1['vector_speed_drift'] = g8['max_speed_drift']
    g9 = dict(orbits=T.circular_and_radial(), ensemble=T.ensemble())
    gates = dict(
        G1_work_free=g1['steering_speed_deviation'] < 1e-9 and g1['vector_speed_drift'] < 1e-9 and g1['index_control_speed_deviation'] > 1e-3,
        G2_composition=probes['steering_difference'] < 1e-14 and probes['vector_difference'] < 1e-14 and abs(probes['cwc1_ratio'] - 2) < 1e-9,
        G5_field_reproduction=repro < 1e-9,
        G5_inward=g5['all_inward'],
        G5_matches_index_bending=g5['relative_difference'] < .02 and g5['weak_regime'] < .02,
        G5_step_refinement=g5['step_refinement_relative'] < 1e-4,
        G5_straight=g5['zero_field_max_bend'] == 0.,
        G6_euler_lagrange=all(v['euler_lagrange_residual'] < 1e-6 and v['parallel_formula_error'] < 1e-6 and v['perpendicular_formula_error'] < 1e-6 for v in g6.values()),
        G6_speed_only_alpha_1=g6['1.0']['speed_change_fraction'] < 1e-8 and all(g6[k]['speed_change_fraction'] > 1e-3 for k in ('0.0', '0.5', '2.0')),
        G7_vector_isotropic_zero=g7['isotropic_mean_abs'] < 3*g7['standard_error']*np.sqrt(2),
        G7_vector_rotating=g7['rotating_error'] < 1e-3,
        G7_support_factor=g7['support_max_difference'] < 2e-3 and g7['isotropic_exact_error'] < 1e-10,
        G8_divergence=g8['divergence_relative'] < 1e-6,
        G8_transverse_focuses=g8['focusing_fraction']['transverse'] > .999,
        G8_vector_does_not=g8['focusing_fraction']['vector_max_abs'] < .1,
        G9_circular=all(g9['orbits'][k]['radius_drift'] < 1e-8 and abs(g9['orbits'][k]['period_over_analytic'] - 1) < 1e-8 for k in ('steering', 'well')),
        G9_radial_unturned=g9['orbits']['radial_steering']['direction_change'] < 1e-12 and g9['orbits']['radial_steering']['speed_change'] < 1e-12,
        G9_conservation=all(v['conservation_drift_per_orbit'] < 1e-8 for v in g9['ensemble']['laws'].values()))
    return dict(stage='exact', gates=gates, numerical_pass=all(gates.values()), field=dict(sha256=field['sha256'], max_abs_phi=field['max_abs_phi']),
                G1=g1, G2=probes, G5=g5, G6=g6, G7=g7, G8=g8, G9=g9, archived_ray_reproduction=repro,
                note='Fixtures: CWC-1 archived 2D field by hash; a solenoidal model field and a spherical steering field, neither an astrophysical source.')

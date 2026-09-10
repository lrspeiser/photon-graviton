"""Render the completed orbit diagnostics with their explicit acceptance gates."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'results.json').read_text())
protocol=json.loads((HERE/'protocol.json').read_text())
g=protocol['acceptance_checks']
lines=['## Orbital results','',
       'Thirty synthetic trajectories (six initial states under five potentials) were integrated over approximately 244 million years. The figures below describe these sampled paths; they are not observed stellar velocities, equilibrium population predictions or statistically independent astronomical tests.','',
       '| Potential | Largest Jacobi drift / (220 km/s)^2 | Tolerance change: position [kpc] | Tolerance change: velocity [km/s] | Retained relative tolerance | Numerical gates |',
       '|---|---:|---:|---:|---:|---|']
passed=True
for name in protocol['models']:
    r=d[name];jac=max(p['max_Jacobi_drift_over_220_squared'] for p in r['probes'])
    ok=(jac<g['max_Jacobi_drift_divided_by_220_squared'] and
        r['tolerance_position_difference_kpc']<g['max_tolerance_refinement_position_difference_kpc'] and
        r['tolerance_velocity_difference_kms']<g['max_tolerance_refinement_velocity_difference_kms'])
    passed &= ok
    lines.append(f"| {name.replace('_',' ')} | {jac:.3g} | {r['tolerance_position_difference_kpc']:.3g} | {r['tolerance_velocity_difference_kms']:.3g} | {r['retained_rtol']:.1g} | {'Pass' if ok else 'FAIL'} |")
frame=d['verification']['inertial_rotating_max_position_difference_kpc']
passed &= frame<g['max_inertial_rotating_coordinate_difference_kpc']
lines+=['',f'A separate integration in inertial coordinates, transformed back to the bar frame over approximately 49 million years, differs by at most {frame:.3g} kpc. '+('The declared numerical orbit gates pass.' if passed else '**At least one declared numerical orbit gate fails; this operator is not ready for stellar inference.**'),'',
        'The original tolerance comparisons, any tighter repeats and the initial interrupted attempt are retained. These finite tests measure numerical sensitivity, not rigorous error bounds for every possible orbit.','',
        '**Illustration of the proposed height effect: maximum sampled absolute height, in kpc.** These six initial states do not carry weights that represent the Galaxy.','',
        '| Initial R [kpc] | Initial z [kpc] | Ordinary matter | Equatorial deposits | Upper/lower deposits | Halo comparison |',
        '|---|---:|---:|---:|---:|---:|']
for i,state in enumerate(d['initial_conditions']):
    R=(state[0]**2+state[1]**2)**.5
    values=[d[m]['probes'][i]['max_abs_z_kpc'] for m in ['ordinary_matter','companion_equatorial','companion_caps','halo_comparison']]
    lines.append(f'| {R:.0f} | {state[2]:.1f} | '+' | '.join(f'{v:.4f}' for v in values)+' |')
bc=d['bar_order_comparison']
height_difference=max(abs(a['max_abs_z_kpc']-b['max_abs_z_kpc']) for a,b in zip(bc['L24_probes'],d['ordinary_matter']['probes']))
occupancy_difference=max(abs(a['sampled_plane_fraction']-b['sampled_plane_fraction']) for a,b in zip(bc['L24_probes'],d['ordinary_matter']['probes']))
lines+=['',f"Changing the bar from order 24 to 40 changes a sampled trajectory position by up to {bc['max_trajectory_position_difference_kpc']:.4g} kpc, the largest sampled height by up to {height_difference:.4g} kpc, and the sampled plane-occupancy fraction by up to {100*occupancy_difference:.3g} percentage points. This is field-truncation sensitivity, not a new physical effect or a measurement-error interval.",'',
        'A different trajectory does not automatically imply that a model better predicts the real bulge. Orbital phases, initial populations, duration, field resolution and selection all matter. Compare distributions under the likelihood contract before interpreting an apparent height or rotation trend as deposited gravity.','',
        '![Illustrative radial–vertical trajectories](orbit-comparison.png)','']
report=HERE/'report.md';text=report.read_text(encoding='utf-8')
before=text.split('## Orbital results')[0]
after=text.split('## What this enables and what remains')[1]
report.write_text(before+'\n'.join(lines)+'\n## What this enables and what remains'+after,encoding='utf-8',newline='\n')
(HERE/'gate-status.json').write_text(json.dumps(dict(orbit_numerical_gates_passed=bool(passed),
    stellar_likelihood_fitted=False,heldout_scores_evaluated=False,field_continuum_convergence_proven=False),indent=2)+'\n',encoding='utf-8',newline='\n')
print('Orbit numerical gates passed:',passed)

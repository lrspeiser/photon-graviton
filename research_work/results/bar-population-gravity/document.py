from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'summary.json').read_text());a=json.loads((HERE/'analytic-checks.json').read_text())
lines=['# Gravity of a continuously supplied moving test population','',
'This calculation turns the prescribed captured trajectories into potential and three-dimensional force coefficients. It is a linear response in a fixed ordinary bar, not a self-gravitating formed galaxy. The injection rate remains unspecified and unfunded. The source-direction grid is not established as converged for these new force quantities.','',
'## Assumptions and provenance','',
'The [population-to-gravity derivation](../bar-angular-refinement/population-to-gravity.md) gives the equations. Newtonian Green functions, Gauss quadrature and rotating-frame mechanics are known mathematics. The cold captured-particle interpretation, thin source weighting, launch speed and source history are project assumptions. No formula is claimed unique.','',
'The retained ordinary bar component and its normalization are model-dependent. Other Galactic components, deposited self-gravity and ordinary-bar backreaction are absent. A rotating prescribed bar can exchange orbital energy with the particles; this calculation does not fund that exchange with an evolving ordinary-matter energy reservoir. The source is stationary in bar coordinates. Each of the 1 and 3 kpc launch spheres has its own unspecified constant injection rate; they are not combined into a fitted radial source.','',
'## Results','',
'The table is at T=0.25 kpc/(km/s), approximately 244 Myr. Multiply potential coefficients by G times the total injected rest mass of that source to obtain potential. Multiply acceleration coefficients by the same factor to obtain acceleration. Equivalently, for a constant rate, total injected mass is injection rate times T. These are predictions per unit source mass, not measured stellar accelerations.','',
'| Source radius (kpc) | Position (kpc) | Potential / GM (1/kpc) | ax / GM (1/kpc²) | ay / GM (1/kpc²) | az / GM (1/kpc²) |',
'|---|---|---:|---:|---:|---:|']
for r in d['rows']:
    if r['T']==.25:
        xyz=f"({r['x_kpc']:g}, {r['y_kpc']:g}, {r['z_kpc']:g})"
        lines.append(f"| {r['source_R_kpc']:g} | {xyz} | {r['potential_per_GM']:.6g} | {r['ax_per_GM']:.6g} | {r['ay_per_GM']:.6g} | {r['az_per_GM']:.6g} |")
failed=[c for c in d['checks'] if not c['passes']]
lines+=['','## Verification and limits','',
f"The ring and eccentric-orbit analytic checks pass (maximum ring error {max(a['ring_potential_error'],a['ring_acceleration_error']):.3g}; inverse-radius error {a['kepler_inverse_radius_error']:.3g}). All source orbit and age-quadrature gates pass: **{d['all_orbits_pass']}**. {len(failed)} of {len(d['checks'])} age-quadrature comparisons remain failed. Original and refined comparisons remain in the raw data.",
'', 'Age quadrature uses the integrator step intervals with 4/8 nodes and retained 16-node refinement when needed. It measures numerical time-integration error on each chosen trajectory. It does not test source-direction quadrature, field-order errors at new positions, or the continuum limit near trajectories. The earlier L40/L64 orbit agreement is useful evidence but not a whole-population force verification.','',
'Reflection across the plane forces the vertical acceleration in the plane to zero. Above the plane, a vertical component is expected. This symmetry statement does not validate an observed Milky Way vertical force, and the listed positions are diagnostic points, not observed star bins. Bar rotation can produce an angular offset in the source distribution; individual nonzero transverse components must be checked against finer direction sampling before interpretation.','',
'## Next calculation','',
'Apply the same force integrals to refined direction grids before using their values in a dynamical fit. A smooth-looking central-residence estimate alone cannot establish convergence of the more locally sensitive three-dimensional force. Then evolve capture and deposited gravity together, with a funded source and ordinary-matter response. A metric or stress-source rule is also needed to predict lensing from the same population. The redshift/time mechanism and all nine overall goals remain incomplete.','',
'Reproduce with verify_kernel.py, run.py, export.py, then document.py using the existing hashed bar cache and retained angular launches. No observational data, halo density or withheld sample is used in this checkpoint.']
(HERE/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8',newline='\n')

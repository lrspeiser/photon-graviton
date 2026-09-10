"""Summarize duration evidence and its limits."""
from pathlib import Path
import json
import numpy as np

HERE = Path(__file__).resolve().parent
d = json.loads((HERE/'results.json').read_text())
a = json.loads((HERE/'assessment.json').read_text())
field = {r['seed_index']: r for r in a['rows']}
baseline = np.array([r['spatial_TV'] for r in a['baseline_72']])
lines = [
    '# Does the orbit library sample a stable population?', '',
    'Six declared probes extend the first low-height and high-height seed in each of three radial strata from 0.25 to 1.0 kpc/(km/s), approximately 244 to 978 million years. This is a duration check, not a new observed-star comparison or a fit of the gravity formula.', '',
    'The probe indices are 0, 16, 24, 40, 48 and 64. They were selected by their existing position in the library, not by the outcome of this extension. Seeds 16 and 40 retain their existing distance-disagreement flags. The same hypothetical empirical-response potential, bar rotation, initial states and numerical thresholds are used throughout.', '',
    '## What is measured?', '',
    'Divide space into fixed bins of galactic radius, signed height and angle relative to the bar. For each orbit, count the fraction of saved time in each bin. Compare the first half of an integration with its second half using the standard total-variation distance: TV = one half of the sum of the absolute differences in bin fractions. A value of 0 means identical binned time fractions; 1 means no overlap. This known statistical distance is not a new physics formula or a significance test.', '',
    'A large change means that these finite segments give different descriptions of where the orbit spends its time. A small change is useful but cannot prove equilibrium: it depends on the bins, duration and dynamics. Adjacent saved samples are correlated, so they are not thousands of independent stars. We also compare radius/height bins with angle summed out, and repeat the spatial comparison after discarding every second saved point.', '',
    '## Numerical checks', '',
    f"{a['numerical_passes']}/6 probes pass the unchanged trajectory-refinement and Jacobi-energy thresholds. {a['field_passes']}/6 pass the 1% sampled-path extra-force refinement check. Every attempted tolerance is recorded, including failures. The occupation values below are descriptive; rows failing numerical checks cannot support a precise dynamical conclusion.", '',
    '| Seed | Distance flag | Numerical pass | Maximum extra-force change |',
    '|---:|---|---|---:|',
]
for r in d['rows']:
    f = field[r['seed_index']]
    change = f"{100*f['max_extra_force_refinement']:.4f}%" if f.get('assessed') else 'not assessed'
    lines.append(f"| {r['seed_index']} | {r['distance_disagreement_flag']} | {r['numerical_pass']} | {change} |")
lines += ['', '## Change between equal early and late halves', '',
          '| Seed | TV at 244 Myr | TV at 489 Myr | TV at 978 Myr | Radius/height only, 978 Myr | Half sampling, 978 Myr |',
          '|---:|---:|---:|---:|---:|---:|']
for r in d['rows']:
    if 'occupation_windows' not in r:
        continue
    w = r['occupation_windows']
    lines.append(f"| {r['seed_index']} | {w[0]['spatial_TV']:.3f} | {w[1]['spatial_TV']:.3f} | {w[2]['spatial_TV']:.3f} | {w[2]['R_z_TV']:.3f} | {w[2]['half_sampling_spatial_TV']:.3f} |")
lines += ['',
    f"Across all original 72 short paths, the median early/late spatial TV is {np.median(baseline):.3f}, with a range {baseline.min():.3f}–{baseline.max():.3f}. This describes the original library at the declared binning, not an acceptance threshold. The six extensions cannot certify the remaining 66 orbits.", '',
    'For the numerically passing seed 24, spatial TV falls from 0.860 to 0.032; longer integration substantially improves this example. Passing seed 48 remains at 0.849. Its unwrapped angle relative to the rotating bar spans only about 0.495 of a turn and ends about 0.056 turns from its starting angle, while its cylindrical radius ranges from 3.77 to 8.01 kpc. Thus counting elapsed time alone does not establish coverage of the bar-relative motion. This is not a classification of its resonances or a claim about the true orbit of the observed star.', '',
    '## Consequences for fitting real stars', '',
    'Orbit integration accuracy and adequate population sampling are different requirements. Longer trajectories can improve time coverage without adding all of the positions and velocities absent from the 72-orbit library. The preceding training-coverage audit still applies: a single successful extended orbit does not repair coverage across the bulge and disk.', '',
    'Before fitting population weights, increase the variety of starting positions and velocities, test duration and library-size sensitivity, and apply equivalent population freedom to the ordinary-matter baseline. Use the same survey-selection and distance/motion-uncertainty treatment. Any regularization or smoothing must be chosen using training data, not tuned after seeing reserved outcomes.', '',
    'The two original long outward paths (seeds 66 and 71) are not among these declared probes. Their outer-domain continuation remains unresolved. This extension cannot be cited as resolving them or as proving the current library stationary. Nor does it resolve the two flagged distance estimates included here.', '',
    'The broader photon-companion hypothesis still needs a common physical account of redshift, event timing, companion retention, capture and the gravity measured by both stellar motions and lensing. This test uses known rotating-frame Hamilton equations and a known QUMOND-style empirical field; it supplies no new photon-to-gravity derivation. The deferred total energy-supply calculation remains unpassed.', '',
    '## Reproduction and integrity', '',
    'Run `run.py`, `analyze.py`, then `report.py`. The first script records input and trajectory hashes, all tolerances and domain failures. The analysis checks those hashes, compares field resolution at all saved extension positions and verifies bin counting independently with digitized cell indices. Large trajectories remain in the ignored data cache. No source catalog, force coefficient or holdout outcome is changed.', '',
]
(HERE/'report.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
print('Duration report written.')

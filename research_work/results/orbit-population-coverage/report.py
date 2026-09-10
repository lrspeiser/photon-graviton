"""Render the measured coverage diagnostic without treating it as a model score."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
d = json.loads((HERE / 'results.json').read_text())
def percent(a, b):
    return f'{100*a/b:.2f}%'

lines = [
    '# Orbit-library coverage before population fitting',
    '',
    'The 72-orbit library is numerically stable but does not yet provide adequate demonstrated coverage for a stellar-population comparison. This is a training-only diagnostic of the available saved paths, not evidence against or for the gravitational formula.',
    '',
    'We compare 27,812 selected training stars after excluding all 72 launch stars. Of these, 27,642 pass the existing positional screen. The flagged and unassessed stars remain in a separate all-training result; none is deleted from the underlying catalog. No validation or final-test observations are opened.',
    '',
    'For each star, we find saved orbit states within a specified three-dimensional spatial radius and ask whether any also lies within a specified three-dimensional velocity radius. We use the same Cartesian bar frame for both. These radii are diagnostic resolution choices, not measurement error bars. The inputs are conditional posterior means; distance and motion uncertainties have not been integrated into a likelihood.',
    '',
    'The search uses 36,000 states: 500 positive-time samples from each of 72 orbits. Initial-time samples are excluded as well as the launch stars themselves. The remaining training stars influenced the original sample scales and seed-selection pool, so this is not an independent holdout.',
    '',
    '## Coverage at several diagnostic resolutions',
    '',
    '| Spatial radius (kpc) | Velocity radius (km/s) | Position covered | Position and velocity covered |',
    '|---:|---:|---:|---:|',
]
for row in d['position_screened_training']:
    lines.append(f"| {row['radius_kpc']} | {row['velocity_radius_kms']:.0f} | {percent(row['spatially_covered'], row['stars'])} | {percent(row['position_and_velocity_covered'], row['stars'])} |")
lines += [
    '',
    'At 0.5 kpc and 50 km/s, 18,699 of 27,642 stars have nearby positions on the saved paths, but only 4,664 have a nearby position and velocity together. A fitted weighting cannot add a missing orbit state. Smooth kernels could spread support beyond these radii, but artificially broad smoothing would blur the physical prediction rather than demonstrate that the library is sufficient.',
    '',
    '## Which regions lack coverage?',
    '',
    'The following table uses 0.5 kpc and 50 km/s throughout. R is distance from the rotation axis; absolute z is distance above or below the disk plane. Counts are observed sample counts, not selection-corrected stellar densities.',
    '',
    '| R (kpc) | Absolute z (kpc) | Stars | Position covered | Position and velocity covered |',
    '|---|---|---:|---:|---:|',
]
for row in d['cells']:
    lines.append(f"| {row['R_range'][0]}–{row['R_range'][1]} | {row['absolute_z_range'][0]}–{row['absolute_z_range'][1]} | {row['stars']:,} | {percent(row['spatially_covered'], row['stars'])} | {percent(row['position_and_velocity_covered'], row['stars'])} |")
coarse = d['half_time_sampling'][4]
fine = d['position_screened_training'][4]
lines += [
    '',
    'The upper bulge and disk regions have particularly sparse velocity coverage. This cannot yet establish that the real bulge needs a different gravity correction: incomplete orbit coverage, uncertain stellar distances, the ordinary-matter model and the hypothesized extra potential can all affect this comparison.',
    '',
    '## Saved-time resolution and verification',
    '',
    f"Keeping every second saved positive-time sample reduces joint coverage at the middle diagnostic resolution from {fine['position_and_velocity_covered']:,} to {coarse['position_and_velocity_covered']:,} stars ({percent(fine['position_and_velocity_covered'], fine['stars'])} to {percent(coarse['position_and_velocity_covered'], coarse['stars'])}). The modest change suggests that simply saving twice as many nearby points is unlikely to remove the large gap. It is not a continuum bound: longer paths, additional starts and denser time sampling still require separate checks.",
    '',
    'Twelve fixed training indices were independently checked against exhaustive distances to all 36,000 states. Search results agree exactly. Larger spatial windows never reduce coverage; thinning the saved states never increases it. Source hashes are checked before and after the calculation. These checks verify the search, not a stellar-population model.',
    '',
    'The nearest-neighbor distances use standard Euclidean geometry. The orbit equations are the known rotating-frame Hamilton equations, and the prescribed extra potential uses known QUMOND-style mathematics with empirical coefficients. No new photon-conversion or capture formula is derived by this diagnostic.',
    '',
    '## Consequence for the research goal',
    '',
    'Do not tune the gravity coefficients to repair this coverage table or present these fractions as predicted stellar agreement. First expand phase-space coverage using a declared training-only rule, test duration/occupation stability, and provide equivalent orbit-population freedom under ordinary matter. Then fit nonnegative population weights with survey selection and correlated distance/motion uncertainties. Library size and regularization choices must be fixed before new validation outcomes are inspected.',
    '',
    'For the next library, record separate coverage of signed height, bar angle and velocity components, not just nine broad radius/height bins. Retain the original results and distance flags. The two long outward paths in the current library also require an explicit outer-field/domain and duration treatment before stationary population weights can be interpreted.',
    '',
    'The larger scientific requirements remain unchanged: a shared prediction of redshift and event timing, physically specified energy transfer and retention, capture into a gravitational source, and a common response for stellar motions and lensing. The total source-energy budget remains deferred rather than passed. Numerical orbit success alone supplies none of these missing causal links.',
    '',
    '## Reproduction',
    '',
    'Run `run.py`, then `report.py` from this folder with the existing catalog and trajectory caches. `results.json` records aggregate counts, input hashes and the diagnostic choices. Individual training-only distances are saved in the ignored cache `research_work/data-cache/orbit-population-coverage/training-coverage.parquet`. Source IDs retain integer precision. No source catalog, force coefficient or holdout role is changed.',
    '',
]
(HERE / 'report.md').write_text('\n'.join(lines), encoding='utf-8', newline='\n')
print('Coverage report written.')

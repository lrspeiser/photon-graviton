"""Audit cached parent measurements; assign no distances or gravity predictions."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.table import Table

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache/cepheid-stars'
t = Table.read(CACHE / 'gaia-dr3-dcep-parent.votable', format='votable')
t.rename_columns(t.colnames, [c.lower() for c in t.colnames])
d = t.to_pandas()
for c in d.select_dtypes(include='object'):
    d[c] = d[c].map(lambda x: x.decode() if isinstance(x, bytes) else x)

def finite(columns):
    return np.isfinite(d[columns].astype(float)).all(axis=1)

flags = {
    'has_sky_pm': finite(['ra', 'dec', 'pmra', 'pmdec']),
    'has_intensity_photometry': finite(['int_average_g', 'int_average_bp', 'int_average_rp']),
    'has_positive_period': d[['pf', 'p1_o', 'p2_o']].gt(0).any(axis=1),
    'has_combined_rv': finite(['radial_velocity', 'radial_velocity_error']) & d.radial_velocity_error.gt(0),
    'has_pulsation_mean_rv': finite(['average_rv', 'average_rv_error']) & d.average_rv_error.gt(0),
    'ruwe_below_1_4': d.ruwe.lt(1.4),
    'combined_rv_transits_at_least_8': d.rv_nb_transits.ge(8),
    'has_astrometric_covariance_fields': finite([c for c in d if c.endswith('_corr')] + ['ra_error','dec_error','parallax_error','pmra_error','pmdec_error']),
}
for k, v in flags.items():
    d['flag_' + k] = v
cutflow = {'parent': len(d)}
keep = np.ones(len(d), dtype=bool)
for k in ['has_sky_pm', 'has_intensity_photometry', 'has_positive_period', 'has_combined_rv', 'ruwe_below_1_4', 'combined_rv_transits_at_least_8', 'has_astrometric_covariance_fields']:
    keep &= flags[k]
    cutflow[k] = int(keep.sum())
d['flag_measurement_candidate'] = keep
giants = pd.read_parquet(ROOT / 'research_work/data-cache/stellar-catalogs/matched-with-gaia-covariance.parquet', columns=['source_id', 'gaia_quality_candidate'])
d['overlaps_prepared_giant_parent'] = d.source_id.isin(giants.source_id)
d['overlaps_quality_giants'] = d.source_id.isin(giants.loc[giants.gaia_quality_candidate, 'source_id'])
assert d.source_id.dtype.kind in 'iu' and d.source_id.is_unique
out = CACHE / 'gaia-dr3-dcep-parent-with-flags.parquet'
d.to_parquet(out, index=False)
overlap = d.loc[d.overlaps_prepared_giant_parent, ['source_id', 'mode_best_classification', 'flag_measurement_candidate', 'overlaps_quality_giants']]
overlap.to_json(HERE / 'overlap.json', orient='records', indent=2)
result = {
    'rows': len(d),
    'mode_counts': {str(k): int(v) for k,v in d.mode_best_classification.value_counts(dropna=False).items()},
    'independent_flag_counts': {k: int(v.sum()) for k,v in flags.items()},
    'sequential_measurement_cutflow': cutflow,
    'overlap_with_prepared_giant_parent': int(d.overlaps_prepared_giant_parent.sum()),
    'overlap_with_quality_giants': int(d.overlaps_quality_giants.sum()),
    'measurement_candidates_overlapping_quality_giants': int((keep & d.overlaps_quality_giants).sum()),
    'giant_parent_rows': len(giants),
    'giant_quality_rows': int(giants.gaia_quality_candidate.sum()),
    'derived_parquet_sha256': hashlib.sha256(out.read_bytes()).hexdigest(),
    'limitations': ['No distance calibration or spatial cuts applied.', 'Not a reconstruction of the published 903-star sample.', 'Covariance field completeness is not a positive-semidefinite covariance check.', 'RV transit count and pulsation clean-epoch count are distinct.', 'Different identifiers do not establish independent survey systematics.', 'No gravity fit or new holdout score.'],
}
(HERE / 'audit.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
print(json.dumps(result, indent=2))

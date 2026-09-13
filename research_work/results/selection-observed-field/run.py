"""Audit observed targeting metadata only; no held-out stellar motions read."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache'
paths = {
    'allstar': CACHE / 'stellar-catalogs/allStarLite-dr17-synspec_rev1.fits',
    'parent': CACHE / 'selection-fields/apogee2Object_300+00.fits',
    'training': CACHE / 'stellar-selection-membership/training-membership.parquet',
}

def strings(a):
    return np.char.strip(np.asarray(a).astype(str))

def frame(table, columns, mask):
    out = {}
    for col in columns:
        a = np.array(table[col][mask])
        out[col] = strings(a) if a.dtype.kind in 'SU' else a.astype(a.dtype.newbyteorder('='))
    return pd.DataFrame(out)

cols = ['APOGEE_ID', 'EXTRATARG', 'J', 'H', 'K', 'AK_TARG',
        'AK_TARG_METHOD', 'MIN_H', 'MAX_H', 'MIN_JK', 'MAX_JK', 'TARGFLAGS']
with fits.open(paths['allstar'], memmap=True) as hd:
    t = hd[1].data
    mask = (strings(t['FIELD']) == '300+00') & (strings(t['TELESCOPE']) == 'lco25m')
    observed = frame(t, cols, mask)
main = observed.loc[observed.EXTRATARG.eq(0)].copy()
assert main.APOGEE_ID.is_unique, 'Resolve duplicate main targets before counting'
with fits.open(paths['parent'], memmap=True) as hd:
    t = hd[1].data
    parent = frame(t, ['APOGEE_ID', 'J', 'H', 'K', 'AK_TARG', 'AK_TARG_METHOD'], slice(None))
assert parent.APOGEE_ID.is_unique
joined = main.merge(parent, on='APOGEE_ID', how='left', suffixes=('_observed', '_parent'), indicator=True, validate='one_to_one')
training = pd.read_parquet(paths['training'], columns=['APOGEE_ID', 'FIELD', 'TELESCOPE', 'status'])
training = training.loc[training.FIELD.eq('300+00') & training.TELESCOPE.eq('lco25m') & training.status.eq('main_red')]
assert set(training.APOGEE_ID).issubset(set(main.APOGEE_ID))
comparison = {}
for c in ['J', 'H', 'K', 'AK_TARG']:
    a, b = joined[c + '_observed'].to_numpy(), joined[c + '_parent'].to_numpy()
    good = np.isfinite(a) & np.isfinite(b)
    comparison[c] = dict(both_finite=int(good.sum()), finite_status_disagreements=int((np.isfinite(a) != np.isfinite(b)).sum()),
                         differences_over_1e_5=int((np.abs(a[good]-b[good]) > 1e-5).sum()),
                         max_abs_difference=float(np.max(np.abs(a[good]-b[good]))) if good.any() else None)
bounds = main.groupby(['MIN_H', 'MAX_H', 'MIN_JK', 'MAX_JK'], dropna=False).size().reset_index(name='observed_main_count')
h_inside = main.H.ge(main.MIN_H) & main.H.le(main.MAX_H)
out = CACHE / 'selection-observed-field'
out.mkdir(exist_ok=True)
main.to_parquet(out / 'observed-main-metadata.parquet', index=False)
def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()
result = dict(scope='Full observed main-sample metadata for one training-selected field; not selection probabilities',
    inputs={str(p.relative_to(ROOT)): sha(p) for p in paths.values()},
    field='300+00', telescope='lco25m', all_observed_rows=len(observed),
    extratarg_counts={str(k): int(v) for k,v in observed.EXTRATARG.value_counts().items()},
    main_rows=len(main), unique_main_ids=main.APOGEE_ID.nunique(), training_main_rows=len(training),
    parent_matches=int(joined._merge.eq('both').sum()), missing_parent_ids=joined.loc[joined._merge.ne('both'), 'APOGEE_ID'].tolist(),
    metadata_comparison=comparison, extinction_method_agreement=int(joined.AK_TARG_METHOD_observed.eq(joined.AK_TARG_METHOD_parent).sum()),
    observed_bounds=bounds.to_dict('records'), h_outside_recorded_inclusive_bounds=main.loc[~h_inside, 'APOGEE_ID'].tolist(),
    documentation='https://www.sdss4.org/dr17/irspec/targets/selection-biases/',
    heldout_kinematic_columns_read=False, completeness_estimated=False)
(HERE / 'results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps({k:v for k,v in result.items() if k != 'inputs'}, indent=2))

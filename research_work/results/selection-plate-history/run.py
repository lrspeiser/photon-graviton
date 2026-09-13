"""Read actual drilled-target metadata for the pilot selection field."""
from pathlib import Path
import hashlib
import json
import re
import requests
import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache/selection-fields'
BASE = 'https://svn.sdss.org/public/data/sdss/platelist/trunk/'
plate_path = CACHE / 'apogee2Plate.fits'
plates = fits.getdata(plate_path, 1)
plates = plates[plates['LOCATION_ID'] == 5540]
manifest, summaries, records = [], [], []
for plate in plates:
    pid = int(plate['PLATE_ID'])
    filename = f'plateHolesSorted-{pid:06d}.par'
    url = BASE + f'plates/{pid//100:04d}XX/{pid:06d}/' + filename
    path = CACHE / filename
    if not path.exists():
        response = requests.get(url, timeout=(15, 45))
        response.raise_for_status()
        assert response.text.startswith('plateId ')
        path.write_bytes(response.content)
    content = path.read_text()
    assert int(re.search(r'^plateId (\d+)', content, re.M)[1]) == pid
    assert int(re.search(r'^designid (\d+)', content, re.M)[1]) == int(plate['DESIGN_ID'])
    body = re.search(r'typedef struct \{(.*?)\} STRUCT1;', content, re.S)[1]
    columns = re.findall(r'\b(?:char|double|int|float|long)\s+(\w+)(?:\[\d+\])?;', body)
    rows = []
    for line in content.splitlines():
        if not line.startswith('STRUCT1 '):
            continue
        tokens = re.findall(r'"[^"]*"|\{[^}]*\}|[^\s]+', line)[1:]
        assert len(tokens) == len(columns), (pid, len(tokens), len(columns))
        row = dict(zip(columns, [v.strip('"') for v in tokens]))
        if row['targettype'].lower() == 'science':
            keep = ['target_ra', 'target_dec', 'targetids', 'tmass_j', 'tmass_h', 'tmass_k',
                    'fiberid', 'assigned', 'conflicted', 'priority', 'apogee2_target1', 'apogee2_target2', 'apogee2_target3']
            rows.append(dict(plate_id=pid, design_id=int(plate['DESIGN_ID']), **{k:row[k] for k in keep}))
    summaries.append(dict(plate_id=pid, design_id=int(plate['DESIGN_ID']), science_rows=len(rows),
                          assigned_science_rows=sum(int(r['assigned']) == 1 for r in rows)))
    records.extend(rows)
    manifest.append(dict(url=url, path=str(path.relative_to(ROOT)), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    print(summaries[-1], flush=True)
table = pd.DataFrame(records)
assert table.targetids.str.strip().str.match(r'^2MASS-J\d{8}[+-]\d{7}$').all()
table['APOGEE_ID'] = table.targetids.str.strip().str.replace('2MASS-J', '2M', regex=False)
bits = table.apogee2_target1.astype('int64')
# Documented APOGEE2_TARGET1 bits; signed int32 values retain these low bits.
table['normal_sample'] = bits.map(lambda v: bool(v & (1 << 14)))
table['cohort'] = bits.map(lambda v: ','.join(name for bit,name in [(11,'short'),(12,'medium'),(13,'long')] if v & (1 << bit)))
main_path = ROOT / 'research_work/data-cache/selection-observed-field/observed-main-metadata.parquet'
main = pd.read_parquet(main_path)
missing = sorted(set(main.APOGEE_ID) - set(table.APOGEE_ID))
assert not missing
long_ids = main.loc[main.MIN_H > 12.7, 'APOGEE_ID'].tolist()
long_rows = table.loc[table.APOGEE_ID.isin(long_ids)]
assert len(long_rows) == 16 and long_rows.normal_sample.all() and long_rows.cohort.eq('long').all()
assert long_rows.groupby('APOGEE_ID').plate_id.nunique().eq(8).all()
for summary in summaries:
    part = table.loc[table.plate_id.eq(summary['plate_id']) & table.normal_sample]
    summary['normal_sample_cohorts'] = part.cohort.value_counts().to_dict()
unique_normal = table.loc[table.normal_sample].drop_duplicates('APOGEE_ID')
assert set(main.APOGEE_ID).issubset(set(unique_normal.APOGEE_ID))
normal_unobserved = sorted(set(unique_normal.APOGEE_ID) - set(main.APOGEE_ID))
allstar_path = ROOT / 'research_work/data-cache/stellar-catalogs/allStarLite-dr17-synspec_rev1.fits'
with fits.open(allstar_path, memmap=True) as hd:
    t = hd[1].data
    selected = (np.char.strip(np.asarray(t['FIELD']).astype(str)) == '300+00') & (np.char.strip(np.asarray(t['TELESCOPE']).astype(str)) == 'lco25m')
    ids = np.char.strip(np.asarray(t['APOGEE_ID'][selected]).astype(str))
    flags = t['EXTRATARG'][selected]
    assert len(set(ids)) == len(ids)
    observed_flags = dict(zip(ids, [int(v) for v in flags]))
not_main_status = pd.Series([str(observed_flags[sid]) if sid in observed_flags else 'absent_from_field_allstar' for sid in normal_unobserved]).value_counts().to_dict()
design_path = CACHE / 'apogee2Design.fits'
design_url = 'https://data.sdss.org/sas/dr17/apogee/target/apogee2Design.fits'
if not design_path.exists():
    response = requests.get(design_url, timeout=(15,45))
    response.raise_for_status()
    design_path.write_bytes(response.content)
designs = fits.getdata(design_path, 1)
designs = designs[designs['LOCATION_ID'] == 5540]
assert len(designs) == 8 and np.all(designs['COHORT_FRACTION'][:,2] == 0)
manifest.append(dict(url=design_url, path=str(design_path.relative_to(ROOT)), sha256=hashlib.sha256(design_path.read_bytes()).hexdigest()))
out = ROOT / 'research_work/data-cache/selection-plate-history'
out.mkdir(exist_ok=True)
table.to_parquet(out / 'drilled-science-metadata.parquet', index=False)
result = dict(scope='Drilled plate metadata, not visit completion or random inclusion probability',
              manifest=manifest, plates=summaries, target_id_examples=table.targetids.head(4).tolist(),
              plate_summary_sha256=hashlib.sha256(plate_path.read_bytes()).hexdigest(),
              observed_metadata_sha256=hashlib.sha256(main_path.read_bytes()).hexdigest(),
              unique_science_targets=int(table.APOGEE_ID.nunique()),
              unique_normal_targets=len(unique_normal), observed_main_matches=len(main),
              observed_main_missing_from_plates=missing,
              normal_drilled_not_in_observed_main_count=len(normal_unobserved),
              normal_drilled_not_in_observed_main_ids=normal_unobserved,
              normal_drilled_not_in_observed_main_status=not_main_status,
              zero_planned_long_fraction_designs=designs['DESIGN_ID'].astype(int).tolist(),
              long_cohort_assignments=long_rows.to_dict('records'),
              bit_documentation='https://www.sdss4.org/dr17/algorithms/bitmasks/#APOGEE2_TARGET1',
              heldout_kinematics_read=False)
(HERE / 'results.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
print(result['target_id_examples'])

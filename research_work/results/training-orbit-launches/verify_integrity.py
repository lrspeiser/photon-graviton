"""Check archived provenance and numerical gates without opening holdouts."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import prepare

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def read(name):
    return json.loads((HERE / name).read_text(encoding='utf-8'))

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

checked = {}
def check(path, expected):
    actual = digest(path)
    assert actual == expected, str(path)
    checked[str(path.relative_to(ROOT))] = actual

prep = read('results.json')
integ = read('integration-results.json')
field = read('field-check.json')
epoch = read('epoch-check.json')
for name, expected in prep['source_hashes'].items():
    check(prepare.paths[name], expected)
for mapping in [prep['code_hashes'], integ['input_field_hashes'],
                field['input_hashes'], epoch['source_hashes']]:
    for path, expected in mapping.items():
        check(ROOT / path, expected)
check(HERE / 'integrate.py', integ['code_sha256'])
for query in epoch['queries']:
    check(HERE / query['file'], query['sha256'])
seeds = read('launches.json')
ids = [row['source_id'] for row in seeds]
assert len(ids) == len(set(ids)) == 72
assert all(row['holdout_role'] == 'training' for row in seeds)
table = pd.read_parquet(prepare.OUT / 'launches.parquet')
assert table.source_id.astype(str).tolist() == ids
with np.load(prepare.OUT / 'full-field-orbits.npz') as paths:
    assert paths['source_id'].astype(str).tolist() == ids
for data in [integ, field, epoch]:
    assert [row['source_id'] for row in data['rows']] == ids
for row in integ['rows']:
    last = row['attempts'][-1]
    assert row['error'] is None
    assert last['position_difference_kpc'] < integ['position_gate_kpc']
    assert last['velocity_difference_kms'] < integ['velocity_gate_kms']
    assert row['Jacobi_drift_over_220_squared'] < integ['Jacobi_gate_over_220_squared']
assert all(row['positions'] == 501 and row['maximum_fraction'] < .01 for row in field['rows'])
assert all(row['separation_arcsec'] <= .5 for row in epoch['rows'])
assert sum(row['parallax_interval_disjoint_5'] for row in seeds) == 3
result = dict(passed=True, seeds=72, sampled_positions=36072,
              checked_hashes=checked, training_only=True,
              population_fit_validated=False, holdouts_opened=False)
(HERE / 'integrity.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
print('Integrity and archived numerical gates passed; no population validation.')

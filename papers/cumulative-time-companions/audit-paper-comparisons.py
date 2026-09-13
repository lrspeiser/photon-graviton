"""Extract existing comparison results and check the manuscript's rounded tables.

No observations are downloaded and no physical parameters are fitted.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sources = {}


def read(relative):
    path = ROOT / 'research_work/results' / relative
    sources[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return json.loads(path.read_text(encoding='utf-8'))


rotation = read('isotropic-galaxy-transfer/model-comparison-results.json')
cluster = read('isotropic-galaxy-transfer/cluster-comparison-detail-results.json')
redshift = read('direct-conversion/results.json')
brightness = read('brightness-distance-consistency/rate-results.json')
paper = (HERE / 'theory-basis.md').read_text(encoding='utf-8')
checked = []
for row in rotation['summary']:
    if row['split'] in ('validation', 'test') and row['model'] in ('baryons', 'companion_third', 'MOND_simple_fitted', 'MOND_simple_fixed', 'NFW_shared_scaling'):
        value = f"{row['RMSE_kms']:.2f}"
        assert value in paper, (row['model'], value)
        checked.append(value)
cluster_table = {}
for name, row in cluster['datasets']['published'].items():
    cluster_table[name] = {'all_chi2': row['full']['all_sum'], 'omitted_sum': row['loo_sum']}
    for value in cluster_table[name].values():
        assert f'{value:.2f}' in paper, (name, value)
        checked.append(f'{value:.2f}')
for name, value in redshift['scales'].items():
    assert f'{value:.3f}' in paper, (name, value)
for name, row in redshift['coarse_oof_statistics'].items():
    assert f"{row['rms_km_s']:.3f}" in paper, name
assert f"{brightness['c_alpha_kms_mpc']:.3f}" in paper
assert '474646' not in paper
output = {
    'version': '1.4', 'evidence_through': '0756dec',
    'scope': 'Snapshot of existing exposed-data diagnostics; no new fitting or full-theory ranking.',
    'source_sha256': sources,
    'rotation': rotation['summary'],
    'rotation_outer': rotation['outer_summary'],
    'cluster': cluster_table,
    'redshift_scales': redshift['scales'],
    'redshift_omitted_region': redshift['coarse_oof_statistics'],
    'brightness_rate': brightness['c_alpha_kms_mpc'],
    'checks': 'Rounded main rotation, Coma, redshift scale and omitted-region values present in manuscript.'
}
(HERE / 'paper-comparison-audit.json').write_text(json.dumps(output, indent=2) + '\n', encoding='utf-8', newline='\n')
print('Comparison extraction and manuscript value checks passed.')

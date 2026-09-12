"""Export existing completed runs; no fit or reserved data access."""
from pathlib import Path
import csv
import json

HERE = Path(__file__).resolve().parent
rows = []
summary = []
for radius in (30, 100, 200):
    source = (HERE.parent / 'coupled-capture-calibration/results.json'
              if radius == 100 else HERE / f'outer{radius}/results.json')
    data = json.loads(source.read_text())
    assert all(c['passes'] for c in data['checks'])
    for fit in data['fits']:
        assert not fit['edge_optimum']
        result = fit['refined']
        summary.append(dict(outer_kpc=radius, p=fit['p'], C=result['C'],
                            rms=result['rms'], mass_Msun=result['mass_Msun'],
                            solar_depth=result['solar_deposit_depth']))
        for i, (prediction, residual) in enumerate(zip(result['predictions'], result['residuals'])):
            rows.append(dict(outer_kpc=radius, p=fit['p'], bin=i,
                             observed_training_proxy_kms=prediction-residual,
                             predicted_kms=prediction, residual_kms=residual,
                             role='exposed training; not validation'))
with (HERE / 'comparison.csv').open('w', encoding='utf8', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
print(json.dumps(summary, indent=2))

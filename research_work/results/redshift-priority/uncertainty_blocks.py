"""Sampling sensitivity of fixed exploratory predictions, not new validation."""
import csv
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
C = 299792.458


def metrics(x):
    return {'rms': float(np.sqrt(np.mean(x*x))), 'bias': float(np.mean(x))}


if __name__ == '__main__':
    source = HERE/'coarse-sky-predictions.csv'
    rows = list(csv.DictReader(source.open()))
    assert len(rows) == 164 and len({r['group_pgc'] for r in rows}) == 164
    region = np.array([int(r['region']) for r in rows])
    error = C*np.array([float(r['constant_oof_z'])-float(r['observed_z']) for r in rows])
    previous = json.loads((HERE/'coarse-sky-results.json').read_text())
    for k,v in metrics(error).items():
        assert abs(v-previous['out_of_fold']['constant'][k]) < 1e-9
    groups = [np.flatnonzero(region == r) for r in sorted(set(region))]
    assert len(groups) == 8
    summary = [{'region': int(r), 'n': len(index), **metrics(error[index]),
                'median_absolute': float(np.median(abs(error[index]))),
                'max_absolute': float(np.max(abs(error[index])))}
               for r,index in zip(sorted(set(region)),groups)]
    rng = np.random.default_rng(2026090917)
    samples = {mode: [] for mode in ('individual_groups','whole_sky_regions')}
    for _ in range(10000):
        samples['individual_groups'].append(metrics(error[rng.integers(0,len(error),len(error))]))
        idx = np.concatenate([groups[j] for j in rng.integers(0,len(groups),len(groups))])
        samples['whole_sky_regions'].append(metrics(error[idx]))
    intervals = {mode: {metric: np.percentile([s[metric] for s in values],[2.5,97.5]).tolist()
                         for metric in ('rms','bias')} for mode,values in samples.items()}
    ranks = []
    for n in (3,7,8,19,64):
        for coverage in (Fraction(4,5),Fraction(9,10),Fraction(19,20)):
            raw = (n+1)*coverage
            k = -(-raw.numerator//raw.denominator)
            ranks.append({'calibration_units': n, 'nominal_coverage': float(coverage),
                          'required_rank': k, 'threshold': 'finite order statistic' if k<=n else '+infinity'})
    out = {'status': 'exploratory exposed-sample sensitivity, not a predictive coverage guarantee',
           'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
           'protocol_sha256': hashlib.sha256((HERE/'uncertainty-block-protocol.md').read_bytes()).hexdigest(),
           'groups':len(rows), 'sky_regions':len(groups), 'fixed_prediction_metrics_kms':metrics(error),
           'region_metrics_kms':summary, 'conditional_bootstrap_percentile_95_kms':intervals,
           'replicates':10000, 'seed':2026090917,
           'split_conformal_rank_diagnostic_not_applied_to_data':ranks}
    (HERE/'uncertainty-block-results.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k!='split_conformal_rank_diagnostic_not_applied_to_data'},indent=2))

"""Independent foreground map lookup and exposed-residual diagnostic."""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT/'redshift_paper/all_164_groups.csv'
FIELDS = ('refPixelValueSandF','refPixelValueSFD','meanValueSandF','meanValueSFD','stdSandF','stdSFD')


def fetch(row):
    url = 'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?' + urllib.parse.urlencode(
        {'locstr':f"{row['ra_deg']} {row['dec_deg']} equ j2000",'regSize':2})
    out = {'group_pgc':row['group_pgc'],'ra_deg':float(row['ra_deg']),
           'dec_deg':float(row['dec_deg']),'url':url}
    try:
        raw = urllib.request.urlopen(url,timeout=40).read()
        root = ET.fromstring(raw)
        assert root.attrib.get('status') == 'ok'
        matches = [b for b in root.findall('result') if 'Reddening' in (b.findtext('desc') or '')]
        assert len(matches) == 1
        stats = matches[0].find('statistics')
        out['statistics_mag'] = {k:float(stats.findtext(k).strip().split()[0]) for k in FIELDS}
        out['refCoordinate'] = stats.findtext('refCoordinate').strip()
        out['response_sha256'] = hashlib.sha256(raw).hexdigest()
    except Exception as exc:
        out['error'] = f'{type(exc).__name__}: {exc}'
    return out


def corr(x,y):
    return {'pearson':float(np.corrcoef(x,y)[0,1]),'spearman':float(spearmanr(x,y).statistic)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser();parser.add_argument('--fetch',action='store_true');args = parser.parse_args()
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == '8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
    rows = list(csv.DictReader(SOURCE.open()))
    snapshot = HERE/'foreground-dust-metadata.json'
    if args.fetch:
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            records = list(pool.map(fetch,rows))
        snapshot.write_text(json.dumps({'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
            'protocol_sha256':hashlib.sha256((HERE/'foreground-dust-protocol.md').read_bytes()).hexdigest(),
            'records':records},indent=2)+'\n',newline='\n')
    records = json.loads(snapshot.read_text())['records']
    assert len(records) == 164 and {r['group_pgc'] for r in records} == {r['group_pgc'] for r in rows}
    original = {r['group_pgc']:r for r in rows}
    for r in records:
        assert r['ra_deg'] == float(original[r['group_pgc']]['ra_deg'])
        assert r['dec_deg'] == float(original[r['group_pgc']]['dec_deg'])
        if 'error' not in r:
            returned = list(map(float,r['refCoordinate'].split()[:2]))
            assert abs(returned[0]-r['ra_deg']) <= 1e-6 and abs(returned[1]-r['dec_deg']) <= 1e-6
    failures = [r['group_pgc'] for r in records if 'error' in r]
    assert not failures, f'Resolve requests before analysis: {failures}'
    records = {r['group_pgc']:r for r in records}
    prediction_source = HERE/'coarse-sky-predictions.csv'
    predictions = list(csv.DictReader(prediction_source.open()))
    assert {r['group_pgc'] for r in predictions} == set(records)
    ebv = np.array([records[r['group_pgc']]['statistics_mag']['refPixelValueSandF'] for r in predictions])
    residual = 299792.458*np.array([float(r['constant_oof_z'])-float(r['observed_z']) for r in predictions])
    region = np.array([int(r['region']) for r in predictions])
    assert np.isfinite(ebv).all() and np.isfinite(residual).all()
    def demean(x):
        return x-np.array([np.mean(x[region==r]) for r in region])
    results = {'n':len(ebv),'primary_input':'Schlafly-Finkbeiner reference-pixel E(B-V), mag',
        'input_range_mag':[float(ebv.min()),float(ebv.max())],
        'correlations':{'signed_residual':corr(ebv,residual),'absolute_residual':corr(ebv,abs(residual)),
                        'region_demeaned_signed':corr(demean(ebv),demean(residual)),
                        'region_demeaned_absolute':corr(demean(ebv),demean(abs(residual)))},
        'regions':[{'region':int(r),'n':int(sum(region==r)),
                    'mean_ebv_mag':float(ebv[region==r].mean()),'mean_residual_kms':float(residual[region==r].mean())}
                   for r in sorted(set(region))],
        'metadata_sha256':hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        'prediction_sha256':hashlib.sha256(prediction_source.read_bytes()).hexdigest(),
        'status':'descriptive exposed-sample audit; no distance correction or physical environmental term fitted'}
    (HERE/'foreground-dust-results.json').write_text(json.dumps(results,indent=2)+'\n',newline='\n')
    print(json.dumps(results,indent=2))

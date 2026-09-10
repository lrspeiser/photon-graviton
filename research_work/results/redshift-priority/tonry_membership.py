"""Source-name/position join; keep positional evidence distinct from provenance weights."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from elves_radio_quality import separation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
HASHES = {'good':'301884b80963ef0b69b547ddca86caff1c893830ccd587e9a0064e68616a84c7',
          'poor':'c3dda46762ee68c7424e70263bbe79ea1ef285bcdc5d12671ca7eceba55c2c28'}


def parse(fields):
    return {'source_name':fields[0], 'ra_deg':float(fields[1]), 'dec_deg':float(fields[2])}


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('cache',type=Path);args=ap.parse_args()
    sources=[];counts={}
    for category,expected in HASHES.items():
        raw=(args.cache/f'tonry-table.{category}').read_bytes()
        assert hashlib.sha256(raw).hexdigest()==expected
        lines=raw.decode().splitlines();assert lines[0].startswith('Pref_name')
        counts[category]=0
        for line in lines[1:]:
            if not line.strip():continue
            fields=line.split();assert len(fields)==22
            record=parse(fields)
            assert record==parse(fields[:3]+['MASKED']*(len(fields)-3))
            sources.append({**record,'source_class':category});counts[category]+=1
    assert len({r['source_name'] for r in sources})==len(sources)
    master=ROOT/'redshift_paper/all_164_groups.csv'
    # Explicit hash below avoids depending on row order from other diagnostic tables.
    assert hashlib.sha256(master.read_bytes()).hexdigest()=='8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
    rows=list(csv.DictReader(master.open()));assert len(rows)==164
    predictions={r['group_pgc']:r for r in csv.DictReader((HERE/'coarse-sky-predictions.csv').open())}
    results=[]
    for row in rows:
        point={'ra_deg':float(row['ra_deg']),'dec_deg':float(row['dec_deg'])}
        matches=[]
        for source in sources:
            angle=separation(point,source)
            if angle<=30:matches.append({**source,'separation_arcsec':angle})
        p=predictions[row['group_pgc']]
        results.append({'pgc':row['pgc'],'group_pgc':row['group_pgc'],'region':int(p['region']),
            'matches':matches,'status':'positional candidate only' if matches else 'no match within fixed radius',
            'existing_distance_mpc':float(row['catalog_distance_mpc']),
            'existing_residual_kms':299792.458*(float(p['constant_oof_z'])-float(p['observed_z']))})
    summary={}
    for matched in (True,False):
        subset=[r for r in results if bool(r['matches'])==matched]
        e=np.array([r['existing_residual_kms'] for r in subset]);d=[r['existing_distance_mpc'] for r in subset]
        summary['matched' if matched else 'unmatched']={'n':len(subset),'rms_kms':float(np.sqrt(np.mean(e*e))),
            'bias_kms':float(e.mean()),'distance_range_mpc':[min(d),max(d)],
            'region_counts':{str(i):sum(r['region']==i for r in subset) for i in range(8)}}
    out={'source_urls':{k:f'https://www.ifa.hawaii.edu/~jt/SBF/table.{k}' for k in HASHES},
         'source_hashes':HASHES,'source_counts':counts,'source_records':sources,'target_count':len(results),
         'multiple_match_targets':sum(len(r['matches'])>1 for r in results),'summary':summary,'records':results,
         'protocol_sha256':hashlib.sha256((HERE/'tonry-membership-protocol.md').read_bytes()).hexdigest(),
         'status':'positional source membership, not confirmed CF4 contribution weights; no correction fitted'}
    (HERE/'tonry-membership.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('records','source_records')},indent=2))

"""Fetch explicitly allowlisted SIMBAD source metadata, never outcome columns."""
import concurrent.futures
import hashlib
import json
import urllib.parse
import urllib.request
from pathlib import Path

HERE=Path(__file__).resolve().parent
ENDPOINT='https://simbad.cds.unistra.fr/simbad/sim-tap/sync'


def query(adql):
    url=ENDPOINT+'?'+urllib.parse.urlencode({'request':'doQuery','lang':'adql','format':'json','query':adql})
    raw=urllib.request.urlopen(url,timeout=25).read()
    data=json.loads(raw);columns=[m['name'] for m in data['metadata']]
    assert not set(columns)&{'rvz_radvel','rvz_redshift','rvz_err','velvalue','meanerror'}
    return {'query':adql,'response_sha256':hashlib.sha256(raw).hexdigest(),
            'records':[dict(zip(columns,row)) for row in data['data']]}


def fetch(target):
    ra,dec=target['ra_deg'],target['dec_deg']
    adql=f"SELECT TOP 100 b.oid, b.main_id, b.ra, b.dec, b.otype, b.rvz_type, b.rvz_nature, b.rvz_qual, b.rvz_wavelength, b.rvz_bibcode, i.ids FROM basic AS b LEFT JOIN ids AS i ON b.oid=i.oidref WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{ra},{dec},0.008333333333333333))"
    result={'target_name':target['target_name'],'query_center_deg':[ra,dec]}
    try:
        result['counterparts']=query(adql)
        result['truncated']=len(result['counterparts']['records'])>=100
        result['measurement_metadata']=[]
        for counterpart in result['counterparts']['records']:
            if counterpart.get('rvz_bibcode'):
                oid=int(counterpart['oid'])
                detail=query(f'SELECT oidref, veltype, nature, wdomain, bibcode, origin, remarks FROM mesVelocities WHERE oidref={oid}')
                result['measurement_metadata'].append(detail)
    except Exception as exc:
        result['error']=str(exc)
    return result


if __name__=='__main__':
    features=json.loads((HERE/'elves-field-feature-audit.json').read_text())
    decisions=json.loads((HERE/'elves-identity-decisions.json').read_text())
    pending={r['target_name'] for r in decisions['records'] if r['decision']=='pending_freshness_and_host_audit'}
    targets=[r for r in features['records'] if r['target_name'] in pending]
    assert len(targets)==26
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records=list(pool.map(fetch,targets))
    out={'endpoint':ENDPOINT,'protocol_sha256':hashlib.sha256((HERE/'elves-velocity-metadata-protocol.md').read_bytes()).hexdigest(),
         'feature_sha256':hashlib.sha256((HERE/'elves-field-feature-audit.json').read_bytes()).hexdigest(),
         'targets':26,'request_errors':sum('error' in r for r in records),
         'targets_with_counterparts':sum(bool(r.get('counterparts',{}).get('records')) for r in records),
         'status':'Current source metadata only; no target outcome columns queried and no automatic counterpart selection.',
         'records':records}
    (HERE/'elves-velocity-metadata.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k!='records'},indent=2))
    for target in records:
        print(target['target_name'],[(r['main_id'],r['otype'],r['rvz_type'],r['rvz_bibcode']) for r in target.get('counterparts',{}).get('records',[])])

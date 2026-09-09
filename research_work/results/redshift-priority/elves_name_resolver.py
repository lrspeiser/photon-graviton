"""Resolve current pending names; retain only a strict metadata projection."""
import concurrent.futures
import hashlib
import json
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

from elves_radio_quality import separation

HERE = Path(__file__).resolve().parent
ALLOWED = {'oname', 'alias', 'otype', 'jradeg', 'jdedeg'}


def fetch(name):
    url = 'https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/SNV?' + urllib.parse.quote(name)
    out = {'target_name': name, 'url': url}
    try:
        raw = urllib.request.urlopen(url, timeout=25).read()
        root = ET.fromstring(raw)
        out['response_sha256'] = hashlib.sha256(raw).hexdigest()
        out['returned_tag_names'] = sorted({e.tag for e in root.iter()})
        out['resolvers'] = [
            {'resolver_name': resolver.attrib.get('name'),
             'metadata': [{'tag': e.tag, 'value': e.text}
                          for e in resolver.iter() if e.tag in ALLOWED]}
            for resolver in root.iter('Resolver')]
    except Exception as exc:
        out['error'] = f'{type(exc).__name__}: {exc}'
    return out


if __name__ == '__main__':
    inputs = ['elves-quarantine-screen.json', 'elves-field-feature-audit.json',
              'elves-velocity-metadata.json', 'elves-name-resolver-protocol.md']
    blobs = {n: (HERE/n).read_bytes() for n in inputs}
    screen, features, counterparts = [json.loads(blobs[n]) for n in inputs[:3]]
    names = [r['target_name'] for r in screen['records']
             if r['decision'] == 'pending_freshness_and_host_audit']
    assert len(names) == 24 and len(set(names)) == 24
    features = {r['target_name']: r for r in features['records']}
    counterparts = {r['target_name']: r for r in counterparts['records']}
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(fetch, names))
    for record in records:
        name = record['target_name']
        for resolver in record.get('resolvers', []):
            fields = {m['tag']: m['value'] for m in resolver['metadata']}
            if 'jradeg' not in fields or 'jdedeg' not in fields:
                continue
            point = {'ra_deg': float(fields['jradeg']), 'dec_deg': float(fields['jdedeg'])}
            resolver['publisher_separation_arcsec'] = separation(point, features[name])
            resolver['prior_cone_comparisons'] = sorted([
                {'oid': c['oid'], 'main_id': c['main_id'], 'otype': c['otype'],
                 'separation_arcsec': separation(point, {'ra_deg': c['ra'], 'dec_deg': c['dec']})}
                for c in counterparts[name]['counterparts']['records']],
                key=lambda x: x['separation_arcsec'])
    out = {'targets': len(records), 'request_errors': sum('error' in r for r in records),
           'targets_with_resolver_records': sum(bool(r.get('resolvers')) for r in records),
           'input_hashes': {n: hashlib.sha256(b).hexdigest() for n,b in blobs.items()},
           'status': 'metadata-only name-resolution evidence; eligibility unchanged', 'records': records}
    (HERE/'elves-name-resolver.json').write_text(json.dumps(out, indent=2)+'\n', newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('records','input_hashes')}, indent=2))
    for r in records:
        print(r['target_name'], [(s['resolver_name'], s['metadata'], s.get('publisher_separation_arcsec'))
                                 for s in r.get('resolvers', [])], r.get('error', ''))

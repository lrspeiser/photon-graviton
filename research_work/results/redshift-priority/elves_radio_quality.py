"""Identifier, position and flag-only join to Yu 2022; no outcomes parsed."""
import argparse
import hashlib
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE_SHA = '816431e3254754d86fc2fccba6fd76b0347d0357f959b57bbfa34ae9d4283225'
SLICES = ((4, 10), (11, 20), (21, 29), (149, 158))


def parse(line):
    return {'agc': int(line[4:10]), 'ra_deg': float(line[11:20]),
            'dec_deg': float(line[21:29]), 'notes': line[149:158].strip()}


def separation(a, b):
    r1, d1, r2, d2 = map(math.radians, (a['ra_deg'], a['dec_deg'], b['ra_deg'], b['dec_deg']))
    h = math.sin((d2-d1)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r2-r1)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1, max(0, h)))))*3600


def run(raw):
    assert hashlib.sha256(raw).hexdigest() == SOURCE_SHA
    lines = raw.decode('ascii').splitlines()
    assert len(lines) == 29958
    radio = []
    for line in lines:
        padded = line.ljust(158)
        masked = ['X']*158
        for lo, hi in SLICES:
            masked[lo:hi] = padded[lo:hi]
        row = parse(padded)
        assert row == parse(''.join(masked))
        assert set(row['notes']) <= set('12rpc ,'), 'Unrecognized flags'
        radio.append(row)
    assert len({r['agc'] for r in radio}) == len(radio)
    inputs = ['elves-field-feature-audit.json', 'elves-velocity-metadata.json',
              'elves-velocity-metadata-summary.json', 'elves-radio-quality-protocol.md']
    blobs = {name: (HERE/name).read_bytes() for name in inputs}
    features = {r['target_name']: r for r in json.loads(blobs[inputs[0]])['records']}
    metadata = {r['target_name']: r for r in json.loads(blobs[inputs[1]])['records']}
    identities = json.loads(blobs[inputs[2]])['records']
    records = []
    for identity in identities:
        name = identity['target_name']
        ids = set(identity['identity_supported_oid'])
        aliases = [c['ids'] or '' for c in metadata[name]['counterparts']['records']
                   if len(ids) == 1 and c['oid'] in ids]
        agcs = set()
        for alias in aliases:
            for kind, number in re.findall(r'\b(AGC|UGC)\s+(\d+)\b', alias):
                number = int(number)
                if kind == 'AGC' or number < 100000:
                    agcs.add(number)
        matches = []
        for r in radio:
            angle = separation(features[name], r)
            if angle <= 30:
                matches.append({**r, 'separation_arcsec': angle,
                                'explicit_supported_alias': r['agc'] in agcs,
                                'overlapping_hi_flag': 'c' in r['notes'],
                                'projected_unconfused_neighbor_flag': 'p' in r['notes'],
                                'edge_low_weight_flag': 'r' in r['notes'],
                                'code2_optical_redshift_association': '2' in r['notes']})
        records.append({'target_name': name, 'supported_agc_aliases': sorted(agcs),
                        'matches': matches,
                        'status': 'quality metadata only; adopted ELVES measurement unresolved'})
    return {'source_url': 'https://cdsarc.cds.unistra.fr/ftp/J/ApJS/261/21/table2.dat.gz',
            'uncompressed_source_sha256': SOURCE_SHA,
            'source_rows': len(radio), 'masked_outcome_invariance_rows': len(radio),
            'input_hashes': {n: hashlib.sha256(b).hexdigest() for n, b in blobs.items()},
            'targets': len(records),
            'targets_with_positional_match': sum(bool(r['matches']) for r in records),
            'targets_with_explicit_supported_match': sum(any(m['explicit_supported_alias']
                                                           for m in r['matches']) for r in records),
            'records': records}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('raw_table', type=Path)
    args = parser.parse_args()
    out = run(args.raw_table.read_bytes())
    (HERE/'elves-radio-quality.json').write_text(json.dumps(out, indent=2)+'\n', newline='\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'records'}, indent=2))
    for r in out['records']:
        print(r['target_name'], [(m['agc'], m['notes'], round(m['separation_arcsec'], 2),
                                  m['explicit_supported_alias']) for m in r['matches']])

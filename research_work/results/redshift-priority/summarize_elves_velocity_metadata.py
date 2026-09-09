"""Replay counterpart support and source references without requesting outcomes."""
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASIC = {'oid', 'main_id', 'ra', 'dec', 'otype', 'rvz_type', 'rvz_nature',
         'rvz_qual', 'rvz_wavelength', 'rvz_bibcode', 'ids'}
MEASUREMENT = {'oidref', 'veltype', 'nature', 'wdomain', 'bibcode', 'origin', 'remarks'}


def normalized(value):
    return re.sub(r'[^a-z0-9]', '', value.lower())


def summarize(metadata, decisions):
    known = {r['target_name']: set(r['confirmed_pgc_aliases'])
             for r in decisions['records']}
    pending = {r['target_name'] for r in decisions['records']
               if r['decision'] == 'pending_freshness_and_host_audit'}
    names = [r['target_name'] for r in metadata['records']]
    assert len(names) == len(set(names)) and set(names) == pending
    rows = []
    for target in metadata['records']:
        # Failure or truncation must be resolved, never silently dropped.
        assert 'error' not in target and not target['truncated']
        counterparts = target['counterparts']['records']
        for c in counterparts:
            assert set(c) == BASIC, 'Unexpected basic columns'
        oids = {c['oid'] for c in counterparts}
        measurement_rows = [m for bundle in target['measurement_metadata']
                            for m in bundle['records']]
        for m in measurement_rows:
            assert set(m) == MEASUREMENT, 'Unexpected measurement columns'
            assert m['oidref'] in oids
        supported = []
        evidence = []
        name = target['target_name']
        for c in counterparts:
            aliases = c.get('ids') or ''
            pgcs = {int(p) for p in re.findall(r'(?:LEDA|PGC)\s+(\d+)', aliases)}
            shared = sorted(pgcs & known[name])
            direct = [a for a in aliases.split('|')
                      if normalized(re.sub(r'^\[[^\]]+\]\s*', '', a)) == normalized(name)]
            if shared or direct:
                supported.append(c)
                evidence.append({'oid': c['oid'], 'shared_confirmed_pgc': shared,
                                 'literal_name_aliases': direct})
        supported_oids = {c['oid'] for c in supported}
        references = sorted({m['bibcode'] for m in measurement_rows
                             if m['oidref'] in supported_oids and m['bibcode']})
        rows.append({
            'target_name': name,
            'counterpart_count': len(counterparts),
            'identity_supported_oid': [c['oid'] for c in supported],
            'identity_evidence': evidence,
            'identity_supported_basic_metadata': [
                {k: c[k] for k in ('main_id', 'rvz_type', 'rvz_nature',
                                   'rvz_bibcode', 'rvz_qual', 'rvz_wavelength')}
                for c in supported],
            'supported_measurement_references': references,
            'status': ('supported counterpart; adopted-table measurement still unresolved'
                       if len(supported) == 1 else 'counterpart identity unresolved or ambiguous')})
    preferred = Counter(c['rvz_bibcode'] or 'NO_BASIC_REFERENCE'
                        for r in rows if len(r['identity_supported_oid']) == 1
                        for c in r['identity_supported_basic_metadata'])
    return {'targets': len(rows),
            'targets_with_one_identity_supported_counterpart': sum(
                len(r['identity_supported_oid']) == 1 for r in rows),
            'multi_object_cones': sum(r['counterpart_count'] > 1 for r in rows),
            'preferred_reference_counts_for_supported_counterparts': dict(sorted(preferred.items())),
            'measurement_metadata_rows': sum(len(m['records']) for t in metadata['records']
                                             for m in t['measurement_metadata']),
            'adopted_measurements_verified': 0,
            'certified_fresh_targets': 0,
            'records': rows}


if __name__ == '__main__':
    inputs = ['elves-velocity-metadata.json', 'elves-identity-decisions.json']
    blobs = [(HERE / name).read_bytes() for name in inputs]
    output = summarize(*(json.loads(blob) for blob in blobs))
    output['input_hashes'] = {name: hashlib.sha256(blob).hexdigest()
                              for name, blob in zip(inputs, blobs)}
    (HERE / 'elves-velocity-metadata-summary.json').write_text(
        json.dumps(output, indent=2) + '\n', newline='\n')
    print(json.dumps({k: v for k, v in output.items() if k != 'records'}, indent=2))

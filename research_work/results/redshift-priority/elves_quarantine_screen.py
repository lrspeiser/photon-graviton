"""Conservative stage-2 screen; no candidate outcome values are read."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
COMMIT = '3e71904'


def norm(value):
    return re.sub(r'[^a-z0-9]', '', value.lower())


def main():
    names = ['elves-identity-decisions.json', 'elves-velocity-metadata-summary.json',
             'elves-velocity-metadata.json', 'elves-search-exposure.json',
             'elves-quarantine-screen-protocol.md']
    blobs = {n: (HERE/n).read_bytes() for n in names}
    initial, support, metadata, exposure = (json.loads(blobs[n]) for n in names[:4])
    assert len(initial['records']) == 29
    assert len({r['target_name'] for r in initial['records']}) == 29
    pending = {r['target_name'] for r in initial['records']
               if r['decision'] == 'pending_freshness_and_host_audit'}
    assert {r['target_name'] for r in support['records']} == pending
    assert {r['target_name'] for r in metadata['records']} == pending
    metadata = {r['target_name']: r for r in metadata['records']}
    patterns = {r['target_name']: set() for r in initial['records']}
    for target in support['records']:
        oids = target['identity_supported_oid']
        if len(oids) != 1:
            continue
        counterpart = next(c for c in metadata[target['target_name']]['counterparts']['records']
                           if c['oid'] == oids[0])
        for alias in (counterpart['ids'] or '').split('|'):
            alias = re.sub(r'^\[[^\]]+\]\s*', '', alias).strip()
            if re.match(r'^(NGC|UGC|IC|DDO|AGC|KKH|KDG|KKSG|LV |dw\d|Sex|Sextans|LEDA|PGC)', alias, re.I):
                patterns[target['target_name']].add(norm(alias))
    paths = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', COMMIT], cwd=ROOT, text=True).splitlines()
    paths = [p for p in paths if Path(p).suffix.lower() in {'.md','.py','.csv','.json','.txt','.dat','.tex','.bib'}]
    # One batch process; parse declared byte lengths rather than delimiters inside files.
    output = subprocess.check_output(['git', 'cat-file', '--batch'], cwd=ROOT,
                                     input=''.join(f'{COMMIT}:{p}\n' for p in paths).encode())
    offset = 0
    hits = {n: [] for n in patterns}
    count = 0
    for path in paths:
        end = output.index(b'\n', offset)
        header = output[offset:end].split()
        assert header[1] == b'blob'
        size = int(header[2]); raw = output[end+1:end+1+size]
        offset = end+1+size+1
        try:
            text = raw.decode('utf-8-sig')
        except UnicodeDecodeError:
            continue
        count += 1; compact = norm(text)
        for target, aliases in patterns.items():
            found = sorted(a for a in aliases if len(a) >= 5 and a in compact)
            if found:
                hits[target].append({'source': path, 'sha256': hashlib.sha256(raw).hexdigest(),
                                     'matched_normalized_aliases': found})
    assert offset == len(output)
    exposed = {norm(n) for n in exposure['target_names_displayed']}
    records = []
    for old in initial['records']:
        target = old['target_name']; reasons = list(old['reasons'])
        possible_overlap = [p for p in old['positional_candidates']
                            if p['historical_object_hit'] or p['historical_group_hit']]
        assert all(0 <= p['separation_arcsec'] <= 30 for p in possible_overlap)
        if possible_overlap:
            reasons.append('precautionary exclusion: positional coincidence with historical object/group; identity may remain unconfirmed')
        if norm(target) in exposed:
            reasons.append('target distance displayed in search before final freeze; no pair scored')
        if hits[target]:
            reasons.append('supported catalog alias found in historical tracked source; conservative source-availability exclusion')
        records.append({'target_name': target, 'initial_decision': old['decision'],
                        'decision': 'exclude_from_fresh_sample' if reasons else 'pending_freshness_and_host_audit',
                        'reasons': reasons, 'positional_overlap_evidence': possible_overlap,
                        'expanded_historical_hits': hits[target]})
    out = {'stage': 2, 'historical_commit': COMMIT, 'historical_text_files_scanned': count,
           'candidate_count': len(records), 'excluded_count': sum(bool(r['reasons']) for r in records),
           'remaining_pending': sum(not r['reasons'] for r in records), 'certified_fresh_targets': 0,
           'newly_excluded_names': [r['target_name'] for r in records if r['reasons'] and r['initial_decision'] != 'exclude_from_fresh_sample'],
           'input_hashes': {n: hashlib.sha256(b).hexdigest() for n,b in blobs.items()},
           'expanded_alias_patterns': {n: sorted(a) for n,a in patterns.items()}, 'records': records}
    (HERE/'elves-quarantine-screen.json').write_text(json.dumps(out, indent=2)+'\n', newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('records','expanded_alias_patterns','input_hashes')}, indent=2))
    print('historical hits', {n: [h['source'] for h in v] for n,v in hits.items() if v})


if __name__ == '__main__':
    main()

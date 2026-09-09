"""Identity-only audit. Never reads CF4 outcome columns or declares a blind sample."""
import csv
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUTPUTS = {'object-exposure-registry.json', 'object-exposure-summary.json'}


def normalize(key, value):
    key = key.lower()
    value = str(value).strip()
    if key in ('pgc', 'group_pgc', '1pgc', 'excluded_group_pgc'):
        if re.fullmatch(r'\d+(?:\.0+)?', value) and int(float(value)) > 0:
            return ('group_pgc' if key != 'pgc' else 'pgc', str(int(float(value))))
    elif key == 'cid' and re.fullmatch(r'\d+', value):
        return ('des_cid', str(int(value)))
    elif key in ('name', 'galaxy'):
        compact = re.sub(r'\s+', '', value).upper()
        m = re.fullmatch(r'(NGC|UGC|IC|PGC|LEDA)0*(\d+)([A-Z]?)', compact)
        if m:
            prefix, number, suffix = m.groups()
            if prefix in ('PGC', 'LEDA') and not suffix:
                return ('pgc', str(int(number)))
            return ('galaxy_name', prefix + str(int(number)) + suffix)
        if key == 'galaxy' or re.fullmatch(r'CGCG\d+-\d+', compact):
            if compact and compact not in ('NONE', 'NAN'):
                return ('galaxy_name', compact)
    return None


def main():
    paths = subprocess.check_output(['git', 'ls-files'], cwd=ROOT, text=True).splitlines()
    registry = defaultdict(set)
    sources = {}
    errors = []
    keys = {'pgc', 'group_pgc', '1pgc', 'excluded_group_pgc', 'cid', 'galaxy', 'name'}

    for rel in paths:
        path = ROOT / rel
        if path.suffix not in ('.csv', '.json') or path.name in OUTPUTS or path.name == 'maser-aliases.json':
            continue
        entries = set()

        def accept(key, value):
            if key.lower() not in keys:
                return
            if isinstance(value, list):
                for item in value:
                    accept(key, item)
            elif isinstance(value, (str, int, float)):
                item = normalize(key, value)
                if item:
                    entries.add(item)

        def walk(value):
            if isinstance(value, dict):
                for key, child in value.items():
                    accept(key, child)
                    walk(child)
            elif isinstance(value, list):
                for child in value:
                    walk(child)

        try:
            if path.suffix == '.csv':
                with path.open(encoding='utf-8-sig', newline='') as handle:
                    for row in csv.DictReader(handle):
                        for key, value in row.items():
                            if key:
                                accept(key, value)
            else:
                walk(json.loads(path.read_text(encoding='utf-8-sig')))
        except Exception as exc:
            errors.append({'source': rel, 'error': str(exc)})
        if entries:
            sources[rel] = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                            'identities': len(entries),
                            'status': 'available in tracked research; conservatively exclude, not proof every label was inspected'}
            for item in entries:
                registry[item].add(rel)

    # Fixed-width ID fields only. No velocities or distance labels are parsed.
    table = ROOT / 'redshift_paper_sources/time-redshift-expanded/table2.dat'
    groups = defaultdict(set)
    for line in table.read_text().splitlines():
        a, b = line[:7].strip(), line[8:15].strip()
        if a.isdigit() and b.isdigit() and int(b) > 0:
            groups[str(int(a))].add(str(int(b)))
    excluded_groups = {ident for (namespace, ident) in registry if namespace == 'group_pgc'}
    for namespace, ident in registry:
        if namespace == 'pgc':
            excluded_groups.update(groups.get(ident, set()))
    aliases_path = HERE / 'maser-aliases.json'
    aliases = json.loads(aliases_path.read_text()) if aliases_path.exists() else []
    maser = []
    for row in aliases:
        pgcs = row['pgc_aliases']
        matched_groups = sorted({g for p in pgcs for g in groups.get(str(p), set())}, key=int)
        hits = sorted({s for p in pgcs for s in registry.get(('pgc', str(p)), set())})
        maser.append({'name': row['name'], 'pgc_aliases': pgcs,
                      'cf4_groups': matched_groups, 'pgc_source_hits': hits,
                      'previously_excluded_groups': [g for g in matched_groups if g in excluded_groups],
                      'status': 'six-maser labels already evaluated; never fresh regardless of historical alias hits'})
        excluded_groups.update(matched_groups)
    records = [{'namespace': ns, 'id': ident, 'sources': sorted(files),
                'cf4_groups': sorted(groups.get(ident, set()), key=int) if ns == 'pgc' else []}
               for (ns, ident), files in sorted(registry.items())]
    out = {'scope': 'Tracked CSV and JSON recognized identity fields, plus CF4 ID-only group joins. Not a complete historical-exposure certification.',
           'limitations': ['Unstructured documents, PDFs, archives and unnamed arrays not identity-scanned.',
                          'General galaxy-name to PGC and supernova-host aliases remain unresolved.',
                          'Full raw catalog availability does not establish per-row label exposure.',
                          'Absence from this registry never certifies freshness.'],
           'group_source': {'path': table.relative_to(ROOT).as_posix(), 'sha256': hashlib.sha256(table.read_bytes()).hexdigest()},
           'sources': sources, 'records': records, 'conservative_excluded_group_pgc': sorted(map(int, excluded_groups)),
           'maser_overlap': maser, 'errors': errors}
    (HERE / 'object-exposure-registry.json').write_text(json.dumps(out, indent=2) + '\n', newline='\n')
    counts = {ns: sum(r['namespace'] == ns for r in records) for ns in sorted({r['namespace'] for r in records})}
    summary = {'source_count': len(sources), 'identity_counts': counts,
               'excluded_group_count': len(excluded_groups), 'errors': errors, 'maser_overlap': maser}
    (HERE / 'object-exposure-summary.json').write_text(json.dumps(summary, indent=2) + '\n', newline='\n')
    print(json.dumps(summary, indent=2))
    assert not errors, errors


if __name__ == '__main__':
    main()

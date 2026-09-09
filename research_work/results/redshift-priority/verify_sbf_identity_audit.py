"""Replay the saved source-identifier audit without downloading outcome data."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
audit=json.loads((HERE/'sbf-2025-identity-audit.json').read_text())
assert hashlib.sha256((HERE/'sbf-2025-audit-protocol.md').read_bytes()).hexdigest()==audit['protocol_sha256']
# Preserve the exact pre-audit exclusion snapshot, rather than letting new names
# inserted after this audit manufacture an apparent historical overlap.
registry_bytes=subprocess.check_output(['git','show','539fe59:research_work/results/redshift-priority/object-exposure-registry.json'],cwd=ROOT)
assert hashlib.sha256(registry_bytes).hexdigest()==audit['registry_sha256']
registry=json.loads(registry_bytes)
objects={r['id'] for r in registry['records'] if r['namespace']=='pgc'}
excluded=set(registry['conservative_excluded_group_pgc'])
source=ROOT/registry['group_source']['path']
assert hashlib.sha256(source.read_bytes()).hexdigest()==registry['group_source']['sha256']
groups={}
for line in source.read_text().splitlines():
    a,b=line[:7].strip(),line[8:15].strip()
    if a.isdigit() and b.isdigit():
        groups[int(a)]=int(b)
assert len(audit['rows'])==16
for row in audit['rows']:
    pgcs=sorted(set(map(int,re.findall(r'\b(?:LEDA|PGC)\s+(\d+)\b',row['identifier_section']))))
    assert pgcs==row['pgc_aliases'] and pgcs
    assert [p for p in pgcs if str(p) in objects]==row['prior_pgc_hits']
    assert sorted({groups[p] for p in pgcs if p in groups and groups[p]>0})==row['cf4_groups']
    assert sorted({groups[p] for p in pgcs if p in groups and groups[p] in excluded})==row['prior_group_hits']
assert sum(bool(r['prior_pgc_hits']) for r in audit['rows'])==audit['prior_object_overlap']==1
assert sum(bool(r['prior_group_hits']) for r in audit['rows'])==audit['prior_group_overlap']==16
print('Verified 16 primary-source aliases and historical group matches against the frozen pre-audit registry.')

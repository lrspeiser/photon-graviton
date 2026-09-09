"""Reproduce conservative decisions from saved identity evidence, not outcomes."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
initial=json.loads((HERE/'elves-identity-audit.json').read_text())
supplement=json.loads((HERE/'elves-supplemental-aliases.json').read_text())
historical=json.loads((HERE/'elves-historical-name-scan.json').read_text())
features=json.loads((HERE/'elves-field-feature-audit.json').read_text())
assert initial['features_sha256']==hashlib.sha256((HERE/'elves-field-feature-audit.json').read_bytes()).hexdigest()
raw=subprocess.check_output(['git','show',initial['historical_registry_commit']+':research_work/results/redshift-priority/object-exposure-registry.json'],cwd=ROOT)
assert hashlib.sha256(raw).hexdigest()==initial['historical_registry_sha256']
registry=json.loads(raw);excluded=set(registry['conservative_excluded_group_pgc'])
table=ROOT/initial['cf4_source']['path']
assert hashlib.sha256(table.read_bytes()).hexdigest()==initial['cf4_source']['sha256']
groups={}
for line in table.read_text().splitlines():
    a,b=line[:7].strip(),line[8:15].strip()
    if a.isdigit() and b.isdigit():groups[int(a)]=int(b)


def norm(name):return re.sub(r'[^a-z0-9]','',name.lower())


result=[]
for row in initial['records']:
    name=row['target_name'];pgcs=set(row['pgc_aliases'])
    accepted=[]
    for evidence in supplement:
        if evidence['target_name']!=name:continue
        # Require the staged target's literal name in the returned identifiers.
        # Positional proximity or similar alternate catalog names alone do not pass.
        ids=evidence.get('identifier_section','')
        if norm(name) in norm(ids) and evidence['pgc_aliases']:
            parsed=set(map(int,re.findall(r'\b(?:LEDA|PGC)\s+(\d+)\b',ids)))
            assert parsed==set(evidence['pgc_aliases'])
            pgcs.update(parsed);accepted.append(evidence['query'])
    group_hits=sorted({groups[p] for p in pgcs if p in groups and groups[p] in excluded})
    reasons=[]
    if group_hits:reasons.append('confirmed alias shares previously excluded group')
    if row['historical_pgc_hits'] or row['direct_name_hit']:reasons.append('known historical object')
    if historical['hits'][name]:reasons.append('historical source availability under a named alias; conservative exclusion')
    result.append({'target_name':name,'confirmed_pgc_aliases':sorted(pgcs),'accepted_supplemental_queries':accepted,
        'excluded_group_hits':group_hits,'historical_name_source_hits':historical['hits'][name],
        'decision':'exclude_from_fresh_sample' if reasons else 'pending_freshness_and_host_audit',
        'reasons':reasons,'positional_candidates':row['cf4_position_candidates']})
assert len(result)==29
assert {r['target_name'] for r in result}=={r['target_name'] for r in features['records'] if r['provisional_eligible']}
excluded_names=[r['target_name'] for r in result if r['reasons']]
assert set(excluded_names)=={'UGC05797','dw1046p1244','NGC4592'}
out={'candidate_count':29,'excluded_count':len(excluded_names),'remaining_pending':29-len(excluded_names),
    'historical_text_files_scanned':historical['text_files_scanned'],
    'confirmed_pgc_targets':sum(bool(r['confirmed_pgc_aliases']) for r in result),
    'certified_fresh_targets':0,'excluded_names':excluded_names,
    'input_hashes':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in ['elves-identity-audit.json','elves-supplemental-aliases.json','elves-historical-name-scan.json']},'records':result}
(HERE/'elves-identity-decisions.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('records','input_hashes')},indent=2))

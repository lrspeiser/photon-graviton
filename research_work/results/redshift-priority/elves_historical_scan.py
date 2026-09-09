"""Conservative historical name scan; source availability is not label exposure."""
import hashlib
import json
import re
import subprocess
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
COMMIT='3e71904'  # Before ELVES source acquisition and feature staging.
audit=json.loads((HERE/'elves-identity-audit.json').read_text())
supplement=json.loads((HERE/'elves-supplemental-aliases.json').read_text())


def compact(value):
    return re.sub(r'[^a-z0-9]','',value.lower())


patterns={}
for row in audit['records']:
    name=row['target_name']
    queries={compact(name)}
    evidence=[row]+[x for x in supplement if x['target_name']==name]
    for record in evidence:
        for pgc in record.get('pgc_aliases',[]):
            queries.update((f'pgc{pgc}',f'leda{pgc}'))
        # Preserve named catalog aliases, not numeric substrings or bibliography.
        ids=record.get('identifier_section','').splitlines()[1:]
        for line in ids:
            for item in re.split(r'\s{2,}',line.strip()):
                item=re.sub(r'^NAME\s+','',item)
                item=re.sub(r'^\[[^\]]+\]\s*','',item)
                if re.match(r'^(NGC|UGC|IC|DDO|AGC|KKH|KDG|KKSG|LV |dw\d|Sex|Sextans|Leo |LeG|M96-DF)',item,re.I):
                    queries.add(compact(item))
    patterns[name]=sorted(q for q in queries if len(q)>=5)

files=subprocess.check_output(['git','ls-tree','-r','--name-only',COMMIT],cwd=ROOT,text=True).splitlines()
hits={name:[] for name in patterns}
count=0
for rel in files:
    if Path(rel).suffix.lower() not in {'.md','.py','.csv','.json','.txt','.dat','.tex','.bib'}:
        continue
    raw=subprocess.check_output(['git','show',f'{COMMIT}:{rel}'],cwd=ROOT)
    try:text=raw.decode('utf-8-sig')
    except UnicodeDecodeError:continue
    normalized=compact(text);count+=1
    for name,queries in patterns.items():
        found=[q for q in queries if q in normalized]
        if found:hits[name].append({'source':rel,'sha256':hashlib.sha256(raw).hexdigest(),'matched_normalized_aliases':found})
result={'historical_commit':COMMIT,'text_files_scanned':count,
    'scope':'Conservative punctuation-insensitive alias substring scan of tracked historical text. Hits require interpretation; no-hit is not freshness certification. PDFs, archives, external conversations and coordinate-only matches are not covered.',
    'targets_with_hits':sum(bool(v) for v in hits.values()),'patterns':patterns,'hits':hits}
(HERE/'elves-historical-name-scan.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps({'files':count,'targets_with_hits':result['targets_with_hits'],'hit_names':[k for k,v in hits.items() if v]},indent=2))

"""Archive pinned KCWI radial products without fitting or transforming them."""
from pathlib import Path
import urllib.request,json,hashlib
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
REPO='TDCOSMO/TDCOSMO2025_public';COMMIT='d7f38db341f68be1df0d9ac1fc528c45113f94cf'
CACHE=ROOT/'research_work/generated/slacs-resolved-data';CACHE.mkdir(parents=True,exist_ok=True)
api='https://api.github.com/repos/'+REPO+'/git/trees/'+COMMIT+'?recursive=1'
tree=json.load(urllib.request.urlopen(api));assert not tree['truncated']
choices=[r for r in tree['tree'] if r['type']=='blob' and (r['path']=='ExternalLenses/SLACS/readme.md' or (r['path'].startswith('ExternalLenses/SLACS/slacs_kcwi_data/') and r['path'].endswith('.csv')))]
roles={r['Name']:r['role'] for r in json.loads((HERE.parent/'lensing-data-readiness/lens-observations-and-image-models.json').read_text(encoding='utf-8'))}
files=[];systems={}
for entry in choices:
 relative=Path(entry['path']);path=CACHE/relative
 assert path.resolve().is_relative_to(CACHE.resolve())
 url='https://raw.githubusercontent.com/'+REPO+'/'+COMMIT+'/'+entry['path']
 raw=urllib.request.urlopen(url,timeout=60).read()
 assert len(raw)==entry['size']
 assert hashlib.sha1(('blob '+str(len(raw))+'\0').encode()+raw).hexdigest()==entry['sha']
 path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
 files.append(dict(source_path=entry['path'],url=url,bytes=len(raw),git_blob_sha1=entry['sha'],sha256=hashlib.sha256(raw).hexdigest(),cache_path=str(path.relative_to(ROOT))))
 if relative.suffix=='.csv':
  identifier=relative.parent.name;name=identifier.removeprefix('SDSS')
  systems.setdefault(identifier,dict(Name=name,role=roles.get(name,'unassigned'),products=[]))['products'].append(relative.name)
assert len(systems)==14 and len(files)==29
assert all(len(r['products'])==2 for r in systems.values())
counts={role:sum(s['role']==role for s in systems.values()) for role in sorted(set(s['role'] for s in systems.values()))}
out=dict(repository=REPO,commit=COMMIT,scope='Archive and identifier join only. Radial values/covariance not fitted; no new observational score.',files=files,systems=list(systems.values()),role_counts=counts)
(HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(files=len(files),systems=len(systems),role_counts=counts,training=[r['Name'] for r in systems.values() if r['role']=='training'])))

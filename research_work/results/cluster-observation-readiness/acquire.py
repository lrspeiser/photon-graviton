"""Acquire actual Coma spectroscopy; no conversion of redshift to distance/mass."""
from pathlib import Path
import json,hashlib,gzip
from datetime import datetime,timezone
import requests
import numpy as np
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
CACHE=ROOT/'research_work/generated/coma-kang2025';CACHE.mkdir(parents=True,exist_ok=True)
base='https://cdsarc.cds.unistra.fr/ftp/J/ApJS/278/51/'
files=[]
for name in ['ReadMe','table2.dat','refs.dat']:
    path=CACHE/name
    url=base+('table2.dat.gz' if name=='table2.dat' else name)
    if not path.exists():
        r=requests.get(url,timeout=60);r.raise_for_status()
        assert b'<html' not in r.content[:1000].lower()
        path.write_bytes(gzip.decompress(r.content) if name=='table2.dat' else r.content)
    files.append(dict(name=name,url=url,bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),hash_applies_to='decompressed local file'))
rows=[]
for line in (CACHE/'table2.dat').read_text().splitlines():
    assert len(line)>=79
    rows.append(dict(seq=int(line[:5]),id=line[6:25].strip(),ra=float(line[26:36]),dec=float(line[37:46]),
        rmag=float(line[47:53]),point=int(line[54]),z=float(line[56:65]),ez=float(line[66:74]),reference=int(line[75:77]),member=int(line[78])))
assert len(rows)==37758
assert len(set(v['seq'] for v in rows))==len(rows)
assert all(np.isfinite([v['ra'],v['dec'],v['rmag'],v['z'],v['ez']]).all() for v in rows)
assert all(v['ez']>=0 and v['point'] in [0,1] and v['member'] in [0,1] for v in rows)
counts=dict(total=len(rows),unique_SDSS_IDs=len(set(v['id'] for v in rows)),point=sum(v['point'] for v in rows),extended=sum(1-v['point'] for v in rows),
    published_member_flag=sum(v['member'] for v in rows),zero_redshift_error=sum(v['ez']==0 for v in rows),this_work_reference=sum(v['reference']==1 for v in rows))
assert counts['point']==10863 and counts['extended']==26895 and counts['published_member_flag']==1826
result=dict(acquired_utc=datetime.now(timezone.utc).isoformat(),source='Kang et al. 2025, ApJS 278,51; CDS J/ApJS/278/51',files=files,counts=counts,
    sky_extent_deg=dict(ra=[min(v['ra'] for v in rows),max(v['ra'] for v in rows)],dec=[min(v['dec'] for v in rows),max(v['dec'] for v in rows)]),
    use='Exploratory acquisition and schema validation only; not a held-out validation set.',
    missing=['individual background shape/shear measurements and calibration','lensing covariance and masks','distances independent of the proposed propagation law','bolometric luminosities and source histories','baryonic spatial model','membership and selection audit'],
    prohibitions=['Do not treat published membership as model-independent.','Do not convert redshift into distance by an expansion law.','Do not treat zero redshift errors as exact measurements.','Do not fit against published halo mass as if it were measured shear.'])
(OUT/'acquisition.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))

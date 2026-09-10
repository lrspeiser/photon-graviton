"""Check the 72 provisional training seeds at their actual 2MASS epochs."""
from pathlib import Path
import hashlib
import io
import json
import re
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache/training-orbit-launches'
seeds=pd.read_parquet(CACHE/'launches.parquet')
ids=seeds.APOGEE_ID.str.removeprefix('2M').tolist()
assert len(ids)==72 and all(re.fullmatch(r'\d{8}[+-]\d{7}',s) for s in ids)
frames=[];manifest=[]
for start in range(0,len(ids),36):
    selected=sorted(ids)[start:start+36]
    query='SELECT designation,ra,dec,jdate,err_maj,err_min,err_ang FROM fp_psc WHERE designation IN ('+','.join("'"+s+"'" for s in selected)+')'
    path=HERE/f'2mass-seed-epochs-{start//36}.csv'
    if not path.exists():
        response=requests.post('https://irsa.ipac.caltech.edu/TAP/sync',data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query},timeout=45)
        response.raise_for_status()
        parsed=pd.read_csv(io.StringIO(response.text),dtype={'designation':str})
        assert set(parsed.designation)==set(selected) and parsed.designation.is_unique
        path.write_text(response.text,encoding='utf-8',newline='\n')
    frame=pd.read_csv(path,dtype={'designation':str});frames.append(frame)
    manifest.append(dict(query=query,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
joined=pd.concat(frames,ignore_index=True)
assert set(joined.designation)==set(ids) and joined.designation.is_unique
source=ROOT/'research_work/data-cache/stellar-association-window/training-screen.parquet'
astrometry=pd.read_parquet(source,columns=['source_id','ra','dec','pmra','pmdec'])
d=seeds.assign(designation=ids).merge(joined,on='designation',validate='one_to_one')
d=d.merge(astrometry,on='source_id',validate='one_to_one',suffixes=('_2mass','_gaia'))
epoch=Time(d.jdate.to_numpy(),format='jd').jyear
dt=epoch-2016.
gaia=SkyCoord(d.ra_gaia.to_numpy()*u.deg,d.dec_gaia.to_numpy()*u.deg)
propagated=gaia.spherical_offsets_by(d.pmra.to_numpy()*dt*u.mas,d.pmdec.to_numpy()*dt*u.mas)
measured=SkyCoord(d.ra_2mass.to_numpy()*u.deg,d.dec_2mass.to_numpy()*u.deg)
residual=propagated.separation(measured).arcsec
rows=[dict(source_id=str(row.source_id),designation=row.designation,epoch_jyear=float(e),
           separation_arcsec=float(s),position_review_flag=bool(s>.5),
           two_mass_major_error_arcsec=float(row.err_maj)) for row,e,s in zip(d.itertuples(),epoch,residual)]
result=dict(scope='Training-seed actual-epoch association check; not a formal match probability.',
            endpoint='https://irsa.ipac.caltech.edu/TAP/sync',catalog='fp_psc',requested=72,returned=len(d),
            separation_threshold_arcsec=.5,flagged=int((residual>.5).sum()),maximum_separation_arcsec=float(residual.max()),
            rows=rows,queries=manifest,holdouts_opened=False,
            approximation='Gaia mean proper motion propagated on the sphere; no perspective acceleration, annual parallax or binary orbit fit.',
            source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),source,CACHE/'launches.parquet']})
(HERE/'epoch-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print('Actual-epoch flags:',result['flagged'],'max separation:',result['maximum_separation_arcsec'],flush=True)

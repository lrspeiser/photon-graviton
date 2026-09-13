"""Acquire documented DR17 targeting products for a preselected training field."""
from pathlib import Path
import json,hashlib
from urllib.parse import urljoin
import requests,numpy as np,pandas as pd
from astropy.io import fits
H=Path(__file__).resolve().parent;R=H.parents[2];C=R/'research_work/data-cache/selection-fields';C.mkdir(exist_ok=True)
base='https://data.sdss.org/sas/dr17/apogee/target/'
files=['allField.fits','allDesign.fits','apogee2Plate.fits','apogee2Object/apogee2Object_300%2B00.fits']
manifest=[]
for name in files:
 url=urljoin(base,name);path=C/name.split('/')[-1].replace('%2B','+')
 if not path.exists():
  response=requests.get(url,timeout=(15,60));response.raise_for_status();assert response.content.startswith(b'SIMPLE')
  path.write_bytes(response.content)
 with fits.open(path,memmap=True) as hd:
  hd.verify('exception');extensions=[]
  for h in hd[1:]:
   if isinstance(h,fits.BinTableHDU):extensions.append(dict(name=h.name,rows=len(h.data),columns=h.columns.names,checksum_present='CHECKSUM' in h.header,checksum_valid=int(h.verify_checksum()) if 'CHECKSUM' in h.header else None))
 manifest.append(dict(url=url,path=str(path.relative_to(R)),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),extensions=extensions));print(json.dumps(dict(file=path.name,rows=[v['rows'] for v in extensions])),flush=True)
memberpath=R/'research_work/data-cache/stellar-selection-membership/training-membership.parquet';m=pd.read_parquet(memberpath);selected=m[m.status.eq('main_red')&m.FIELD.eq('300+00')&m.TELESCOPE.eq('lco25m')];assert len(selected)==521
with fits.open(C/'apogee2Object_300+00.fits') as hd:
 t=hd[1].data;ids=np.char.strip(np.asarray(t['APOGEE_ID']).astype(str));present=np.isin(selected.APOGEE_ID.to_numpy(dtype=str),ids)
 locids=np.unique(t['LOCATION_ID']).astype(int).tolist()
 info=dict(parent_rows=len(t),unique_parent_ids=len(np.unique(ids)),training_main_stars=len(selected),matched_training_ids=int(present.sum()),unmatched_training_ids=selected.loc[~present,'APOGEE_ID'].tolist(),location_ids=locids)
with fits.open(C/'allField.fits') as hd:
 t=hd[1].data;print('FIELD_COLUMNS '+str(t.names),flush=True)
 namecol=next(k for k in ['FIELD_NAME','FIELD'] if k in t.names);mask=np.char.strip(np.asarray(t[namecol]).astype(str))=='300+00';info['field_summary_rows']=int(mask.sum())
 info['field_metadata']=[{k:(v[k].tolist() if isinstance(v[k],np.ndarray) else v[k].item() if isinstance(v[k],np.generic) else v[k]) for k in t.names} for v in t[mask]]
with fits.open(C/'allDesign.fits') as hd:
 t=hd[1].data;selected_designs=t[np.isin(t['LOCATION_ID'],locids)]
 fields=['DESIGN_ID','LOCATION_ID','RADIUS','NUMBER_OF_VISITS','COHORT_SHORT_VERSION','COHORT_MEDIUM_VERSION','COHORT_LONG_VERSION','COHORT_FRACTION','COHORT_MIN_H','COHORT_MAX_H','NUMBER_OF_SELECTION_BINS','BIN_DEREDDENED_MIN_JK_COLOR','BIN_DEREDDENED_MAX_JK_COLOR','PLATERUN']
 info['designs']=[{k:(v[k].tolist() if isinstance(v[k],np.ndarray) else v[k].item() if isinstance(v[k],np.generic) else v[k]) for k in fields} for v in selected_designs]
result=dict(scope='One high-count training field selected by membership count; targeting metadata acquisition, not completeness or gravity inference',manifest=manifest,membership_sha256=hashlib.sha256(memberpath.read_bytes()).hexdigest(),pilot=info,documentation='https://www.sdss4.org/dr17/irspec/spectro_data/',final_test_kinematics_opened=False)
(H/'results.json').write_text(json.dumps(result,indent=2,default=str)+'\n',encoding='utf-8',newline='\n');print(json.dumps(info,default=str))

"""Recover targeting metadata for the existing training sample only."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from astropy.io import fits
H=Path(__file__).resolve().parent;R=H.parents[2];C=R/'research_work/data-cache'
paths=dict(training=C/'bulge-velocity-components/training-components.parquet',parent=C/'stellar-catalogs/matched-with-gaia-covariance.parquet',allstar=C/'stellar-catalogs/allStarLite-dr17-synspec_rev1.fits')
def text(v):return v.decode().strip() if isinstance(v,bytes) else str(v).strip()
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
training=pd.read_parquet(paths['training'],columns=['source_id','holdout_role','mean_R_kpc','mean_z_kpc']);assert len(training)==27884 and training.holdout_role.eq('training').all()
parent=pd.read_parquet(paths['parent'],columns=['source_id','APOGEE_ID','FIELD','TELESCOPE'],filters=[('source_id','in',training.source_id.tolist())]);assert parent.source_id.is_unique
keys=['APOGEE_ID','FIELD','TELESCOPE']
for k in keys:parent[k]=parent[k].map(text)
base=training.merge(parent,on='source_id',validate='one_to_one');assert len(base)==len(training)
columns=['EXTRATARG','TARGFLAGS','SURVEY','PROGRAMNAME','J','H','K','AK_TARG','AK_TARG_METHOD','MIN_H','MAX_H','MIN_JK','MAX_JK','APOGEE_TARGET1','APOGEE_TARGET2','APOGEE2_TARGET1','APOGEE2_TARGET2','APOGEE2_TARGET3']
with fits.open(paths['allstar'],memmap=True) as hd:
 table=hd[1].data;ids=np.char.strip(np.asarray(table['APOGEE_ID']).astype(str));use=np.isin(ids,base.APOGEE_ID.to_numpy(dtype=str));index=np.flatnonzero(use)
 values={k:[text(v) for v in table[k][index]] for k in keys}
 for k in columns:
  a=np.array(table[k][index]);values[k]=[text(v) for v in a] if a.dtype.kind in 'SUO' else a.astype(a.dtype.newbyteorder('='))
 values['allstar_row']=index
 candidates=base.merge(pd.DataFrame(values),on=keys,how='left',validate='one_to_many')
records=[]
for sid,g in candidates.groupby('source_id',sort=False):
 valid=g.EXTRATARG.dropna();status='unmatched' if len(valid)==0 else ('ambiguous_flags' if valid.nunique()>1 else ('main_red' if int(valid.iloc[0])==0 else 'non_main'))
 records.append(dict(source_id=int(sid),status=status,matches=len(valid),EXTRATARG=int(valid.iloc[0]) if len(valid)>0 and valid.nunique()==1 else None))
audit=base.merge(pd.DataFrame(records),on='source_id',validate='one_to_one');out=C/'stellar-selection-membership';out.mkdir(exist_ok=True);candidates.to_parquet(out/'training-targeting-candidates.parquet',index=False);audit.to_parquet(out/'training-membership.parquet',index=False)
regions={'all':np.ones(len(audit),bool),'bulge_plane':(audit.mean_R_kpc<3.5)&(audit.mean_z_kpc.abs()<.2),'bulge_offplane':(audit.mean_R_kpc<3.5)&audit.mean_z_kpc.abs().between(.5,1.5),'disk_plane':(audit.mean_R_kpc>=5)&(audit.mean_z_kpc.abs()<.2),'disk_offplane':(audit.mean_R_kpc>=5)&audit.mean_z_kpc.abs().between(.5,1.5)}
summary={name:dict(n=int(mask.sum()),status_counts=audit.loc[mask,'status'].value_counts().to_dict()) for name,mask in regions.items()}
result=dict(scope='Targeting-metadata recovery for exposed training stars; no selection probability or gravitational fit',inputs={str(p.relative_to(R)):sha(p) for p in paths.values()},outputs={str(p.relative_to(R)):sha(p) for p in out.glob('*.parquet')},matching_keys=keys,recovered_columns=columns,region_summary=summary,stars_with_multiple_matches=int((audit.matches>1).sum()),documentation=['https://www.sdss4.org/dr17/irspec/targets/selection-biases/','https://www.sdss4.org/dr17/irspec/targettingbits/'],final_test_opened=False)
(H/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(summary))

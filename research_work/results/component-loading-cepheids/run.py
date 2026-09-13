"""Fixed component-deposit prediction on previously exposed Cepheid bins."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent;R=H.parents[2]
source=H.parent/'component-attached-loading/results.json';d=json.loads(source.read_text(encoding='utf-8'))
for name,value in d['hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
eta=next(v['coefficient_by_component'] for v in d['cases'] if v['resolution']=='finer' and v['fit_subset']=='inner')
validation=H.parent/'cepheid-validation/results.json';prior=json.loads(validation.read_text(encoding='utf-8'));catalog=R/'research_work/data-cache/cepheid-stars/cepheid-validation-selected.parquet'
assert hashlib.sha256(catalog.read_bytes()).hexdigest()==prior['selected_catalog_sha256']
stars=pd.read_parquet(catalog);assert len(stars)==167 and set(stars.role)=={'validation'}
binpath=H.parent/'cepheid-validation/validation-bins.json';bins=json.loads(binpath.read_text(encoding='utf-8'));radius=stars.R_kpc.to_numpy();membership=stars.bin.to_numpy()
for row in bins:assert int((membership==row['bin']).sum())==row['n']
code=H.parent/'cepheid-validation/components.py';spec=importlib.util.spec_from_file_location('cep_components',code);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);f=mod.f
assert eta['nuclei']==eta['softened_centre']==0
outputs=[]
for refined in [False,True]:
 baryon=f.Baryons(refined);vbar2=-radius*f.force(baryon,radius,0)[:,0]
 bar=f.CachedAxisymmetric(R/'research_work/data-cache/bar-field/bar-L64.npz')
 added=eta['bar']*f.force(bar,radius,0)[:,0]
 for i in range(4):
  if eta['disk_'+str(i)]>0:added+=eta['disk_'+str(i)]*f.force(mod.DiskPart(refined,[i]),radius,0)[:,0]
 total=vbar2-radius*added;assert np.all(total>0)
 predicted=[]
 for row in bins:
  mask=membership==row['bin']
  if not row.get('moment_valid'):continue
  vc=float(np.sqrt(total[mask].mean()));ordinary=float(np.sqrt(vbar2[mask].mean()))
  if refined:assert abs(ordinary-row['original_ordinary_vc_kms'])<1e-6
  predicted.append(dict(bin=row['bin'],n=row['n'],radius_mean=row['R_mean_kpc'],inferred_jeans_speed=row['jeans_proxy_kms'],predicted_balance_speed=vc,ordinary_speed=ordinary,original_empirical_speed=row['original_completion_vc_kms'],residual_kms=vc-row['jeans_proxy_kms']))
 residual=np.array([v['residual_kms'] for v in predicted]);outputs.append(dict(refined=refined,rms_kms=float(np.sqrt(np.mean(residual**2))),bias_kms=float(residual.mean()),underpredicted_bins=int((residual<0).sum()),bins=predicted))
print(json.dumps([{k:v for k,v in o.items() if k!='bins'} for o in outputs]))
refinement=max(abs(a['predicted_balance_speed']-b['predicted_balance_speed']) for a,b in zip(outputs[0]['bins'],outputs[1]['bins']))
paths=[Path(__file__),source,validation,binpath,catalog,code,H.parent/'conservative-field-completion/run.py',R/'research_work/data-cache/bar-field/bar-L64.npz',R/'research_work/data-cache/bar-field/nuclei-L16.npz']
out=dict(scope='Fixed new source model on already exposed Cepheid validation subset; approximate Jeans comparison, not fresh validation',hashes={str(v.relative_to(R)):hashlib.sha256(v.read_bytes()).hexdigest() for v in paths},coefficients=eta,selected_stars=167,final_test_opened=False,parameters_fitted_to_these_velocities=False,maximum_refinement_kms=refinement,results=outputs)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

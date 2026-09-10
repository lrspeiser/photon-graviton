"""Evaluate frozen models on the reserved validation role, never the test role."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from components import Components, Model, f

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
protocol=json.loads((HERE/'protocol.json').read_text())
for name,digest in protocol['frozen_sha256'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
spec=importlib.util.spec_from_file_location('frozen_reduction',HERE.parent/'cepheid-common-frame/run.py')
red=importlib.util.module_from_spec(spec);spec.loader.exec_module(red)
CACHE=ROOT/'research_work/data-cache/cepheid-stars'
split=pd.read_parquet(CACHE/'cepheid-common-frame-split.parquet')
# Filtering in the parquet reader limits subsequent numerical work to reserved validation IDs.
ids=split.loc[split.role.eq('validation') & split.measurement_and_mode_eligible,'source_id'].tolist()
d=pd.read_parquet(CACHE/'gaia-dr3-dcep-parent-with-flags.parquet',filters=[('source_id','in',ids)])
assert len(d)==len(ids)==433 and d.source_id.is_unique
assert not set(d.source_id)&set(split.loc[split.role.eq('test'),'source_id'])
cal=red.P['calibration'];is_f=d.mode_best_classification.eq('FUNDAMENTAL')
period=np.where(is_f,d.pf,d.p1_o)
aa=np.where(is_f,cal['fundamental']['intercept'],cal['first_overtone']['intercept'])
bb=np.where(is_f,cal['fundamental']['slope'],cal['first_overtone']['slope'])
w=d.int_average_g-cal['wesenheit_color_coefficient']*(d.int_average_bp-d.int_average_rp)
d['distance_kpc']=10**((w-aa-bb*np.log10(period)-10)/5)
c=red.coords(d)
select=(c[:,0]>=6)&(c[:,0]<=18)&(abs(c[:,1])<=30)&(abs(c[:,2])<=.5)&(abs(c[:,5])<=100)
d=d.loc[select].copy();c=c[select]
membership=np.clip(np.floor(c[:,0]-6).astype(int),0,11)
variance=red.variances(d)
bins=red.summarize_bins(c,variance,membership)
for j,name in enumerate(['R_kpc','phi_deg','z_kpc','vR_kms','vphi_kms','vz_kms']):d[name]=c[:,j]
d['bin']=membership;d['role']='validation'
d.to_parquet(CACHE/'cepheid-validation-selected.parquet',index=False)
outputs={}
for refined in [False,True]:
    print('Building validation predictions, refined:',refined,flush=True)
    components=Components(refined)
    predictions={}
    for label,scales in [('original',np.ones(3)),('balanced',np.array(protocol['balanced_component_scales']))]:
        b=Model(components,scales);extra=f.Completion(b,refined)
        ordinary=-c[:,0]*f.force(b,c[:,0],0)[:,0]
        total=ordinary-c[:,0]*f.force(extra,c[:,0],0)[:,0]
        assert np.all(ordinary>0) and np.all(total>0)
        for suffix,values in [('ordinary',ordinary),('completion',total)]:
            predictions[label+'_'+suffix]=[float(np.sqrt(values[membership==r['bin']].mean())) if r.get('moment_valid') else None for r in bins]
    outputs['refined' if refined else 'coarse']=predictions
scores={};refinement={}
for model,prediction in outputs['refined'].items():
    valid=[(r,value) for r,value in zip(bins,prediction) if value is not None]
    residual=np.array([value-r['jeans_proxy_kms'] for r,value in valid])
    scores[model]=dict(bins=len(valid),rms_kms=float(np.sqrt(np.mean(residual**2))),
        bias_kms=float(residual.mean()),underpredicted_bins=int((residual<0).sum()),
        median_prediction_to_proxy=float(np.median([value/r['jeans_proxy_kms'] for r,value in valid])))
    refinement[model]=max(abs(a-b) for a,b in zip(outputs['coarse'][model],prediction) if a is not None)
    for row,value in zip(bins,prediction):row[model+'_vc_kms']=value
    assert refinement[model]<.1
result=dict(role='Reserved validation stars from an already examined Gaia catalog; conditional Jeans-pipeline test, not independent astrophysical evidence',
    pre_spatial_candidates=len(ids),selected_validation_stars=len(d),
    selected_modes=d.mode_best_classification.value_counts().to_dict(),
    evaluated_bins=[r['bin'] for r in bins if r.get('moment_valid')],
    unscored_bins=[dict(bin=r['bin'],n=r['n'],reason='Below five stars' if r['n']<5 else 'Nonphysical corrected moment under declared error model') for r in bins if not r.get('moment_valid')],
    scores=scores,max_coarse_refined_change_kms=refinement,
    primary_model='original_completion',secondary_model='balanced_completion',
    parameters_fitted_on_validation=False,final_test_outcomes_opened=False,
    validation_now_exposed=True,
    selected_catalog_sha256=hashlib.sha256((CACHE/'cepheid-validation-selected.parquet').read_bytes()).hexdigest(),
    protocol_sha256=hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
    limitations=protocol['limitations'])
for name,obj in [('results',result),('validation-bins',bins)]:
    (HERE/f'{name}.json').write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2),flush=True)

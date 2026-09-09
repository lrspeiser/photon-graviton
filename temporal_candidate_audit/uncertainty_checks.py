from pathlib import Path
import pandas as pd,numpy as np,json
P=Path(__file__).resolve().parent;O=P/'results'
d=pd.read_csv(O/'SPARC_predictions.csv');loss={}
for m,g in d.groupby('model'):
 loss[m]=g.assign(loss=np.log10(g.predicted_kms/g.observed_kms)**2).groupby('galaxy').loss.mean().sort_index().values
rng=np.random.default_rng(12345);ix=rng.integers(0,31,(5000,31));out={}
for a,b in [('density_p1','ordinary'),('density_free_p','ordinary'),('density_p1','power'),('density_free_p','power'),('power','RAR')]:
 v=np.sqrt(loss[a][ix].mean(axis=1))-np.sqrt(loss[b][ix].mean(axis=1));out[a+'_minus_'+b]={'delta_log_RMSE_dex':float(np.sqrt(loss[a].mean())-np.sqrt(loss[b].mean())),'descriptive_bootstrap95':np.quantile(v,[.025,.975]).tolist()}
(O/'SPARC_uncertainty.json').write_text(json.dumps(out,indent=2))
# Verify declared data sizes, disjoint identities, and finite predictions.
r=pd.read_csv(O/'redshift_predictions.csv');s=r[r.model=='exponential'];sp=json.loads((P/'data/sparc_frozen.json').read_text())['split']
checks={'redshift_164_unique_groups':s.group.nunique()==164,'redshift_split_counts':s.groupby('split').size().to_dict(),'redshift_no_group_overlap':s.groupby('group').split.nunique().max()==1,'redshift_no_sky_tile_overlap':s.groupby('tile').split.nunique().max()==1,'SPARC_no_galaxy_overlap':len(sum(sp.values(),[]))==len(set(sum(sp.values(),[]))),'all_numeric_predictions_finite':bool(np.isfinite(r[['predicted_cz','sigma_cz']]).all().all() and np.isfinite(d[['predicted_kms']]).all().all()),'SPARC_new_parameter_count':{'density_p1':2,'density_free_p':3},'test_status':'Reused historical splits, not new blind observations'}
(O/'verification.json').write_text(json.dumps(checks,indent=2,default=lambda x:bool(x)))
print(json.dumps({'bootstrap':out,'verification':checks},indent=2,default=lambda x:bool(x)))

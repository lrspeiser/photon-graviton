"""Exact global-mass homogeneity diagnostic for the frozen field completion.

All observations here were previously exposed. These are nuisance fits, not
new held-out evidence or a mass-prior posterior.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import minimize_scalar,brentq

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'conservative-field-completion'
INPUT=SOURCE/'predictions-refined.json'
PARAM=SOURCE/'results-refined.json'
raw=json.loads(INPUT.read_text());par=json.loads(PARAM.read_text())
P=par['parameters']['p']
data={}
for obs in ['rotation','vertical']:
    ordinary=[r for r in raw if r['observable']==obs and r['model']=='ordinary']
    total=[r for r in raw if r['observable']==obs and r['model']=='conservative_completion']
    assert len(ordinary)==len(total)
    assert [(r['R_kpc'],r['z_kpc'],r['observed']) for r in ordinary]==[(r['R_kpc'],r['z_kpc'],r['observed']) for r in total]
    power=2 if obs=='rotation' else 1
    b=np.array([r['predicted']**power for r in ordinary])
    c=np.array([r['predicted']**power for r in total])-b
    assert np.all(b>0) and np.all(c>0)
    data[obs]=dict(b=b,c=c,power=power,seen=np.array([r['observed'] for r in total]),rows=total)

def prediction(obs,lam):
    d=data[obs]
    return (lam*d['b']+lam**P*d['c'])**(1/d['power'])

def metrics(lam):
    result={}
    for obs,d in data.items():
        pred=prediction(obs,lam)
        result[obs]=dict(rmse=float(np.sqrt(np.mean((pred-d['seen'])**2))),
            bias=float(np.mean(pred-d['seen'])),median_ratio=float(np.median(pred/d['seen'])),
            log_mse=float(np.mean(np.log(pred/d['seen'])**2)))
    result['equal_observable_log_score']=.5*(result['rotation']['log_mse']+result['vertical']['log_mse'])
    return result

BOUNDS=(.3,2.)
cases=[]
for objective in ['rotation','vertical','joint']:
    def score(lam):
        m=metrics(lam)
        return m['equal_observable_log_score'] if objective=='joint' else m[objective]['rmse']**2
    fit=minimize_scalar(score,bounds=BOUNDS,method='bounded',options={'xatol':1e-12})
    assert fit.success
    grid=np.linspace(*BOUNDS,1701)
    assert fit.fun<=min(score(x) for x in grid)+1e-9
    cases.append(dict(calibration_objective=objective,lambda_mass=float(fit.x),metrics=metrics(fit.x)))
identity=metrics(1.)
for obs in data:
    assert abs(identity[obs]['rmse']-par['scores'][obs+'/conservative_completion']['rmse'])<1e-10

# One mass multiplier needed by each row; these values are diagnostics, never
# adopted as independent per-star or per-region fitting parameters.
row_roots=[]
for obs,d in data.items():
    for i,row in enumerate(d['rows']):
        root=brentq(lambda lam:prediction(obs,lam)[i]-d['seen'][i],.01,10.,xtol=1e-13)
        row_roots.append(dict(observable=obs,R_kpc=row['R_kpc'],z_kpc=row['z_kpc'],
            individually_required_lambda=root,adopted_as_model_parameter=False))
out=dict(classification='Exposed-data global ordinary-mass sensitivity; not fresh validation or full baryon uncertainty',
    power=P,lambda_search_bounds=list(BOUNDS),fixed_geometry=True,
    scaling='a_total(lambda)=lambda*a_b + lambda^p*a_extra',
    component_scope='All ordinary matter scaled together, including stars, gas and central component; no claim this is an allowed astrophysical prior',
    baseline=identity,cases=cases,
    median_individually_required_lambda={obs:float(np.median([r['individually_required_lambda'] for r in row_roots if r['observable']==obs])) for obs in data},
    input_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [INPUT,PARAM]},
    holdout_scores_opened=False)
rows=[]
for case in [dict(calibration_objective='original',lambda_mass=1.)]+cases:
    for obs,d in data.items():
        pred=prediction(obs,case['lambda_mass'])
        for i,row in enumerate(d['rows']):
            rows.append(dict(case=case['calibration_objective'],lambda_mass=case['lambda_mass'],
                observable=obs,R_kpc=row['R_kpc'],z_kpc=row['z_kpc'],observed=row['observed'],predicted=float(pred[i])))
profile=[dict(lambda_mass=float(lam),metrics=metrics(lam)) for lam in np.linspace(*BOUNDS,171)]
for name,obj in [('results.json',out),('predictions.json',rows),('row-mass-diagnostics.json',row_roots),('profile.json',profile)]:
    (HERE/name).write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
print(json.dumps(out,indent=2))

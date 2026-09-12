"""Diagnose a scatter-boundary fit without modifying the running batch."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.optimize import minimize
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import population_log_likelihood
p=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'));results=[]
for case in p['cases']:
    path=ROOT/'research_work/generated/timing-coverage-calibration'/(case['label']+'.npz')
    assert hashlib.sha256(path.read_bytes()).hexdigest()==case['array_sha256']
    with np.load(path) as a:x=a['log_width'];logs=a['event_log_likelihoods'];z=a['redshift']
    logs=logs-logs.max(axis=1)[:,None];saved=case['fit'];profile=[]
    for sigma in p['sigma_grid']:
        objective=lambda q:-population_log_likelihood(x,logs,z,[q[0],q[1],np.log(sigma)])
        trials=[minimize(objective,q,method='L-BFGS-B',bounds=[np.log([5,100]),[-2,3]],options={'maxiter':500,'ftol':1e-12,'gtol':1e-6}) for q in [[saved['a'],saved['b']],[np.log(30),0],[np.log(30),1]]]
        valid=[v for v in trials if v.success];assert valid
        fit=min(valid,key=lambda v:v.fun)
        difference=2*(fit.fun-saved['negative_centered_log_likelihood']);assert difference>=-2e-6
        profile.append({'sigma':sigma,'a':float(fit.x[0]),'b':float(fit.x[1]),'profile_lr_relative_to_saved_best':float(max(0,difference)),'successful_starts':len(valid)})
    assert abs(profile[0]['profile_lr_relative_to_saved_best'])<2e-6
    result={'label':case['label'],'truth_b':case['truth_b'],'saved_b':saved['b'],'saved_sigma':saved['sigma'],'optimizer_success':saved['optimizer_success'],'boundary_is_scatter_floor':bool(abs(saved['sigma']-.03)<1e-8 and saved['interior_solution']==False),'profile':profile,'nominal_b95_would_contain_truth_without_validity_screen':bool(case['profile_lr_at_truth']<=3.841458820694124),'injected_sigma_profile_lr':next(v['profile_lr_relative_to_saved_best'] for v in profile if v['sigma']==.1)}
    results.append(result)
out={'scope':p['scope'],'results':results,'source_hashes':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [HERE/'run.py',HERE/'protocol.json',HERE.parent/'timing-population/estimator.py']}}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(results))

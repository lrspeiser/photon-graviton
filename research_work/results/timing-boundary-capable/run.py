"""Development refit with physical sigma=0 boundary; frozen baseline unchanged."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr
H=Path(__file__).resolve().parent;R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-continuous-scatter'))
from integration import averaged_likelihood

def outside(x,mu,sigma):
    if sigma==0:return ((mu<x[0])|(mu>x[-1])).astype(float)
    return ndtr((x[0]-mu)/sigma)+ndtr((mu-x[-1])/sigma)

def fit(x,logs,z,old):
    L=np.exp(logs-logs.max(axis=1)[:,None]);lz=np.log1p(z)
    def objective(q):
        try:v=averaged_likelihood(x,L,q[0]+q[1]*lz,q[2])
        except ArithmeticError:return 1e12
        if np.any(v<=0):return 1e12
        return float(-np.log(v).sum())
    def constraint(q):return 1e-6-outside(x,q[0]+q[1]*lz,q[2])
    bounds=[np.log([5,100]),[-2,3],[0,.6]]
    starts=[[old['a'],old['b'],s] for s in [0,.03,.1]]+[[np.log(30),0,.1],[np.log(30),1,.1]]
    trials=[]
    for q in starts:
        f=minimize(objective,q,method='SLSQP',bounds=bounds,constraints={'type':'ineq','fun':constraint},options={'ftol':1e-10,'maxiter':600})
        trials.append(dict(parameters=f.x.tolist(),objective=float(f.fun),success=bool(f.success),message=str(f.message),maximum_outside=float(outside(x,f.x[0]+f.x[1]*lz,f.x[2]).max())))
    good=[v for v in trials if v['success'] and v['maximum_outside']<=1.00001e-6 and v['objective']<1e11]
    if not good:raise RuntimeError('No valid optimum')
    best=min(good,key=lambda v:v['objective'])
    return dict(best=best,starts=trials,successful_starts=len(good),successful_objective_range=float(max(v['objective'] for v in good)-best['objective']))

if __name__=='__main__':
    source=H.parent/'timing-coverage-calibration/results.json'
    data=json.loads(source.read_text(encoding='utf-8'))
    cases=[v for v in data['cases'] if not v['valid_fit']]
    assert len(cases)==16
    rows=[]
    for c in cases:
        p=R/'research_work/generated/timing-coverage-calibration'/(c['label']+'.npz')
        assert hashlib.sha256(p.read_bytes()).hexdigest()==c['array_sha256']
        with np.load(p) as d:result=fit(d['log_width'],d['event_log_likelihoods'],d['redshift'],c['fit'])
        rows.append(dict(label=c['label'],array_sha256=c['array_sha256'],original=c['fit'],revised=result))
        print(json.dumps({'case':c['label'],'best':result['best'],'start_objective_range':result['successful_objective_range']}),flush=True)
    out=dict(scope='Development on all 16 exposed invalid baseline cases; not full recalibration or scientific inference',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),integration_sha256=hashlib.sha256((H.parent/'timing-continuous-scatter/integration.py').read_bytes()).hexdigest(),cases=rows)
    (H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

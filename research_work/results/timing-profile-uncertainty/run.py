"""Nominal profile intervals on exposed synthetic timing controls."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.stats import chi2
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import population_log_likelihood

def main():
    protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    source=ROOT/protocol['input'];controls=json.loads((source/'results.json').read_text(encoding='utf-8'))
    cases=[(r['label'],r['fit'],r['truth_b']) for r in controls['runs']]
    selected=controls['refinement']['selected_case']
    truth=next(r['truth_b'] for r in controls['runs'] if r['label']==selected)
    cases.append((selected+'-refined',controls['refinement']['fit'],truth))
    results=[]
    for label,fit,truth in cases:
        path=source/(label+'.npz');assert hashlib.sha256(path.read_bytes()).hexdigest()==controls['array_hashes'][path.name]
        with np.load(path) as a:x=a['log_width'];logs=a['event_log_likelihoods'];z=a['redshift']
        logs=logs-logs.max(axis=1)[:,None];cache={}
        def evaluate(b):
            b=float(b)
            if b not in cache:
                objective=lambda q:-population_log_likelihood(x,logs,z,[q[0],b,q[1]])
                runs=[minimize(objective,start,method='L-BFGS-B',bounds=[np.log([5,100]),np.log([.03,.6])],options={'maxiter':500,'ftol':1e-12,'gtol':1e-6}) for start in [[fit['a'],np.log(fit['sigma'])],[np.log(30),np.log(.2)]]]
                valid=[r for r in runs if r.success]
                if not valid:raise RuntimeError(f'No successful nuisance fit: {label}, b={b}')
                best=min(valid,key=lambda r:r.fun)
                cache[b]={'b':b,'nll':float(best.fun),'a':float(best.x[0]),'sigma':float(np.exp(best.x[1])),'successful_starts':len(valid)}
            return cache[b]['nll']
        reference=fit['negative_centered_log_likelihood']
        agreement=abs(evaluate(fit['b'])-reference);assert agreement<1e-6
        ratio=lambda b:2*(evaluate(b)-reference)
        grid=np.unique(np.r_[np.linspace(*protocol['b_range'],101),fit['b'],0.,1.])
        values=np.array([ratio(b) for b in grid]);assert values.min()>-2e-6
        intervals=[]
        for level in protocol['interval_levels']:
            threshold=float(chi2.ppf(level,1));roots=[]
            for l,h,vl,vh in zip(grid[:-1],grid[1:],values[:-1],values[1:]):
                if (vl-threshold)*(vh-threshold)<0:roots.append(float(brentq(lambda b:ratio(b)-threshold,l,h,xtol=1e-9)))
            boundaries=[float(grid[0])]+roots+[float(grid[-1])];regions=[]
            for low,high in zip(boundaries[:-1],boundaries[1:]):
                if ratio((low+high)/2)<=threshold:
                    regions.append({'low':low,'high':high,'lower_truncated':bool(low==grid[0]),'upper_truncated':bool(high==grid[-1])})
            errors=[abs(ratio(root)-threshold) for root in roots];assert max(errors,default=0)<1e-5
            intervals.append({'nominal_level':level,'threshold':threshold,'regions':regions,'contains_injected_truth':any(v['low']<=truth<=v['high'] for v in regions),'maximum_endpoint_error':max(errors,default=0)})
        assert min(v['nll'] for v in cache.values())>=reference-1e-6
        result={'label':label,'injected_b':truth,'fitted_b':fit['b'],'intervals':intervals,'profile_lr_at_b0':ratio(0),'profile_lr_at_b1':ratio(1),'mle_nll_agreement':agreement,'profile_evaluations':sorted(cache.values(),key=lambda v:v['b'])}
        results.append(result)
        print(json.dumps({k:v for k,v in result.items() if k!='profile_evaluations'}),flush=True)
    out={'scope':protocol['scope'],'cases':results,'nominal_not_calibrated':True,'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'run.py',source/'results.json',HERE.parent/'timing-population/estimator.py']}}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

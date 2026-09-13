"""Fresh-seed calibration of the frozen boundary-capable timing estimator."""
from pathlib import Path
import hashlib, importlib.util, itertools, json, sys
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
from scipy.stats import chi2, beta

H=Path(__file__).resolve().parent; R=H.parents[2]
sys.path.insert(0,str(H.parent/'timing-population'))
from estimator import event_log_likelihood, fit_population
sys.path.insert(0,str(H.parent/'timing-boundary-capable'))
from run import fit, outside, averaged_likelihood
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):
    q=p.with_suffix('.tmp');q.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n');q.replace(p)

protocol=json.loads((H/'protocol.json').read_text())
p=json.loads((H.parent/'timing-coverage-calibration/protocol.json').read_text())
sources=[H/'protocol.json',H/'run.py',H.parent/'timing-coverage-calibration/protocol.json',H.parent/'timing-population/estimator.py',H.parent/'timing-injection/run_injections.py',H.parent/'timing-boundary-capable/run.py',H.parent/'timing-continuous-scatter/integration.py',R/p['input_audit'],R/p['source_inventory']]
hashes={str(x.relative_to(R)):sha(x) for x in sources}
out=R/'research_work/generated/timing-fresh-calibration';out.mkdir(exist_ok=True)
cp=out/'checkpoint.json';state=json.loads(cp.read_text()) if cp.exists() else dict(hashes=hashes,cases=[])
assert state['hashes']==hashes
spec=importlib.util.spec_from_file_location('injection',H.parent/'timing-injection/run_injections.py');injection=importlib.util.module_from_spec(spec);spec.loader.exec_module(injection)
audit=json.loads((R/p['input_audit']).read_text());paths={Path(v['path']).name:R/v['cache_relative_path'] for v in audit['input_manifest']}
for v in audit['input_manifest']:assert sha(R/v['cache_relative_path'])==v['sha256']
heads={str(v['SNID']).strip():v for v in fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1)}
phot=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
inventory=json.loads((R/p['source_inventory']).read_text())['runs'][0]['events'];assert len(inventory)==98
slots=[]
for item in inventory:
    h=heads[item['id']];obs=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]
    t=np.asarray(obs['MJD'])-float(h['PEAKMJD']);e=np.asarray(obs['FLUXCALERR'],dtype=float)
    keep=(np.char.strip(obs['BAND'].astype(str))==item['band'])&(abs(t)<=80)&np.isfinite(t)&np.isfinite(e)&(e>0)
    t=t[keep];e=e[keep];order=np.argsort(t);slots.append(dict(id=item['id'],z=item['z'],t=t[order],e=e[order]/np.median(e)))
z=np.array([s['z'] for s in slots]);lz=np.log1p(z);c=p['calibration'];x=np.linspace(*np.log(p['width_days_range']),c['quadrature']['width_nodes'])
for seed,family,snr,truth in itertools.product(protocol['seeds'],c['families'],c['nominal_snr'],c['truth_b']):
    label=f'{family}-snr{snr}-b{truth}-seed{seed}'
    prior=next((v for v in state['cases'] if v['label']==label),None)
    if prior:
        assert sha(out/(label+'.npz'))==prior['array_sha256'];continue
    rng=np.random.default_rng(seed);logs=[]
    for s in slots:
        scale=np.exp(rng.normal(0,c['intrinsic_log_scatter']))*(1+s['z'])**truth
        offset=rng.uniform(-5,5);e=s['e']/snr
        flux=injection.injected(s['t'],scale,offset,family)+rng.normal(0,e)
        logs.append(event_log_likelihood(s['t'],flux,e,x,c['quadrature'],p))
    logs=np.array(logs);initial=fit_population(x,logs,z,p);fitted=fit(x,logs,z,initial);best=fitted['best'];a,b,sigma=best['parameters']
    L=np.exp(logs-logs.max(axis=1)[:,None])
    def objective(q):
        try:v=averaged_likelihood(x,L,q[0]+truth*lz,q[1])
        except ArithmeticError:return 1e12
        return float(-np.log(v).sum()) if np.all(v>0) else 1e12
    trials=[]
    for q in [[a,sigma],[a,0],[np.log(30),.1]]:
        f=minimize(objective,q,method='SLSQP',bounds=[np.log([5,100]),[0,.6]],constraints={'type':'ineq','fun':lambda q:1e-6-outside(x,q[0]+truth*lz,q[1])},options={'ftol':1e-10,'maxiter':600})
        trials.append(dict(parameters=f.x.tolist(),objective=float(f.fun),success=bool(f.success),outside=float(outside(x,f.x[0]+truth*lz,f.x[1]).max())))
    good=[v for v in trials if v['success'] and v['outside']<=1.00001e-6 and v['objective']<1e11]
    lr=2*(min(v['objective'] for v in good)-best['objective']) if good else None
    valid=bool(good) and lr>=-2e-6 and min(a-np.log(5),np.log(100)-a,b+2,3-b,.6-sigma)>1e-5
    array=out/(label+'.npz');np.savez_compressed(array,log_width=x,event_log_likelihoods=logs,redshift=z)
    row=dict(label=label,seed=seed,family=family,snr=snr,truth_b=truth,fit=fitted,truth_profiles=trials,raw_lr=lr,valid=bool(valid),accepted_95=bool(valid and max(0,lr)<=chi2.ppf(.95,1)),array_sha256=sha(array))
    state['cases'].append(row);save(cp,state)
    print(json.dumps(dict(completed=len(state['cases']),total=160,label=label,b=b,valid=bool(valid),lr=lr)),flush=True)
assert len(state['cases'])==160
summary=[]
for family,snr,truth in itertools.product(c['families'],c['nominal_snr'],c['truth_b']):
    rows=[v for v in state['cases'] if (v['family'],v['snr'],v['truth_b'])==(family,snr,truth)]
    n=len(rows);k=sum(v['accepted_95'] for v in rows);bias=float(np.mean([v['fit']['best']['parameters'][1]-truth for v in rows]));invalid=sum(not v['valid'] for v in rows)
    summary.append(dict(family=family,snr=snr,truth_b=truth,n=n,accepted=k,invalid=invalid,mean_bias=bias,binomial_95=[float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.],screen_pass=bool(abs(bias)<=.1 and k/n>=.85 and invalid==0)))
state['summary']=summary;state['real_flux_used']=False;state['physical_validation']=False
save(H/'results.json',state)

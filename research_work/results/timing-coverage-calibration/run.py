"""Resumable, protocol-hashed synthetic coverage calibration."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,sys,time
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import beta,chi2
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=ROOT/'research_work/generated/timing-coverage-calibration'
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import event_log_likelihood,fit_population,log_width_weights,population_log_likelihood

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n');temporary.replace(path)
def interval(k,n):
    return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]

def main():
    p=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'));c=p['calibration']
    hashes={str(q.relative_to(ROOT)):digest(q) for q in [HERE/'protocol.json',HERE/'run.py',HERE.parent/'timing-population/estimator.py',HERE.parent/'timing-injection/run_injections.py',ROOT/p['input_audit'],ROOT/p['source_inventory']]}
    OUT.mkdir(parents=True,exist_ok=True);checkpoint=OUT/'checkpoint.json'
    state=json.loads(checkpoint.read_text(encoding='utf-8')) if checkpoint.exists() else {'hashes':hashes,'cases':{},'refinements':{}}
    if state['hashes']!=hashes:raise RuntimeError('Refusing to mix results from changed protocol/code/inputs')
    spec=importlib.util.spec_from_file_location('injection',HERE.parent/'timing-injection/run_injections.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    audit=json.loads((ROOT/p['input_audit']).read_text(encoding='utf-8'));paths={Path(v['path']).name:ROOT/v['cache_relative_path'] for v in audit['input_manifest']}
    for v in audit['input_manifest']:assert digest(ROOT/v['cache_relative_path'])==v['sha256']
    hd=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1);pd=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
    heads={str(v['SNID']).strip():v for v in hd};inventory=json.loads((ROOT/p['source_inventory']).read_text(encoding='utf-8'))['runs'][0]['events'];assert len(inventory)==98
    slots=[]
    for item in inventory:
        h=heads[item['id']];obs=pd[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]
        t=np.asarray(obs['MJD'])-float(h['PEAKMJD']);e=np.asarray(obs['FLUXCALERR'],dtype=float)
        keep=(np.char.strip(obs['BAND'].astype(str))==item['band'])&(abs(t)<=80)&np.isfinite(t)&np.isfinite(e)&(e>0)
        t=t[keep];e=e[keep];order=np.argsort(t)
        slots.append({'id':item['id'],'z':item['z'],'t':t[order],'e':e[order]/np.median(e)})
    z=np.array([s['z'] for s in slots]);lz=np.log1p(z);x=np.linspace(*np.log(p['width_days_range']),c['quadrature']['width_nodes'])
    def compute(seed,family,snr,b,settings,label):
        start=time.perf_counter();rng=np.random.default_rng(seed);logs=[];scales=[]
        for slot in slots:
            scale=np.exp(rng.normal(0,c['intrinsic_log_scatter']))*(1+slot['z'])**b;offset=rng.uniform(-5,5);e=slot['e']/snr
            flux=old.injected(slot['t'],scale,offset,family)+rng.normal(0,e)
            logs.append(event_log_likelihood(slot['t'],flux,e,x,settings,p));scales.append(scale)
        logs=np.array(logs);fit=fit_population(x,logs,z,p);centered=logs-logs.max(axis=1)[:,None]
        objective=lambda q:-population_log_likelihood(x,centered,z,[q[0],b,q[1]])
        candidates=[minimize(objective,q,method='L-BFGS-B',bounds=[np.log([5,100]),np.log([.03,.6])],options={'maxiter':500,'ftol':1e-12,'gtol':1e-6}) for q in [[fit['a'],np.log(fit['sigma'])],[np.log(30),np.log(.2)]]]
        valid=[v for v in candidates if v.success];ratio=None
        if valid:
            ratio=2*(min(v.fun for v in valid)-fit['negative_centered_log_likelihood'])
            if ratio < -2e-6:raise RuntimeError('Fixed-truth fit beats unrestricted optimum')
            ratio=float(max(0,ratio))
        np.savez_compressed(OUT/(label+'.npz'),log_width=x,event_log_likelihoods=logs,redshift=z,ids=np.array([s['id'] for s in slots]))
        good=fit['optimizer_success'] and fit['interior_solution'] and bool(valid) and fit['maximum_population_mass_outside_width_bounds']<=p['numerical_gate']['maximum_continuous_population_probability_outside_width_bounds']
        return {'label':label,'seed':seed,'family':family,'snr':snr,'truth_b':b,'fit':fit,'true_sample_slope':float(np.dot(lz-lz.mean(),np.log(scales))/np.sum((lz-lz.mean())**2)),'profile_lr_at_truth':ratio,'valid_fit':bool(good),'coverage':{str(level):bool(good and ratio<=chi2.ppf(level,1)) for level in c['levels']},'elapsed_seconds':time.perf_counter()-start,'array_sha256':digest(OUT/(label+'.npz'))},logs
    cells=list(itertools.product(c['families'],c['nominal_snr'],c['truth_b']))
    for seed in c['seeds']:
        for family,snr,b in cells:
            label=f'{family}-snr{snr}-b{b}-seed{seed}'
            if label in state['cases']:
                assert digest(OUT/(label+'.npz'))==state['cases'][label]['array_sha256'];continue
            result,_=compute(seed,family,snr,b,c['quadrature'],label);state['cases'][label]=result;save(checkpoint,state)
            print(json.dumps({'completed':len(state['cases']),'total':len(cells)*len(c['seeds']),'case':label,'b':result['fit']['b'],'lr_at_truth':result['profile_lr_at_truth'],'valid':result['valid_fit']}),flush=True)
    for family,snr,b in cells:
        members=[v for v in state['cases'].values() if (v['family'],v['snr'],v['truth_b'])==(family,snr,b)]
        worst=max(members,key=lambda v:abs(v['fit']['b']-b));label=worst['label']+'-refined'
        if label in state['refinements']:
            assert digest(OUT/(label+'.npz'))==state['refinements'][label]['result']['array_sha256'];continue
        result,fl=compute(worst['seed'],family,snr,b,{'sobol_power':13,'width_nodes':321,'seed':1402},label)
        with np.load(OUT/(worst['label']+'.npz')) as a:cl=a['event_log_likelihoods']
        cl=cl-cl.max(axis=1)[:,None];fl=fl-fl.max(axis=1)[:,None];fit=result['fit']
        post=fl+log_width_weights(x,z,[fit['a'],fit['b'],np.log(fit['sigma'])]);post=np.exp(post-logsumexp(post,axis=1)[:,None])
        change=float(np.mean(np.sum(post*abs(fl-cl),axis=1)));bd=abs(fit['b']-worst['fit']['b'])
        state['refinements'][label]={'result':result,'b_change':bd,'mean_event_log_change':change,'numerical_pass':bool(bd<=.05 and change<=.1 and result['valid_fit'])};save(checkpoint,state)
        print(json.dumps({'refined':label,'b_change':bd,'log_change':change}),flush=True)
    summary=[]
    for family,snr,b in cells:
        members=[v for v in state['cases'].values() if (v['family'],v['snr'],v['truth_b'])==(family,snr,b)]
        errors=np.array([v['fit']['b']-b for v in members]);n=len(members);cov={}
        for level in c['levels']:
            k=sum(v['coverage'][str(level)] for v in members);cov[str(level)]={'covered':k,'total':n,'fraction':k/n,'exact_binomial_95_interval':interval(k,n)}
        summary.append({'family':family,'snr':snr,'truth_b':b,'mean_bias':float(errors.mean()),'bias_monte_carlo_se':float(errors.std(ddof=1)/np.sqrt(n)),'invalid_fits':sum(not v['valid_fit'] for v in members),'coverage':cov})
    save(HERE/'results.json',{'scope':p['scope'],'hashes':hashes,'summary':summary,'cases':list(state['cases'].values()),'refinements':state['refinements'],'real_flux_used':False,'full_scientific_validation':False})
    print('Calibration batch complete; results.json written.',flush=True)

if __name__=='__main__':main()

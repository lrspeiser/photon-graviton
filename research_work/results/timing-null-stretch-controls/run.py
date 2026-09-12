"""Frozen small null/stretch control set, artificial flux only."""
from pathlib import Path
import hashlib, importlib.util, json, sys, time
import numpy as np
from astropy.io import fits
from scipy.special import logsumexp
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import event_log_likelihood,fit_population,log_width_weights

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
    protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    design=protocol['control_design']
    old=module('old_injections',HERE.parent/'timing-injection/run_injections.py')
    audit=json.loads((ROOT/protocol['input_audit']).read_text(encoding='utf-8'))
    paths={Path(r['path']).name:ROOT/r['cache_relative_path'] for r in audit['input_manifest']}
    for r in audit['input_manifest']:
        assert hashlib.sha256((ROOT/r['cache_relative_path']).read_bytes()).hexdigest()==r['sha256']
    hd=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1)
    pd=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
    heads={str(r['SNID']).strip():r for r in hd}
    inventory=json.loads((ROOT/protocol['source_inventory']).read_text(encoding='utf-8'))['runs'][0]['events']
    assert len(inventory)==98
    slots=[]
    for item in inventory:
        h=heads[item['id']];obs=pd[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]
        t=np.asarray(obs['MJD'])-float(h['PEAKMJD']);err=np.asarray(obs['FLUXCALERR'],dtype=float)
        keep=(np.char.strip(obs['BAND'].astype(str))==item['band'])&(abs(t)<=80)&np.isfinite(t)&np.isfinite(err)&(err>0)
        t=t[keep];err=err[keep];order=np.argsort(t);t=t[order]
        err=err[order]/(design['nominal_snr']*np.median(err))
        slots.append({'id':item['id'],'z':item['z'],'t':t,'err':err})
    z=np.array([s['z'] for s in slots]);lx=np.log1p(z)
    x=np.linspace(*np.log(protocol['width_days_range']),design['quadrature']['width_nodes'])
    def simulate(seed,b):
        rng=np.random.default_rng(seed);events=[];scales=[]
        for slot in slots:
            intrinsic=np.exp(rng.normal(0,design['intrinsic_log_scatter']))
            offset=rng.uniform(-5,5);scale=intrinsic*(1+slot['z'])**b
            flux=old.injected(slot['t'],scale,offset,design['family'])+rng.normal(0,slot['err'])
            events.append(flux);scales.append(scale)
        scales=np.log(scales)
        truth_slope=float(np.dot(lx-lx.mean(),scales)/np.sum((lx-lx.mean())**2))
        return events,truth_slope
    def calculate(label,events,settings):
        started=time.perf_counter();logs=[]
        for i,(slot,flux) in enumerate(zip(slots,events)):
            logs.append(event_log_likelihood(slot['t'],flux,slot['err'],x,settings,protocol))
            if (i+1)%49==0:print(label+f': {i+1}/98',flush=True)
        logs=np.array(logs);fit=fit_population(x,logs,z,protocol)
        np.savez_compressed(HERE/(label+'.npz'),log_width=x,event_log_likelihoods=logs,redshift=z,ids=np.array([s['id'] for s in slots]))
        fit['elapsed_seconds']=time.perf_counter()-started
        return fit,logs
    results=[]
    for b in design['truth_b']:
        for seed in design['seeds']:
            label=f'b{b}-seed{seed}';events,truth_slope=simulate(seed,b)
            fit,logs=calculate(label,events,design['quadrature'])
            result={'label':label,'truth_b':b,'seed':seed,'true_finite_sample_slope':truth_slope,'fit':fit,'bias_relative_to_injected_b':fit['b']-b,'difference_from_true_sample_slope':fit['b']-truth_slope}
            results.append(result)
            (HERE/'partial-results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8',newline='\n')
            print(json.dumps(result),flush=True)
    worst=max(results,key=lambda r:abs(r['bias_relative_to_injected_b']))
    events,_=simulate(worst['seed'],worst['truth_b'])
    fine,fl=calculate(worst['label']+'-refined',events,design['refinement'])
    with np.load(HERE/(worst['label']+'.npz')) as a:cl=a['event_log_likelihoods']
    cl=cl-cl.max(axis=1)[:,None];fl=fl-fl.max(axis=1)[:,None]
    post=fl+log_width_weights(x,z,[fine['a'],fine['b'],np.log(fine['sigma'])])
    post=np.exp(post-logsumexp(post,axis=1)[:,None])
    log_change=float(np.mean(np.sum(post*abs(fl-cl),axis=1)))
    b_change=abs(fine['b']-worst['fit']['b'])
    medians={str(b):float(np.median([r['fit']['b'] for r in results if r['truth_b']==b])) for b in design['truth_b']}
    bias_pass=all(abs(medians[str(b)]-b)<=design['maximum_absolute_median_bias_per_truth'] for b in design['truth_b'])
    separation=medians['1']-medians['0']
    numerical_pass=b_change<=design['maximum_refinement_b_change'] and log_change<=design['maximum_refinement_mean_event_log_change']
    optimizer_pass=all(r['optimizer_success'] and r['interior_solution'] and r['maximum_population_mass_outside_width_bounds']<=protocol['numerical_gate']['maximum_continuous_population_probability_outside_width_bounds'] for r in [v['fit'] for v in results]+[fine])
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'run.py',HERE.parent/'timing-population/estimator.py']}
    out={'scope':protocol['scope'],'runs':results,'median_b_by_truth':medians,'median_separation':separation,'median_bias_gate_pass':bool(bias_pass),'refinement':{'selected_case':worst['label'],'fit':fine,'b_change':b_change,'mean_event_log_change':log_change,'gate_pass':bool(numerical_pass)},'optimizer_gate_pass':bool(optimizer_pass),'control_gate_pass':bool(bias_pass and numerical_pass and optimizer_pass and separation>=design['minimum_median_separation']),'real_flux_used':False,'interval_coverage_established':False,'source_hashes':hashes,'array_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.glob('*.npz')}}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'control_gate_pass':out['control_gate_pass'],'median_b_by_truth':medians,'refinement_b_change':b_change,'refinement_log_change':log_change}),flush=True)

if __name__=='__main__':main()

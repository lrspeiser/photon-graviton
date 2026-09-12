"""Frozen numerical convergence pilot on synthetic DES-cadence light curves."""
from pathlib import Path
import argparse, importlib.util, hashlib, json, os, time
import numpy as np
from astropy.io import fits
from scipy.special import logsumexp
from estimator import batched_log_flux,event_log_likelihood,fit_population,log_width_weights,population_log_likelihood

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'timing-population'


def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    global OUT
    parser=argparse.ArgumentParser()
    parser.add_argument('--protocol',type=Path,default=HERE/'protocol.json')
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    OUT=args.output
    protocol=json.loads(args.protocol.read_text())
    scalar=load_module('scalar_component',ROOT/'research_work/results/timing-likelihood/likelihood.py')
    previous=load_module('old_injection',ROOT/'research_work/results/timing-injection/run_injections.py')
    rng=np.random.default_rng(8123);s=rng.uniform(0,1,(3,4,9));err=rng.uniform(.05,.3,9);y=rng.normal(.5,.2,9)
    batch=batched_log_flux(y,err,s)
    scalar_error=max(abs(batch[i,j]-scalar.log_flux_likelihood(y,err,s[i,j])) for i in range(3) for j in range(4))
    assert scalar_error<1e-8
    x=np.linspace(np.log(2),np.log(256),161);z=np.array([.02,.3,.8]);constant=np.array([2.,-1.,.5])
    neutral=max(abs(population_log_likelihood(x,np.repeat(constant[:,None],len(x),axis=1),z,[a,b,np.log(sig)])-constant.sum())
                for a,b,sig in [(np.log(30),0,.1),(np.log(5),-2,.6),(np.log(100),3,.03)])
    assert neutral<1e-10
    try:batched_log_flux([1],[0],[[1]])
    except ValueError:pass
    else:raise AssertionError('Invalid errors were accepted')
    audit=json.loads((ROOT/protocol['input_audit']).read_text())
    paths={Path(r['path']).name:ROOT/r['cache_relative_path'] for r in audit['input_manifest']}
    for item in audit['input_manifest']:
        assert hashlib.sha256((ROOT/item['cache_relative_path']).read_bytes()).hexdigest()==item['sha256']
    hd=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1);pd=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
    heads={str(row['SNID']).strip():row for row in hd}
    old=json.loads((ROOT/protocol['source_inventory']).read_text())
    inventory=old['runs'][0]['events'];assert len(inventory)==98 and not old['feasibility_gate_pass']
    inject=protocol['pilot_injection'];rng=np.random.default_rng(inject['seed']);slots=[]
    for item in inventory:
        h=heads[item['id']];obs=pd[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]
        t=np.asarray(obs['MJD'])-float(h['PEAKMJD']);err=np.asarray(obs['FLUXCALERR'],dtype=float)
        keep=(np.char.strip(obs['BAND'].astype(str))==item['band'])&(abs(t)<=80)&np.isfinite(t)&np.isfinite(err)&(err>0)
        t=t[keep];err=err[keep];order=np.argsort(t);t=t[order];err=err[order]/(inject['nominal_snr']*np.median(err))
        intrinsic=np.exp(rng.normal(0,inject['intrinsic_log_scatter']));offset=rng.uniform(*inject['offset_range_days'])
        scale=intrinsic*(1+item['z'])**inject['b']
        y=previous.injected(t,scale,offset,inject['family'])+rng.normal(0,err)
        slots.append({'id':item['id'],'z':item['z'],'t':t,'err':err,'flux':y,'scale':scale})
    redshift=np.array([s['z'] for s in slots]);all_results={};arrays={};OUT.mkdir(parents=True,exist_ok=True)
    for label,settings in protocol['quadrature'].items():
        start=time.perf_counter();logw=np.linspace(*np.log(protocol['width_days_range']),settings['width_nodes']);event_logs=[]
        for i,slot in enumerate(slots):
            event_logs.append(event_log_likelihood(slot['t'],slot['flux'],slot['err'],logw,settings,protocol))
            if (i+1)%14==0:print(f'{label}: {i+1}/98 events integrated',flush=True)
        event_logs=np.asarray(event_logs);result=fit_population(logw,event_logs,redshift,protocol)
        result['elapsed_seconds']=time.perf_counter()-start;result['event_count']=len(slots)
        all_results[label]=result;arrays[label]=(logw,event_logs)
        np.savez_compressed(OUT/(label+'-event-likelihoods.npz'),log_width=logw,event_log_likelihoods=event_logs,redshift=redshift,ids=np.array([s['id'] for s in slots]))
        print(label+': '+json.dumps(result),flush=True)
    comparisons=[]
    for left,right in [('base','refined'),('refined','independent_scramble')]:
        xl,ll=arrays[left];xr,lr=arrays[right]
        interp=np.array([np.interp(xr,xl,row-row.max()) for row in ll]);centered=lr-lr.max(axis=1)[:,None]
        fit=all_results[right];params=[fit['a'],fit['b'],np.log(fit['sigma'])]
        post=centered+log_width_weights(xr,redshift,params);post=np.exp(post-logsumexp(post,axis=1)[:,None])
        change=float(np.mean(np.sum(post*abs(interp-centered),axis=1)))
        comparisons.append({'left':left,'right':right,'absolute_b_change':abs(all_results[left]['b']-fit['b']),
                            'mean_posterior_weighted_centered_event_log_change':change})
    gate=protocol['numerical_gate']
    passed=(all(c['absolute_b_change']<=gate['maximum_b_change_between_runs'] and
                c['mean_posterior_weighted_centered_event_log_change']<=gate['maximum_mean_absolute_centered_event_log_likelihood_change'] for c in comparisons)
            and all(r['optimizer_success'] and r['interior_solution'] and r['maximum_population_mass_outside_width_bounds']<=gate['maximum_continuous_population_probability_outside_width_bounds'] for r in all_results.values()))
    logz=np.log1p(redshift);true_log_scale=np.log([s['scale'] for s in slots])
    true_slope=float(np.dot(logz-logz.mean(),true_log_scale)/np.sum((logz-logz.mean())**2))
    result={'scope':protocol['scope'],'protocol_sha256':hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
            'checks':{'scalar_batch_max_log_error':float(scalar_error),'constant_likelihood_max_error':float(neutral),'invalid_error_rejected':True,'all_original_cadence_slots':len(slots)},
            'injection':inject,'slope_of_true_injected_widths_in_finite_sample':true_slope,'runs':all_results,'comparisons':comparisons,
            'numerical_pilot_gate_pass':bool(passed),'calibrated_estimator_established':False,'real_flux_used':False,
            'interval_coverage_established':False,'next_step':'Record expanded calibration protocol, separate quadrature/prior sensitivities and test bias/coverage before real data.' if passed else 'Numerical gate failed: refine or replace integration under a new recorded protocol before calibration or real-flux inference.',
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'run.py',HERE/'estimator.py',ROOT/protocol['input_audit'],ROOT/protocol['source_inventory']]},
            'array_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob('*-event-likelihoods.npz'))}}
    (OUT/'timing-population-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'numerical_gate_pass':passed,'comparisons':comparisons,'true_finite_sample_slope':true_slope}),flush=True)


if __name__=='__main__':main()

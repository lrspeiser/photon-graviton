"""Post-failure attribution using known injection truth; not a new inference fit."""
from pathlib import Path
import importlib.util,json,hashlib,os
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize_scalar,brentq

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'timing-diagnosis'
spec=importlib.util.spec_from_file_location('injections',ROOT/'research_work/results/timing-injection/run_injections.py')
injections=importlib.util.module_from_spec(spec);spec.loader.exec_module(injections)

def slope(x,y):return float(np.dot(x-x.mean(),y)/np.sum((x-x.mean())**2))

def main():
    source=ROOT/'research_work/results/timing-injection/timing-injection-results.json'
    prior=json.loads(source.read_text(encoding='utf-8'))
    protocol=json.loads((ROOT/'research_work/results/timing-injection/protocol.json').read_text(encoding='utf-8'))
    assert not prior['feasibility_gate_pass']
    audit=json.loads((ROOT/'research_work/results/timing-foundation/input-audit.json').read_text(encoding='utf-8'))
    paths={Path(r['path']).name:ROOT/r['cache_relative_path'] for r in audit['input_manifest']}
    for r in audit['input_manifest']:
        assert hashlib.sha256((ROOT/r['cache_relative_path']).read_bytes()).hexdigest()==r['sha256']
    hd=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1);pd=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
    heads={str(r['SNID']).strip():r for r in hd};rows=[];replays=0;max_error=0.
    for run in prior['runs']:
        family=run['family'];b=run['injected_propagation_b'];rng=np.random.default_rng(run['seed'])
        shape=lambda u:float(injections.injected(np.array(u),1,0,family))
        peak=minimize_scalar(lambda u:-shape(u),bounds=(-10,20),method='bounded').x
        half=shape(peak)/2
        l=brentq(lambda u:shape(u)-half,-200,peak);r=brentq(lambda u:shape(u)-half,peak,200)
        truth=[];fitted=[];xs=[];accepted=[];geometry=[];status=[]
        for j,event in enumerate(run['events']):
            header=heads[event['id']];obs=pd[int(header['PTROBS_MIN'])-1:int(header['PTROBS_MAX'])]
            t=np.asarray(obs['MJD'])-float(header['PEAKMJD']);err=np.asarray(obs['FLUXCALERR'],dtype=float)
            mask=(np.char.strip(obs['BAND'].astype(str))==event['band'])&(abs(t)<=80)&np.isfinite(t)&np.isfinite(err)&(err>0)
            t=t[mask];err=err[mask];order=np.argsort(t);t=t[order];err=err[order];err=err/(20*np.median(err))
            intrinsic=np.exp(rng.normal(0,.1));offset=rng.uniform(-5,5);scale=intrinsic*(1+event['z'])**b
            flux=injections.injected(t,scale,offset,family)+rng.normal(0,err)
            if j<2:
                width,st=injections.width_fit(t,flux,err,protocol)
                assert st==event['status']
                if width is not None:assert abs(width/event['width_days']-1)<1e-10
                replays+=1
            truth.append((r-l)*scale);xs.append(np.log1p(event['z']))
            fitted.append(event['width_days'] if event['width_days'] is not None else np.nan)
            accepted.append(event['status']=='accepted');status.append(event['status'])
            geometry.append(offset+l*scale>=min(t) and offset+r*scale<=max(t))
        x=np.array(xs);truth=np.array(truth);fitted=np.array(fitted);ok=np.array(accepted);covered=np.array(geometry)
        all_true=slope(x,np.log(truth));accepted_true=slope(x[ok],np.log(truth[ok]));measured=slope(x[ok],np.log(fitted[ok]))
        assert abs(measured-run['recovered_b'])<1e-10
        sample=all_true-b;selection=accepted_true-all_true;estimator=measured-accepted_true
        closure=abs((measured-b)-(sample+selection+estimator));max_error=max(max_error,closure)
        residual=np.log(fitted[ok]/truth[ok]);weight=(x[ok]-x[ok].mean())/np.sum((x[ok]-x[ok].mean())**2)
        contribution=weight*residual
        assert abs(sum(contribution)-estimator)<1e-10
        indices=np.where(ok)[0];top=np.argsort(abs(contribution))[-3:][::-1]
        rows.append({'family':family,'injected_b':b,'seed':run['seed'],
                     'slope_of_all_true_widths':all_true,'slope_of_true_widths_in_accepted_sample':accepted_true,
                     'slope_of_fitted_accepted_widths':measured,'finite_sample_intrinsic_component':sample,
                     'acceptance_selection_component':selection,'fitted_width_error_component':estimator,
                     'true_crossings_within_cadence_range':int(sum(covered)),
                     'accepted_despite_true_crossing_outside_range':int(sum(ok&~covered)),
                     'rejected_despite_true_crossings_inside_range':int(sum(~ok&covered)),
                     'accepted_median_fitted_over_true_width':float(np.median(np.exp(residual))),
                     'accepted_fraction_with_over_20_percent_width_error':float(np.mean(abs(np.exp(residual)-1)>.2)),
                     'largest_absolute_fitting_slope_contributions':[{'id':run['events'][indices[k]]['id'],
                         'z':run['events'][indices[k]]['z'],'fitted_over_true_width':float(np.exp(residual[k])),
                         'slope_contribution':float(contribution[k])} for k in top]})
    result={'scope':'Post-failure diagnosis of the already exposed synthetic timing trials. No new gate, estimator, real-flux inference or validation claim.',
            'source_result_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'runs':rows,
            'checks':{'replayed_individual_fits':replays,'slope_decomposition_max_absolute_error':max_error},
            'identity':'b_hat-b_truth = (b_true_all-b_truth)+(b_true_accepted-b_true_all)+(b_fitted_accepted-b_true_accepted)',
            'interpretation_limit':'True injected widths are unavailable for real events. Geometric endpoint coverage alone is not sufficient sampling or recoverability. Three seeds do not establish confidence coverage; top-contributor identification does not authorize dropping those events.',
            'next_step':'A revised protocol must address fitted-width instability as well as incomplete coverage, and test retained/rejected samples against known truth before interpreting real-data slopes.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'timing-diagnosis-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'checks':result['checks'],'rows':[{k:v for k,v in row.items() if k!='largest_absolute_fitting_slope_contributions'} for row in rows]},indent=2))

if __name__=='__main__':main()

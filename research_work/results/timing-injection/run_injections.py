"""Cadence feasibility only: real observation times/errors, artificial fluxes."""
from pathlib import Path
import hashlib,json,os
import numpy as np
from astropy.io import fits
from scipy.optimize import least_squares,brentq

ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'timing-injection'

def pulse(t,logamp,peak,logrise,logdelta,baseline):
    rise=np.exp(logrise);fall=rise+np.exp(logdelta)
    shift=rise*np.log(fall/rise-1)
    u=t-peak+shift
    logshape=-u/fall-np.logaddexp(0,-u/rise)
    logmax=-shift/fall-np.logaddexp(0,-shift/rise)
    return baseline+np.exp(logamp+logshape-logmax)

def width_fit(t,f,err,protocol):
    b=protocol['fit_bounds']
    low=[b['log_amplitude'][0],b['peak_days'][0],np.log(b['rise_days'][0]),np.log(b['fall_minus_rise_days'][0]),b['baseline'][0]]
    high=[b['log_amplitude'][1],b['peak_days'][1],np.log(b['rise_days'][1]),np.log(b['fall_minus_rise_days'][1]),b['baseline'][1]]
    trials=[]
    for rise in [3,7,15]:
        p=[0,float(np.clip(t[np.argmax(f)],-20,20)),np.log(rise),np.log(20),0]
        fit=least_squares(lambda p:(pulse(t,*p)-f)/err,p,bounds=(low,high),max_nfev=b['max_nfev'],ftol=1e-7,xtol=1e-7,gtol=1e-7)
        if fit.success:trials.append(fit)
    if not trials:return None,'optimizer'
    fit=min(trials,key=lambda x:np.dot(x.fun,x.fun));p=fit.x
    if np.any(p-np.array(low)<1e-4) or np.any(np.array(high)-p<1e-4):return None,'parameter_boundary'
    peak=p[1];fhalf=lambda time:(pulse(np.array(time),*p)-p[4])/np.exp(p[0])-.5
    try:
        left=brentq(fhalf,peak-2000,peak);right=brentq(fhalf,peak,peak+2000)
    except ValueError:return None,'crossing'
    if left<min(t) or right>max(t):return None,'coverage'
    return float(right-left),'accepted'

def injected(t,scale,offset,family):
    u=(t-offset)/scale
    sigma=np.where(u<0,8.,18.)
    f=np.exp(-.5*(u/sigma)**2)
    if family=='split_gaussian_shoulder':f+=.2*np.exp(-.5*((u-20)/10)**2)
    return f

def main():
    protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    audit=json.loads((ROOT/protocol['input_audit']).read_text(encoding='utf-8'))
    paths={Path(r['path']).name:ROOT/r['cache_relative_path'] for r in audit['input_manifest']}
    for r in audit['input_manifest']:
        assert hashlib.sha256((ROOT/r['cache_relative_path']).read_bytes()).hexdigest()==r['sha256']
    hd=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1);pd=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
    slots=[]
    for row in hd:
        z=float(row['REDSHIFT_HELIO']);peak=float(row['PEAKMJD'])
        if row['FAKE']!=0 or row['SNTYPE']!=1 or not np.isfinite(z+peak) or z<=0:continue
        choices=sorted((abs(center/(1+z)/450-1),band) for band,center in protocol['filter_centers_nm'].items())
        mismatch,band=choices[0]
        if mismatch>.05:continue
        obs=pd[int(row['PTROBS_MIN'])-1:int(row['PTROBS_MAX'])]
        t=np.asarray(obs['MJD'])-peak;err=np.asarray(obs['FLUXCALERR'],dtype=float)
        mask=(np.char.strip(obs['BAND'].astype(str))==band)&(abs(t)<=80)&np.isfinite(t)&np.isfinite(err)&(err>0)
        t=t[mask];err=err[mask]
        if len(t)<8 or sum(t<0)<2 or sum(t>0)<3:continue
        order=np.argsort(t);t=t[order];err=err[order]
        # Never read the real flux column: artificial flux has declared SNR.
        err=err/(20*np.median(err))
        slots.append({'id':str(row['SNID']).strip(),'z':z,'band':band,'t':t,'err':err})
    print('Eligible cadence slots:',len(slots),flush=True)
    assert len(slots)>2,'Too few slots even for a diagnostic regression'
    rows=[]
    for family in protocol['injections']['families']:
        for b in protocol['injections']['propagation_exponents']:
            for seed in protocol['injections']['seeds']:
                rng=np.random.default_rng(seed);widths=[];zs=[];rejections={};event_rows=[]
                for s in slots:
                    intrinsic=np.exp(rng.normal(0,.1));offset=rng.uniform(-5,5)
                    scale=intrinsic*(1+s['z'])**b
                    noiseless=injected(s['t'],scale,offset,family)
                    # Source evolution e=b and propagation b=0 is exactly identical.
                    equivalent=injected(s['t'],intrinsic*(1+s['z'])**(0+b),offset,family)
                    assert np.array_equal(noiseless,equivalent)
                    f=noiseless+rng.normal(0,s['err'])
                    width,status=width_fit(s['t'],f,s['err'],protocol)
                    event_rows.append({'id':s['id'],'z':s['z'],'band':s['band'],'width_days':width,'status':status})
                    rejections[status]=rejections.get(status,0)+1
                    if width is not None:widths.append(width);zs.append(s['z'])
                if len(widths)>=3:
                    design=np.column_stack([np.ones(len(zs)),np.log1p(zs)])
                    beta=np.linalg.lstsq(design,np.log(widths),rcond=None)[0]
                    bhat=float(beta[1])
                else:bhat=None
                row={'family':family,'injected_propagation_b':b,'seed':seed,'accepted':len(widths),'eligible':len(slots),
                     'accepted_fraction':len(widths)/len(slots),'rejection_counts':rejections,'recovered_b':bhat,'events':event_rows}
                rows.append(row);print(f'{family} b={b} seed={seed}: {len(widths)}/{len(slots)}, recovered={bhat}',flush=True)
    gates=protocol['feasibility_gate'];groups=[]
    for family in protocol['injections']['families']:
        for b in protocol['injections']['propagation_exponents']:
            group=[r for r in rows if r['family']==family and r['injected_propagation_b']==b]
            vals=[r['recovered_b'] for r in group if r['recovered_b'] is not None]
            median=float(np.median(vals)) if vals else None
            groups.append({'family':family,'injected_b':b,'median_recovered_b':median,
                           'bias_gate_pass':len(vals)==len(group) and abs(median-b)<=gates['maximum_absolute_median_exponent_bias_each_family_truth']})
    passed=len(slots)>=gates['minimum_eligible_events'] and all(r['accepted_fraction']>=gates['minimum_accepted_fraction_each_run'] for r in rows) and all(g['bias_gate_pass'] for g in groups)
    result={'scope':protocol['scope'],'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
            'input_audit_sha256':hashlib.sha256((ROOT/protocol['input_audit']).read_bytes()).hexdigest(),
            'eligible_events':len(slots),'eligible_redshift_range':[min(s['z'] for s in slots),max(s['z'] for s in slots)],
            'runs':rows,'group_summaries':groups,'feasibility_gate_pass':passed,
            'real_flux_timing_inference_performed':False,'source_evolution_degeneracy':'Artificial curves depend on b_propagation+e_intrinsic; matching sums give exactly identical noiseless data.',
            'next_step':'Limited cadence gate passed; broader wavelength/SNR/selection/coverage checks remain before real-data inference.' if passed else 'Feasibility gate failed; diagnose bias or acceptance before any real-flux exponent fit. Do not relax frozen gates retrospectively.'}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'timing-injection-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'gate_pass':passed,'groups':groups},indent=2),flush=True)

if __name__=='__main__':main()

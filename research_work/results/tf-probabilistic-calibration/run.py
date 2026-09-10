from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.optimize import minimize,brentq
from scipy.special import logsumexp,lambertw,ndtr
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
calpath=ROOT/'research_work/results/indicator-tf-calibration/cross-validation.csv'
qpath=ROOT/'research_work/results/tf-quality-audit/joined.csv'
rows=[r for r in csv.DictReader(calpath.open()) if r['photometry_p']=='1']
q={int(r['pgc']):r for r in csv.DictReader(qpath.open())}
x=np.array([float(q[int(r['pgc'])]['logWmxi'])-2.5 for r in rows]);sx=np.array([float(q[int(r['pgc'])]['e_logWmxi']) for r in rows])
y=np.array([float(r['observed_absolute_magnitude']) for r in rows]);g=np.array([int(r['group']) for r in rows])
vy=np.array([float(r['sigma_mu'])**2+.05**2+float(q[int(r['pgc'])]['e_Ai'])**2 for r in rows])
assert len(rows)==73 and np.all(np.isfinite(vy))
def components(theta,idx,mix):
    a,b,ls=theta[:3];s=np.exp(ls);f=theta[3] if mix else 0.
    v=vy[idx]+a*a*sx[idx]**2
    return y[idx]-a*x[idx]-b,v+s*s,v+25*s*s,f
def logdensity(r,v1,v2,f):
    l1=-.5*(np.log(2*np.pi*v1)+r*r/v1)
    if f==0:return l1
    l2=-.5*(np.log(2*np.pi*v2)+r*r/v2)
    return logsumexp(np.array([np.log1p(-f)+l1,np.log(f)+l2]),axis=0)
fits=[]
def fit(idx,mix,label):
    bounds=[(-20,-1),(-30,-10),(np.log(.01),np.log(5))]+([(0,.5)] if mix else [])
    candidates=[]
    for s in [.2,.6,1.2]:
        init=[-7.5,-20.8,np.log(s)]+([.15] if mix else [])
        result=minimize(lambda t:-float(np.sum(logdensity(*components(t,idx,mix)))),init,method='L-BFGS-B',bounds=bounds,options=dict(maxiter=1000,ftol=1e-11,gtol=1e-6))
        candidates.append(result)
    valid=[r for r in candidates if r.success and np.isfinite(r.fun)]
    assert valid, (label,[r.message for r in candidates])
    best=min(valid,key=lambda r:r.fun)
    fits.append(dict(model='mixture' if mix else 'gaussian',fold=label,successful_starts=len(valid),start_objectives=[float(r.fun) for r in candidates],best_parameters=best.x.tolist(),
        boundary=any(abs(v-lo)<1e-5 or abs(v-hi)<1e-5 for v,(lo,hi) in zip(best.x,bounds))))
    return best.x
out=[];summary=[];cdf_errors=[]
for mix in [False,True]:
    label='mixture' if mix else 'gaussian';full=fit(np.ones(len(rows),bool),mix,'full')
    selected=[];seen=np.zeros(len(rows),int)
    for group in np.unique(g):
        held=g==group;theta=fit(~held,mix,str(group));res,v1,v2,f=components(theta,held,mix)
        for j,k in enumerate(np.where(held)[0]):
            seen[k]+=1
            def cdf(v):return (1-f)*ndtr(v/np.sqrt(v1[j]))+f*ndtr(v/np.sqrt(v2[j]))
            h68=brentq(lambda z:cdf(z)-.84,0,12*np.sqrt(v2[j]));h95=brentq(lambda z:cdf(z)-.975,0,12*np.sqrt(v2[j]))
            cdf_errors.extend([abs(cdf(h68)-.84),abs(cdf(h95)-.975)])
            mean=theta[0]*x[k]+theta[1];app=10**((float(rows[k]['m'])-mean-25)/5);alpha=.0002488993286382367;d=float(lambertw(alpha*app).real/alpha)
            logp=float(logdensity(res[j],v1[j],v2[j],f))
            responsibility=0. if f==0 else float(np.exp(np.log(f)-.5*(np.log(2*np.pi*v2[j])+res[j]**2/v2[j])-logp))
            item=dict(model=label,pgc=rows[k]['pgc'],name=rows[k]['name'],group=int(group),prediction=mean,residual=float(res[j]),log_density=logp,
                residual_halfwidth68=h68,residual_halfwidth95=h95,covered68=bool(abs(res[j])<=h68),covered95=bool(abs(res[j])<=h95),
                distance_ratio=d/float(rows[k]['distance']),broad_responsibility=responsibility)
            selected.append(item);out.append(item)
    assert np.all(seen==1)
    summary.append(dict(model=label,slope=float(full[0]),intercept=float(full[1]),intrinsic_core_mag=float(np.exp(full[2])),broad_fraction=float(full[3]) if mix else 0.,
        heldout_log_score=sum(r['log_density'] for r in selected),coverage68=sum(r['covered68'] for r in selected),coverage95=sum(r['covered95'] for r in selected),n=len(selected),
        median_halfwidth68=float(np.median([r['residual_halfwidth68'] for r in selected])),rms_mag=float(np.sqrt(np.mean([r['residual']**2 for r in selected]))),
        median_abs_fractional_distance_discrepancy=float(np.median([abs(r['distance_ratio']-1) for r in selected]))))
assert max(cdf_errors)<1e-10 and sum(r['pgc']=='40809' for r in out)==2
with (HERE/'predictions.csv').open('w',encoding='utf8',newline='') as file:
    w=csv.DictWriter(file,fieldnames=list(out[0]),lineterminator='\n');w.writeheader();w.writerows(out)
result=dict(n=73,groups=len(np.unique(g)),models=summary,max_cdf_error=max(cdf_errors),optimizer_fits=len(fits),fits_with_any_failed_start=sum(r['successful_starts']<3 for r in fits),boundary_fits=sum(r['boundary'] for r in fits),
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [calpath,qpath]},status='Conditional plug-in residual predictions, not population selection model or new holdout')
paired={label:{r['pgc']:r for r in out if r['model']==label} for label in ['gaussian','mixture']}
gain={pgc:paired['mixture'][pgc]['log_density']-r['log_density'] for pgc,r in paired['gaussian'].items()}
result['posthoc_score_concentration_diagnostic']=dict(total_gain=sum(gain.values()),ngc4424_gain=gain['40809'],remaining_72_gain=sum(v for k,v in gain.items() if k!='40809'),interpretation='Accounting diagnostic only; no rows removed or refitted')
for name,value in [('results.json',result),('optimizer.json',fits)]:
    (HERE/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(result,indent=2))

"""Fixed exact-third capture with training-only conservative redistribution."""
import json,io,zipfile,hashlib,sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
P=Path(__file__).resolve().parent;ROOT=P.parents[2];BASE=ROOT/'temporal_candidate_audit/data'
FOLLOWUP='--compact-followup' in sys.argv
base=json.loads((P/'third-radiation-retention-results.json').read_text());cp=base['models']['attenuated'];assert cp['q']==1/3
for fn,h in base['input_sha256'].items():assert hashlib.sha256((BASE/fn).read_bytes()).hexdigest()==h
splits=json.loads((BASE/'sparc_frozen.json').read_text())['split'];props={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)!=19:continue
    try:props[f[0]]=(float(f[11]),float(f[7]))
    except ValueError:pass
data=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as z:
    for split,names in splits.items():
        for name in names:
            t=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))));vb=t[:,3]*abs(t[:,3])+.5*t[:,4]*abs(t[:,4])+.7*t[:,5]*abs(t[:,5]);ok=np.isfinite(t).all(axis=1)&(t[:,0]>0)&(t[:,1]>0)&(t[:,2]>0)&(vb>0);t,vb=t[ok],vb[ok]
            rd,L=props[name];X=L/rd**2;eta=X**(1/3)/(1+X**(1/3))
            data.append(dict(name=name,split=split,r=t[:,0],v=t[:,1],vb=vb,rd=rd,eta=eta))
assert len(data)==149 and sum(len(d['v']) for d in data)==3150
def prepare(n,order):
    mu,w=leggauss(order);tables=[]
    for d in data:
        r=np.geomspace(d['rd']*1e-8,max(d['r'][-1]*(16.01 if FOLLOWUP else 4.01),d['rd']*100),n);a=cp['scale_to_disk']*d['rd'];x=r/a
        t=x[:,None]*mu;B2=1+x[:,None]**2*(1-mu**2);B=np.sqrt(B2)
        tau=cp['k0_per_kpc']*a*(t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3))
        J=.5*np.sum(np.exp(-np.maximum(tau,0))*w,axis=1);rho=2*cp['C_Msun_kpc3']*d['eta']*J/(1+x*x)**2
        m=4*np.pi*(rho[0]*r[0]**3/3+cumulative_trapezoid(rho*r*r,r,initial=0))
        tables.append(PchipInterpolator(np.log(r),m,extrapolate=False))
    return tables
def setting(kind,p,eta):
    if kind=='baseline':return 1.,1.
    if kind=='shared':return np.exp(p[0]),1.
    if kind in ['partial','compact_partial']:return np.exp(p[1]),p[0]
    return np.exp(p[0]*(1-eta)),1.
def predict(kind,p,tables):
    result=[]
    for d,M in zip(data,tables):
        s,f=setting(kind,p,d['eta']);mass=(1-f)*M(np.log(d['r']))+f*M(np.log(d['r']/s))
        result.append(np.sqrt(d['vb']+4.30091727003628e-6*mass/d['r']))
    return result
def scores(pred):
    return {split:dict(n=sum(d['split']==split for d in data),RMSE_kms=float(np.sqrt(np.mean([np.mean((v-d['v'])**2) for d,v in zip(data,pred) if d['split']==split]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(v/d['v'])**2) for d,v in zip(data,pred) if d['split']==split])))) for split in splits}
coarse=prepare(2049,64);fine=prepare(4097,96);lim=float(np.log(4));models={};predictions=[]
for kind in (['baseline','compact_partial'] if FOLLOWUP else ['baseline','shared','partial','retention_conditioned']):
    partial=kind in ['partial','compact_partial']
    bounds=[(0,1),(-2*lim,0 if FOLLOWUP else lim)] if partial else [(-lim,lim)]
    if partial and not FOLLOWUP:bounds[1]=(-lim,lim)
    starts=[[f,b] for f in ([.005,.02,.1] if FOLLOWUP else [.25,.75,1.]) for b in ([-2.5,-1.4] if FOLLOWUP else [-1.,1.])] if partial else [[-1.],[0.],[1.]]
    def loss(p):return scores(predict(kind,p,coarse))['train']['log_RMS']**2
    opts=[] if kind=='baseline' else [minimize(loss,p,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-14,'gtol':1e-9,'maxiter':300}) for p in starts]
    opt=min(opts,key=lambda o:o.fun) if opts else None;p=opt.x if opt is not None else []
    pred=predict(kind,p,fine);sc=scores(pred);cs=scores(predict(kind,p,coarse));drift=max(abs(sc[k]['RMSE_kms']-cs[k]['RMSE_kms']) for k in splits);assert drift<.05
    models[kind]=dict(parameters=list(p),scores=sc,refinement_max_RMS_change_kms=drift,boundary=False if opt is None else any(min(abs(v-lo),abs(v-hi))<1e-5 for v,(lo,hi) in zip(p,bounds)),starts=[dict(parameters=o.x.tolist(),loss=float(o.fun),success=bool(o.success)) for o in opts])
    for d,v in zip(data,pred):predictions.append(dict(model=kind,galaxy=d['name'],split=d['split'],predicted_kms=v.tolist(),observed_kms=d['v'].tolist()))
    print(kind,json.dumps(models[kind]),flush=True)
assert all(abs(models['baseline']['scores'][k]['RMSE_kms']-cp['finer_scores'][k]['RMSE_kms'])<.05 for k in splits)
prefix='redistribution-compact' if FOLLOWUP else 'redistribution'
(P/(prefix+'-results.json')).write_text(json.dumps(dict(models=models,capture_input_sha256=hashlib.sha256((P/'third-radiation-retention-results.json').read_bytes()).hexdigest(),input_sha256=base['input_sha256']),indent=2)+'\n')
(P/(prefix+'-predictions.json')).write_text(json.dumps(predictions,indent=2)+'\n')

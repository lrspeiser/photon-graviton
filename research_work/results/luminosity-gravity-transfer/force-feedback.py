"""Evolve positive deposition under force feedback, fit and transfer."""
import hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
BASE=ROOT/'temporal_candidate_audit/data'; G=4.30091727003628e-6
previous=json.loads((HERE/'results.json').read_text())
for fn,digest in previous['input_sha256'].items():
    assert hashlib.sha256((BASE/fn).read_bytes()).hexdigest()==digest
splits=json.loads((BASE/'sparc_frozen.json').read_text())['split']
assert list(splits)==['train','validation','test'] and len(splits['train'])==89
names=sum(splits.values(),[]); data=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as z:
    for name in names:
        a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))))
        vb=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5])
        good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb>0)
        a,vb=a[good],vb[good]
        data.append(dict(name=name,R=a[:,0],v=a[:,1],vb=vb))
assert len(data)==149 and sum(len(d['v']) for d in data)==3150

class Calculation:
    def __init__(self,n,steps):
        self.x=np.linspace(0,1,n+1); self.steps=steps
        self.R=np.array([d['R'][-1] for d in data])[:,None]
        self.gb=np.array([np.interp(self.x*d['R'][-1],np.r_[0,d['R']],np.r_[0,d['vb']/d['R']]) for d in data])
        self.x2=self.x**2; self.denom=np.where(self.x2==0,1,self.x2)
    def force(self,D,C):
        enclosed=cumulative_trapezoid(D*self.x2,self.x,axis=1,initial=0)
        return 4*np.pi*G*C*self.R*enclosed/self.denom
    def calculate(self,p,feedback):
        C,g0=10**p[0],10**p[1]
        if feedback:
            D=np.zeros_like(self.gb); dt=1/self.steps
            for _ in range(self.steps):
                g=self.gb+self.force(D,C); k=(g/(g+g0))**2
                gm=self.gb+self.force(D+.5*dt*k,C)
                D+=dt*(gm/(gm+g0))**2
        else: D=(self.gb/(self.gb+g0))**2
        assert np.min(D)>=0 and np.max(D)<=1+1e-12
        extra=self.force(D,C)
        return [np.sqrt(d['vb']+d['R']*np.interp(d['R']/d['R'][-1],self.x,extra[i])) for i,d in enumerate(data)]

def score(pred):
    out={}; k=0
    for split,nn in splits.items():
        vv=pred[k:k+len(nn)]; dd=data[k:k+len(nn)]; k+=len(nn)
        out[split]=dict(n=len(nn),RMSE_kms=float(np.sqrt(np.mean([np.mean((v-d['v'])**2) for v,d in zip(vv,dd)]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(v/d['v'])**2) for v,d in zip(vv,dd)]))))
    return out

coarse=Calculation(256,32); fine=Calculation(512,64)
results=dict(models={},input_sha256=previous['input_sha256']); rows=[]
bounds=[(-2,12),(-2,7)]
for feedback in [False,True]:
    name='feedback' if feedback else 'fixed_baryonic_force'
    def loss(p):
        vv=coarse.calculate(p,feedback)
        return np.mean([np.mean(np.log10(v/d['v'])**2) for v,d in zip(vv[:89],data[:89])])
    opts=[]
    for start in [[7,3],[5,2],[9,4]]:
        opt=minimize(loss,start,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-12,'gtol':1e-7,'maxiter':150})
        opts.append(opt); print(name,'start',start,'result',opt.x,'loss',opt.fun,flush=True)
    opt=min(opts,key=lambda o:o.fun); p=opt.x
    pred=coarse.calculate(p,feedback); finer=fine.calculate(p,feedback)
    scores=score(pred); fs=score(finer)
    drift=max(abs(scores[s]['RMSE_kms']-fs[s]['RMSE_kms']) for s in splits)
    results['models'][name]=dict(C_Msun_kpc3=float(10**p[0]),g0_kms2_per_kpc=float(10**p[1]),scores=scores,finer_scores=fs,refinement_max_RMS_change_kms=drift,refinement_pass=drift<.1,optimizer_success=bool(opt.success),boundary=bool(any(min(abs(v-lo),abs(v-hi))<1e-5 for v,(lo,hi) in zip(p,bounds))),starts=[dict(p=o.x.tolist(),loss=float(o.fun),success=bool(o.success)) for o in opts])
    for d,v in zip(data,finer):
        split=next(s for s in splits if d['name'] in splits[s])
        rows.append(dict(model=name,galaxy=d['name'],split=split,R_kpc=d['R'].tolist(),observed_kms=d['v'].tolist(),predicted_kms=v.tolist()))
for fn,obj in [('force-feedback-results.json',results),('force-feedback-predictions.json',rows)]:
    (HERE/fn).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(results['models'],indent=2),flush=True)

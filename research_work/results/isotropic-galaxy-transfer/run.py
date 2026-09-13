"""Fit attenuated all-angle external capture, then transfer frozen rule."""
import hashlib, io, json, zipfile
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
BASE=ROOT/'temporal_candidate_audit/data'; G=4.30091727003628e-6
previous=json.loads((HERE.parent/'luminosity-gravity-transfer/results.json').read_text())
for fn,digest in previous['input_sha256'].items():
    assert hashlib.sha256((BASE/fn).read_bytes()).hexdigest()==digest
splits=json.loads((BASE/'sparc_frozen.json').read_text())['split']
assert list(splits)==['train','validation','test'] and len(splits['train'])==89
rd={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)!=19: continue
    try: rd[f[0]]=float(f[11])
    except ValueError: pass
data=[]
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as z:
    for name in sum(splits.values(),[]):
        a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))))
        vb=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5])
        good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb>0)
        a,vb=a[good],vb[good]
        data.append(dict(name=name,R=a[:,0],v=a[:,1],vb=vb,rd=rd[name]))
assert len(data)==149 and sum(len(d['v']) for d in data)==3150

class Calculator:
    def __init__(self,n,order):
        self.x=np.linspace(0,1,n+1)
        self.r=np.array([d['R'][-1]*self.x for d in data])
        self.mu,self.w=leggauss(order)
        self.rd=np.array([d['rd'] for d in data])[:,None]
    def predict(self,p,attenuated):
        C=10**p[0]; scale=10**p[-1]; a=scale*self.rd
        y=self.r/a; shape=(1+y*y)**-2
        if attenuated:
            k0=10**p[1]
            t=y[:,:,None]*self.mu
            B2=1+y[:,:,None]**2*(1-self.mu**2); B=np.sqrt(B2)
            # Integral from upstream infinity to the point; reversing mu leaves J unchanged.
            primitive=t/(2*B2*(B2+t*t))+(np.arctan(t/B)+np.pi/2)/(2*B**3)
            tau=np.maximum(0,k0*a[:,:,None]*primitive)
            J=.5*np.sum(np.exp(-tau)*self.w,axis=2)
        else: J=np.ones_like(shape)
        mass=4*np.pi*C*cumulative_trapezoid(shape*J*self.r**2,self.r,axis=1,initial=0)
        pred=[]
        for i,d in enumerate(data):
            M=np.interp(d['R'],self.r[i],mass[i])
            pred.append(np.sqrt(d['vb']+G*M/d['R']))
        return pred

def scores(pred):
    out={}; j=0
    for s,nn in splits.items():
        pairs=list(zip(pred[j:j+len(nn)],data[j:j+len(nn)]));j+=len(nn)
        out[s]=dict(n=len(nn),RMSE_kms=float(np.sqrt(np.mean([np.mean((v-d['v'])**2) for v,d in pairs]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(v/d['v'])**2) for v,d in pairs]))))
    return out

coarse=Calculator(256,32); fine=Calculator(512,64)
results=dict(models={},input_sha256=previous['input_sha256']); rows=[]
for attenuated in [False,True]:
    name='attenuated' if attenuated else 'transparent_control'
    bounds=[(-2,12),(-6,2),(-1,2)] if attenuated else [(-2,12),(-1,2)]
    starts=[[7,-2,0],[7,-1,.5],[6,-4,1]] if attenuated else [[7,0],[7,.5],[6,1]]
    def loss(p):
        pred=coarse.predict(p,attenuated)
        return np.mean([np.mean(np.log10(v/d['v'])**2) for v,d in zip(pred[:89],data[:89])])
    opts=[]
    for start in starts:
        opt=minimize(loss,start,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-12,'gtol':1e-7,'maxiter':150})
        opts.append(opt);print(name,start,opt.x,opt.fun,flush=True)
    opt=min(opts,key=lambda o:o.fun);p=opt.x
    pred=fine.predict(p,attenuated);cs=scores(coarse.predict(p,attenuated));fs=scores(pred)
    drift=max(abs(cs[s]['RMSE_kms']-fs[s]['RMSE_kms']) for s in splits)
    results['models'][name]=dict(C_Msun_kpc3=float(10**p[0]),scale_to_disk=float(10**p[-1]),k0_per_kpc=float(10**p[1]) if attenuated else None,scores=cs,finer_scores=fs,refinement_max_RMS_change_kms=drift,refinement_pass=drift<.1,optimizer_success=bool(opt.success),boundary=bool(any(min(abs(v-lo),abs(v-hi))<1e-5 for v,(lo,hi) in zip(p,bounds))),starts=[dict(parameters=o.x.tolist(),loss=float(o.fun),success=bool(o.success)) for o in opts])
    for d,v in zip(data,pred):
        rows.append(dict(model=name,galaxy=d['name'],split=next(s for s in splits if d['name'] in splits[s]),R_kpc=d['R'].tolist(),observed_kms=d['v'].tolist(),predicted_kms=v.tolist()))
for fn,obj in [('results.json',results),('predictions.json',rows)]:
    (HERE/fn).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(results['models'],indent=2),flush=True)

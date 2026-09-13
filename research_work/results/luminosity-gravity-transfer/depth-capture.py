"""Shared depth-dependent deposition: actual galaxy fit and frozen transfer."""
import hashlib
import io
import json
import zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]; BASE=ROOT/'temporal_candidate_audit/data'
G=4.30091727003628e-6 # kpc (km/s)^2 / Msun
previous=json.loads((HERE/'results.json').read_text())
for fn,digest in previous['input_sha256'].items():
    assert hashlib.sha256((BASE/fn).read_bytes()).hexdigest()==digest
splits=json.loads((BASE/'sparc_frozen.json').read_text())['split']
gal={}
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as z:
    for name in sum(splits.values(),[]):
        a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))))
        vb=a[:,3]*abs(a[:,3])+.5*a[:,4]*abs(a[:,4])+.7*a[:,5]*abs(a[:,5])
        good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb>0)
        a,vb=a[good],vb[good]; R=a[:,0]
        assert np.all(np.diff(R)>0)
        r=np.unique(np.r_[0,np.geomspace(R[0]*1e-5,R[-1],1200),R])
        gb=np.interp(r,np.r_[0,R],np.r_[0,vb/R])
        integ=cumulative_trapezoid(gb,r,initial=0)
        W=vb[-1]+integ[-1]-integ
        assert np.all(W>0)
        gal[name]=dict(R=R,v=a[:,1],vb=vb,r=r,W=W,index=np.searchsorted(r,R))
assert len(gal)==149 and sum(len(d['v']) for d in gal.values())==3150

def predict(d,p):
    C,W0=10**p[0],10**p[1]
    rho=C*(d['W']/(d['W']+W0))**4
    mass=4*np.pi*cumulative_trapezoid(rho*d['r']**2,d['r'],initial=0)[d['index']]
    extra=G*mass/d['R']
    return np.sqrt(d['vb']+extra),extra/d['vb']
def loss(p):
    return np.mean([np.mean(np.log10(predict(gal[n],p)[0]/gal[n]['v'])**2) for n in splits['train']])
bounds=[(-2,12),(1,7)]
opts=[minimize(loss,p,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-14,'gtol':1e-9}) for p in [[7,4],[5,3],[9,5]]]
opt=min(opts,key=lambda o:o.fun); p=opt.x
out=dict(C_Msun_kpc3=float(10**p[0]),W0_kms2=float(10**p[1]),optimizer_success=bool(opt.success),boundary=bool(any(min(abs(v-a),abs(v-b))<1e-5 for v,(a,b) in zip(p,bounds))),scores={},galaxies=[],input_sha256=previous['input_sha256'],starts=[dict(parameters=o.x.tolist(),loss=float(o.fun),success=bool(o.success)) for o in opts])
rows=[]; ratios=[]
for split,names in splits.items():
    mse=[]; logs=[]
    for name in names:
        d=gal[name]; yp,ratio=predict(d,p); ratios.extend(ratio.tolist())
        mse.append(float(np.mean((yp-d['v'])**2))); logs.append(float(np.mean(np.log10(yp/d['v'])**2)))
        out['galaxies'].append(dict(name=name,split=split,RMSE_kms=mse[-1]**.5,log_RMS=logs[-1]**.5))
        rows.append(dict(name=name,split=split,R_kpc=d['R'].tolist(),observed_kms=d['v'].tolist(),predicted_kms=yp.tolist(),extra_to_baryonic_force=ratio.tolist()))
    out['scores'][split]=dict(n=len(names),RMSE_kms=float(np.sqrt(np.mean(mse))),log_RMS=float(np.sqrt(np.mean(logs))))
out['feedback_diagnostic']=dict(points_extra_force_exceeds_baryonic=int(np.sum(np.array(ratios)>1)),points=len(ratios),median_extra_to_baryonic=float(np.median(ratios)))
transition=[float(d['vb'][-1]*d['R'][-1]/10**p[1]) for d in gal.values()]
out['Kepler_W_equals_W0_radius_kpc']=dict(min=min(transition),median=float(np.median(transition)),max=max(transition),warning='Conditional continuation, not observed halo extent')
for fn,obj in [('depth-capture-results.json',out),('depth-capture-predictions.json',rows)]:
    (HERE/fn).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k not in ['galaxies','input_sha256']},indent=2))

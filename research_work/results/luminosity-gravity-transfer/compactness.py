"""Fit fixed luminosity/area proxy, preserving historical galaxy splits."""
import io
import hashlib
import json
import zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
BASE=ROOT/'temporal_candidate_audit/data'
previous=json.loads((HERE/'results.json').read_text())
for fn,digest in previous['input_sha256'].items():
    assert hashlib.sha256((BASE/fn).read_bytes()).hexdigest()==digest
astar=previous['a_star_m_s2']; KPC=3.085677581491367e19
splits=json.loads((BASE/'sparc_frozen.json').read_text())['split']
meta={}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)!=19: continue
    try: meta[f[0]]=(float(f[7]),float(f[11]))
    except ValueError: pass
gal={}
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as z:
    for name in sum(splits.values(),[]):
        x=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+'_rotmod.dat'))))
        vb=x[:,3]*abs(x[:,3])+.5*x[:,4]*abs(x[:,4])+.7*x[:,5]*abs(x[:,5])
        good=np.isfinite(x).all(axis=1)&(x[:,0]>0)&(x[:,1]>0)&(x[:,2]>0)&(vb>0)
        x,vb=x[good],vb[good]
        L,rd=meta[name]
        assert L>0 and rd>0 and len(x)>=5
        gal[name]=dict(R=x[:,0],v=x[:,1],vb=vb,gb=vb*1e6/(x[:,0]*KPC),X=(L/10)/(rd/3)**2)
assert len(gal)==149 and sum(len(d['v']) for d in gal.values())==3150
def predict(d,p):
    ge=10**p[0]*astar*(d['gb']/astar)**p[1]*d['X']**p[2]
    return np.sqrt(d['vb']+ge*d['R']*KPC/1e6)
def loss(p):
    return np.mean([np.mean(np.log10(predict(gal[n],p)/gal[n]['v'])**2) for n in splits['train']])
bounds=[(-4,2),(0,1),(-1,1)]
opts=[minimize(loss,p,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-14,'gtol':1e-9}) for p in [[-.6,.5,0],[-1,0,0],[-1,1,0]]]
opt=min(opts,key=lambda o:o.fun); p=opt.x
out=dict(A=float(10**p[0]),p=float(p[1]),q=float(p[2]),a_star_m_s2=astar,optimizer_success=bool(opt.success),boundary=bool(any(min(abs(v-a),abs(v-b))<1e-5 for v,(a,b) in zip(p,bounds))),scores={},galaxies=[],input_sha256=previous['input_sha256'])
predictions=[]
for split,names in splits.items():
    mse=[]; logs=[]
    for n in names:
        d=gal[n]; yp=predict(d,p)
        mse.append(float(np.mean((yp-d['v'])**2)))
        logs.append(float(np.mean(np.log10(yp/d['v'])**2)))
        out['galaxies'].append(dict(name=n,split=split,X=d['X'],RMSE_kms=mse[-1]**.5,log_RMS=logs[-1]**.5))
        predictions.append(dict(name=n,split=split,R_kpc=d['R'].tolist(),observed_kms=d['v'].tolist(),predicted_kms=yp.tolist()))
    out['scores'][split]=dict(n=len(names),RMSE_kms=float(np.sqrt(np.mean(mse))),log_RMS=float(np.sqrt(np.mean(logs))))
# Independently reproduce the previous q=0 control on identical input selection.
base=previous['models']['baseline']; bp=[np.log10(base['A']),base['p'],0]
for s,names in splits.items():
    score=float(np.sqrt(np.mean([np.mean((predict(gal[n],bp)-gal[n]['v'])**2) for n in names])))
    assert abs(score-base['scores'][s]['RMSE_kms'])<1e-8
for fn,obj in [('compactness-results.json',out),('compactness-predictions.json',predictions)]:
    (HERE/fn).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in out.items() if k not in ['galaxies','input_sha256']},indent=2))

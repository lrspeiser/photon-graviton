"""Actual rotation fit/transfer; choices recorded in protocol.md."""
import hashlib
import io
import json
import zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = ROOT / 'temporal_candidate_audit/data'
KPC = 3.085677581491367e19
old = json.loads((HERE.parent/'joint-galaxy-audit/results.json').read_text())
astar = old['sparc']['parameters']['a_star_m_s2']
splits = json.loads((BASE/'sparc_frozen.json').read_text())['split']
meta = {}
for line in (BASE/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f = line.split()
    if len(f) != 19:
        continue
    try:
        meta[f[0]] = dict(inc=float(f[5]), L=float(f[7]), rd=float(f[11]), Q=int(f[17]))
    except ValueError:
        pass
gal = {}
with zipfile.ZipFile(BASE/'Rotmod_LTG.zip') as z:
    for fn in sorted(z.namelist()):
        name = fn.replace('_rotmod.dat', '')
        m = meta[name]
        if m['Q'] > 2 or m['inc'] < 30 or m['rd'] <= 0:
            continue
        x = np.atleast_2d(np.loadtxt(io.BytesIO(z.read(fn))))
        vb2 = x[:,3]*abs(x[:,3]) + .5*x[:,4]*abs(x[:,4]) + .7*x[:,5]*abs(x[:,5])
        good = np.isfinite(x).all(axis=1)&(x[:,0]>0)&(x[:,1]>0)&(x[:,2]>0)&(vb2>0)
        if good.sum() < 5:
            continue
        x, vb2 = x[good], vb2[good]
        assert m['L'] > 0
        gal[name] = dict(R=x[:,0], v=x[:,1], vb2=vb2, gb=vb2*1e6/(x[:,0]*KPC), Lratio=m['L']/10)
assert len(gal)==149 and sum(len(d['v']) for d in gal.values())==3150
assert set(gal)==set(sum(splits.values(), []))

def velocity(d, p):
    q = p[2] if len(p)==3 else 0.
    ge = 10**p[0]*astar*(d['gb']/astar)**p[1]*d['Lratio']**q
    return np.sqrt(d['vb2']+ge*d['R']*KPC/1e6)

def loss(p):
    return np.mean([np.mean(np.log10(velocity(gal[n],p)/gal[n]['v'])**2) for n in splits['train']])

results = dict(a_star_m_s2=astar, models={}, input_sha256={})
predictions = []
for model in ['baseline', 'luminosity']:
    bounds = [(-4,2),(0,1)] + ([(-1,1)] if model=='luminosity' else [])
    starts = [[-.6,.5],[-1,0],[-1,1]]
    if model=='luminosity':
        starts = [p+[0.] for p in starts]
    opts = [minimize(loss,p,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-14,'gtol':1e-9}) for p in starts]
    opt = min(opts,key=lambda o:o.fun)
    p = opt.x
    out = dict(A=float(10**p[0]),p=float(p[1]),q=float(p[2]) if len(p)==3 else 0.,
               optimizer_success=bool(opt.success), boundary=bool(any(min(abs(v-lo),abs(v-hi))<1e-5 for v,(lo,hi) in zip(p,bounds))), scores={},galaxies=[])
    for split,names in splits.items():
        mse, logs = [], []
        for n in names:
            d = gal[n]; v = velocity(d,p)
            mse.append(float(np.mean((v-d['v'])**2)))
            logs.append(float(np.mean(np.log10(v/d['v'])**2)))
            out['galaxies'].append(dict(name=n,split=split,RMSE_kms=mse[-1]**.5,log_RMS=logs[-1]**.5,L_3p6_1e9_Lsun=d['Lratio']*10))
            for r,obs,yp in zip(d['R'],d['v'],v):
                predictions.append(dict(model=model,galaxy=n,split=split,R_kpc=float(r),observed_kms=float(obs),predicted_kms=float(yp)))
        out['scores'][split] = dict(n=len(names),RMSE_kms=float(np.sqrt(np.mean(mse))),log_RMS=float(np.sqrt(np.mean(logs))))
    results['models'][model] = out
for fn in ['SPARC_Lelli2016c.mrt','Rotmod_LTG.zip','sparc_frozen.json']:
    results['input_sha256'][fn] = hashlib.sha256((BASE/fn).read_bytes()).hexdigest()
for split in splits:
    assert abs(results['models']['baseline']['scores'][split]['RMSE_kms']-old['sparc']['scores']['power'][split]['galaxy_weighted_RMSE_kms']) < .001
for fn,obj in [('results.json',results),('predictions.json',predictions)]:
    (HERE/fn).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:{a:b for a,b in v.items() if a!='galaxies'} for k,v in results['models'].items()},indent=2))

from pathlib import Path
import json,hashlib,importlib.util,sys
import numpy as np
import pandas as pd
from scipy.special import roots_legendre,eval_legendre
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
outer=float(sys.argv[1]);assert outer in [30.,200.]
OUT=HERE/f'outer{outer:g}';OUT.mkdir(exist_ok=True)
reference=json.loads((HERE.parent/'coupled-capture-calibration/results.json').read_text())
src=HERE.parent/'conservative-field-completion/run.py'
spec=importlib.util.spec_from_file_location('gravity',src);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
prior_path=HERE.parent/'milky-way-depth-capture/results.json';prior=json.loads(prior_path.read_text())
star_path=ROOT/'research_work/data-cache/cepheid-stars/cepheid-common-frame-train-selected.parquet'
stars=pd.read_parquet(star_path,columns=['role','R_kpc','bin']);assert stars.role.eq('train').all() and len(stars)==542
rs=stars.R_kpc.to_numpy();bins=stars.bin.to_numpy();observed=np.array(prior['radial_observed']);base=np.array(prior['ordinary_radial_prediction'])**2
print('Building fixed refined ordinary field',flush=True);baryons=g.Baryons(True)

class Model:
    def __init__(self,nr,na,lmax):
        self.nr=nr;self.na=na;self.lmax=lmax
        self.r=r=np.geomspace(1e-4,outer,nr);mu,w=roots_legendre(na);self.ells=ells=np.arange(0,lmax+1,2)
        self.leg=eval_legendre(ells[:,None],mu[None,:]).T;self.angular=self.leg*w[:,None]*(2*ells+1)[None,:]/2
        self.wb=np.empty((nr,na))
        for start in range(0,nr,8):
            rr=np.broadcast_to(r[start:start+8,None],(len(r[start:start+8]),na));mm=np.broadcast_to(mu,rr.shape)
            self.wb[start:start+len(rr)]=-baryons.evaluate(rr.ravel(),mm.ravel())[0].reshape(rr.shape)
        self.dr=np.r_[np.diff(r)[0]/2,(r[2:]-r[:-2])/2,np.diff(r)[-1]/2]
        big=np.maximum(r[:,None],r[None,:]);ratio=np.minimum(r[:,None],r[None,:])/big
        self.kernel=np.stack([-4*np.pi*g.G/(2*l+1)*ratio**l/big*(r*r*self.dr)[None,:] for l in ells])
    def potential(self,rho):return np.einsum('lij,jl->il',self.kernel,rho@self.angular,optimize=False)
    def run(self,p,C,steps):
        rho=np.zeros_like(self.wb);h=1/steps
        def rate(rho):
            W=np.maximum(self.wb-self.potential(rho)@self.leg.T,0.)
            return C*(W/(W+40000))**p
        for i in range(steps):
            midpoint=rho+h/2*rate(rho);rho+=h*rate(midpoint)
            assert rho.min()>=0 and rho.max()<=C*(i+1)/steps*(1+1e-12)
        spline=CubicSpline(np.log(self.r),self.potential(rho),axis=0);basis=eval_legendre(self.ells,0.)
        extra=np.sum(spline(np.log(rs),1)*basis,axis=1);assert np.all(extra>0)
        predicted=np.sqrt(base+np.array([extra[bins==i].mean() for i in range(12)]))
        return dict(C=C,p=p,radial_nodes=self.nr,angular_nodes=self.na,lmax=self.lmax,time_steps=steps,
            rms=float(np.sqrt(np.mean((predicted-observed)**2))),predictions=predicted.tolist(),residuals=(predicted-observed).tolist(),
            mass_Msun=float(4*np.pi*np.sum((rho@self.angular)[:,0]*self.r**2*self.dr)),
            solar_deposit_depth=float(-np.sum(spline(np.log(8.2))*basis)),density_over_C_range=[float(rho.min()/C),float(rho.max()/C)])

coarse=Model(400,64,24);fits=[]
for p in [2,3]:
    old_C=next(f['best']['C'] for f in reference['fits'] if f['p']==p)
    transferred=coarse.run(p,old_C,128)
    print(json.dumps(dict(stage='transferred',outer=outer,p=p,result=transferred)),flush=True)
    cache={}
    def trial(logC):
        key=float(logC)
        if key not in cache:cache[key]=coarse.run(p,10**key,128)
        return cache[key]
    grid=np.linspace(4,8.5,13);scan=[]
    for x in grid:
        r=trial(x);scan.append(r);print(json.dumps(dict(stage='scan',p=p,C=r['C'],rms=r['rms'])),flush=True)
    candidates=list(scan);optimizers=[]
    for i in range(1,len(grid)-1):
        if scan[i]['rms']<=scan[i-1]['rms'] and scan[i]['rms']<=scan[i+1]['rms']:
            fit=minimize_scalar(lambda x:trial(x)['rms']**2,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-5})
            assert fit.success;record=trial(fit.x);candidates.append(record);optimizers.append(dict(logC=float(fit.x),success=bool(fit.success),evaluations=int(fit.nfev)))
    best=min(candidates,key=lambda r:r['rms']);edge=bool(best['C']==scan[0]['C'] or best['C']==scan[-1]['C'])
    sensitivity=[coarse.run(p,best['C']*factor,128) for factor in [.99,1.01]]
    timecheck=coarse.run(p,best['C'],256)
    fits.append(dict(p=p,transferred=transferred,scan=scan,optimizers=optimizers,edge_optimum=edge,best=best,sensitivity=sensitivity,timecheck=timecheck))
    print(json.dumps(dict(stage='fit',p=p,best=best)),flush=True)
    (OUT/'calibration.json').write_text(json.dumps(fits,indent=2)+'\n',encoding='utf8',newline='\n')
del coarse
fine=Model(800,128,48);checks=[]
for fit in fits:
    best=fit['best'];refined=fine.run(fit['p'],best['C'],256);fit['refined']=refined
    dt=float(max(abs(np.array(best['predictions'])-fit['timecheck']['predictions'])))
    ds=float(max(abs(np.array(best['predictions'])-refined['predictions'])))
    checks.append(dict(p=fit['p'],time_prediction_difference=dt,combined_prediction_difference=ds,passes=dt<.5 and ds<.5))
    print(json.dumps(dict(stage='refined',p=fit['p'],refined=refined,check=checks[-1])),flush=True)
out=dict(outer_kpc=outer,fits=fits,checks=checks,observed_training_proxy=observed.tolist(),ordinary_prediction=np.sqrt(base).tolist(),
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [src,prior_path,star_path]},
    scope='Boundary sensitivity of coupled growth, including transferred and refitted amplitudes; no new holdout or physical boundary law.')
(OUT/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')

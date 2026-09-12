from pathlib import Path
import json,importlib.util,hashlib,sys
import numpy as np
import pandas as pd
from scipy.special import roots_legendre,eval_legendre
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
path=HERE.parent/'conservative-field-completion/run.py'
spec=importlib.util.spec_from_file_location('gravity',path);g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
prior_path=HERE.parent/'milky-way-depth-capture/results.json';prior=json.loads(prior_path.read_text())
stars=pd.read_parquet(ROOT/'research_work/data-cache/cepheid-stars/cepheid-common-frame-train-selected.parquet',columns=['role','R_kpc','bin']);assert stars.role.eq('train').all()
rs=stars.R_kpc.to_numpy();bins=stars.bin.to_numpy();observed=np.array(prior['radial_observed']);base=np.array(prior['ordinary_radial_prediction'])**2
print('Building ordinary field',flush=True);baryons=g.Baryons(True)
refine='--refine' in sys.argv
old=json.loads((HERE/'results.json').read_text()) if refine else {}
runs=[r for r in old.get('runs',[]) if r['time_steps']!=256]
grids=[(400,64,24,[256]),(800,128,48,[256])] if refine else [(400,64,24,[64]),(800,128,48,[64,128])]
for nr,na,lmax,steps_list in grids:
    r=np.geomspace(1e-4,100,nr);mu,w=roots_legendre(na);ells=np.arange(0,lmax+1,2)
    leg=eval_legendre(ells[:,None],mu[None,:]).T;angular=leg*w[:,None]*(2*ells+1)[None,:]/2
    wb=np.empty((nr,na))
    for start in range(0,nr,8):
        rr=np.broadcast_to(r[start:start+8,None],(len(r[start:start+8]),na));mm=np.broadcast_to(mu,rr.shape)
        wb[start:start+len(rr)]=-baryons.evaluate(rr.ravel(),mm.ravel())[0].reshape(rr.shape)
    dr=np.r_[np.diff(r)[0]/2,(r[2:]-r[:-2])/2,np.diff(r)[-1]/2]
    small=np.minimum(r[:,None],r[None,:]);big=np.maximum(r[:,None],r[None,:]);ratio=small/big
    kernel=np.stack([-4*np.pi*g.G/(2*l+1)*ratio**l/big*(r*r*dr)[None,:] for l in ells])
    def potential(rho):return np.einsum('lij,jl->il',kernel,rho@angular,optimize=False)
    def summary(rho):
        ph=potential(rho);spline=CubicSpline(np.log(r),ph,axis=0);basis=eval_legendre(ells,0.)
        grad=np.sum(spline(np.log(rs),1)*basis,axis=1)/rs
        extra=rs*grad;assert np.all(extra>0)
        extra_bins=np.array([extra[bins==i].mean() for i in range(12)])
        pred=np.sqrt(base+extra_bins)
        mass=float(4*np.pi*np.sum((rho@angular)[:,0]*r*r*dr))
        solar=float(-np.sum(spline(np.log(8.2))*basis))
        return dict(predictions=pred.tolist(),rms=float(np.sqrt(np.mean((pred-observed)**2))),mass_Msun=mass,solar_deposit_depth=solar)
    for p in [2,3]:
        if refine and p!=3:continue
        previous=next(a for a in prior['runs'] if a['refined'] and a['p']==p and a['outer_kpc']==100)
        C=previous['density_scale_Msun_kpc3']
        control=summary(C*(wb/(wb+40000))**p)
        control_error=float(max(abs(np.array(control['predictions'])-previous['radial_predictions'])))
        assert control_error<.5
        def rate(rho):
            total=np.maximum(wb-potential(rho)@leg.T,0.)
            return C*(total/(total+40000))**p
        for steps in steps_list:
            rho=np.zeros_like(wb);h=1/steps;history=[]
            for i in range(steps):
                midpoint=rho+h/2*rate(rho);rho+=h*rate(midpoint)
                assert rho.min()>=0 and rho.max()<=C*(i+1)/steps*(1+1e-12)
                if (i+1)%(steps//8)==0:
                    snap=summary(rho);history.append(dict(exposure=(i+1)/steps,mass_Msun=snap['mass_Msun'],solar_deposit_depth=snap['solar_deposit_depth']))
            final=summary(rho)
            record=dict(p=p,radial_nodes=nr,angular_nodes=na,lmax=lmax,time_steps=steps,C=C,control=control,control_max_prediction_difference=control_error,final=final,
                density_over_C_range=[float(rho.min()/C),float(rho.max()/C)],history=history)
            runs.append(record);print(json.dumps({k:v for k,v in record.items() if k not in ['control','history']}),flush=True)
checks=[r for r in old.get('checks',[]) if r.get('stage')!='refined'] if refine else []
for p in [2,3]:
    if refine and p!=3:continue
    coarse=next(r for r in runs if r['p']==p and r['radial_nodes']==400 and r['time_steps']==(256 if refine else 64))
    fine=next(r for r in runs if r['p']==p and r['radial_nodes']==800 and r['time_steps']==(256 if refine else 128))
    timecontrol=next(r for r in runs if r['p']==p and r['radial_nodes']==800 and r['time_steps']==(128 if refine else 64))
    dv=float(max(abs(np.array(coarse['final']['predictions'])-fine['final']['predictions'])))
    dm=abs(coarse['final']['mass_Msun']/fine['final']['mass_Msun']-1)
    dt=float(max(abs(np.array(timecontrol['final']['predictions'])-fine['final']['predictions'])))
    checks.append(dict(p=p,stage='refined' if refine else 'initial',spacetime_prediction_difference=dv,mass_fraction_difference=dm,time_prediction_difference=dt,passes=dv<2 and dm<.01 and dt<2))
(HERE/'results.json').write_text(json.dumps(dict(runs=runs,checks=checks,prior_sha256=hashlib.sha256(prior_path.read_bytes()).hexdigest(),scope='Frozen-amplitude chronological capture feedback; no new observational fit, support or energy supply proof.'),indent=2)+'\n',encoding='utf8',newline='\n')

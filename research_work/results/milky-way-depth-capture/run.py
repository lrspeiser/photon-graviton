from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
import pandas as pd
from scipy.special import roots_legendre
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
source=HERE.parent/'conservative-field-completion/run.py'
spec=importlib.util.spec_from_file_location('poisson',source);f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
train_path=ROOT/'research_work/data-cache/cepheid-stars/cepheid-common-frame-train-selected.parquet'
bins_path=HERE.parent/'cepheid-common-frame/training-bins.json';vertical_path=HERE.parent/'milky-way-capture/inputs.json'
stars=pd.read_parquet(train_path,columns=['role','R_kpc','bin']);assert stars.role.eq('train').all() and len(stars)==542
bins=json.loads(bins_path.read_text());assert len(bins)==12 and all(r['moment_valid'] for r in bins)
vertical=json.loads(vertical_path.read_text())['bovy']['rows'];assert len(vertical)==43
R=stars.R_kpc.to_numpy();membership=stars.bin.to_numpy();obs=np.array([b['jeans_proxy_kms'] for b in bins])
vr=np.array([v['R_kpc'] for v in vertical]);vobs=np.array([v['Kz_over_2piG_Msun_pc2'] for v in vertical])
print('Building refined ordinary matter',flush=True);baryons=f.Baryons(True)
b2=-R*f.force(baryons,R,0)[:,0]
base=np.array([np.sqrt(b2[membership==i].mean()) for i in range(12)])
vbase=f.force(baryons,vr,1.1)
runs=[]
for refined in [False,True]:
    for outer in [100.,200.]:
        r=np.geomspace(1e-4,outer,1200 if refined else 600);ells=np.arange(0,49 if refined else 25,2)
        mu,w=roots_legendre(192 if refined else 96);leg,_=f.basis(ells,mu);angular=leg*(w[:,None]*(2*ells+1)[None,:]/2)
        shape=np.empty((len(r),len(mu)))
        for start in range(0,len(r),8):
            rr=np.broadcast_to(r[start:start+8,None],(len(r[start:start+8]),len(mu)));mm=np.broadcast_to(mu,rr.shape)
            W=np.maximum(-baryons.evaluate(rr.ravel(),mm.ravel())[0].reshape(rr.shape),0.)
            shape[start:start+len(rr)]=W/(W+40000)
        for p in [2,3]:
            coeff=(1e6*shape**p)@angular
            inn,out=f.scaled_integrals(r,coeff*r[:,None],coeff*r[:,None],ells);fac=4*np.pi*f.G/(2*ells+1)
            field=f.Expansion(r,ells,-fac*(inn+out),fac*((ells+1)*inn-ells*out)/r[:,None])
            extra2=-R*f.force(field,R,0)[:,0];assert np.all(extra2>0)
            extra_bin=np.array([extra2[membership==i].mean() for i in range(12)])
            base_bin=base**2
            objective=lambda scale:float(np.mean((np.sqrt(base_bin+scale*extra_bin)-obs)**2))
            fit=minimize_scalar(objective,bounds=(0,1000),method='bounded',options={'xatol':1e-9});assert fit.success
            scale=float(fit.x);assert 1e-6<scale<999.999
            pred=np.sqrt(base_bin+scale*extra_bin)
            totalvertical=vbase+scale*f.force(field,vr,1.1)
            vp=abs(totalvertical[:,1])/(2*np.pi*f.G*1e6)
            phi_b=float(baryons.evaluate(np.array([8.2]),np.array([0.]))[0][0]);phi_d=float(field.evaluate(np.array([8.2]),np.array([0.]))[0][0])*scale
            grid=[]
            for rad in [2.,5.,8.2,15.]:
                for z in [.3,1.,3.]:
                    bforce=f.force(baryons,rad,z)[0];dforce=scale*f.force(field,rad,z)[0]
                    grid.append(dict(R_kpc=rad,z_kpc=z,ordinary_aR=float(bforce[0]),ordinary_az=float(bforce[1]),deposit_aR=float(dforce[0]),deposit_az=float(dforce[1]),units='(km/s)^2/kpc'))
            record=dict(refined=refined,outer_kpc=outer,p=p,density_scale_Msun_kpc3=scale*1e6,
                radial_training_rms=float(np.sqrt(objective(scale))),vertical_development_rms=float(np.sqrt(np.mean((vp-vobs)**2))),vertical_bias=float(np.mean(vp-vobs)),
                deposit_mass_Msun=float(4*np.pi*inn[-1,0]*r[-1]*scale),solar_deposit_to_baryon_depth=float(phi_d/phi_b),
                radial_predictions=pred.tolist(),vertical_predictions=vp.tolist(),force_grid=grid)
            runs.append(record);print(json.dumps({k:v for k,v in record.items() if not isinstance(v,list)}),flush=True)
ref=[]
for outer in [100.,200.]:
    for p in [2,3]:
        lo=next(r for r in runs if not r['refined'] and r['outer_kpc']==outer and r['p']==p);hi=next(r for r in runs if r['refined'] and r['outer_kpc']==outer and r['p']==p)
        dr=float(max(abs(np.array(lo['radial_predictions'])-hi['radial_predictions'])));dv=float(max(abs(np.array(lo['vertical_predictions'])-hi['vertical_predictions'])))
        ref.append(dict(outer=outer,p=p,radial_max_difference=dr,vertical_max_difference=dv,passes=dr<.5 and dv<1.))
vbp=abs(vbase[:,1])/(2*np.pi*f.G*1e6)
out=dict(training_stars=len(stars),training_bins=12,vertical_exposed_rows=43,
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,train_path,bins_path,vertical_path]},
    ordinary_radial_rms=float(np.sqrt(np.mean((base-obs)**2))),ordinary_vertical_rms=float(np.sqrt(np.mean((vbp-vobs)**2))),
    radial_observed=obs.tolist(),ordinary_radial_prediction=base.tolist(),vertical_observed=vobs.tolist(),ordinary_vertical_prediction=vbp.tolist(),
    runs=runs,refinement=ref,vertical_admissibility='Excluded from accepted observational validation: published Kz inferences use a potential family including a dark halo. Retained only as a historical alternate-model diagnostic; not fitted.',scope='Thin-capture one-amplitude training calibration; exposed vertical diagnostic; no new holdout or self-gravity opacity feedback.')
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')

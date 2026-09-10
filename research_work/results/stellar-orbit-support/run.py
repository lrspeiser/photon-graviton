"""Use real training observations to diagnose the orbit library, not rank gravity."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from scipy.special import logsumexp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
ORBIT=HERE.parent/'rotating-bar-orbits'
sys.path.insert(0,str(ORBIT))
# Load under a distinct module name to avoid importing this script as run.
import importlib.util
spec=importlib.util.spec_from_file_location('orbit_driver',ORBIT/'run.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
Field=module.Field;DiskGrid=module.DiskGrid
CACHE=ROOT/'research_work/data-cache/stellar-catalogs'
OUT=ROOT/'research_work/data-cache/stellar-orbit-support';OUT.mkdir(parents=True,exist_ok=True)
P=json.loads((HERE/'protocol.json').read_text());OMEGA=P['bar_pattern_speed_kms_per_kpc']
LABELS=['R_kpc','z_kpc','phi_rad','vR_kms','vphi_kms','vz_kms']

def save(name,value):
    (HERE/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def cartesian(y):
    R,z,phi,vr,vp,vz=y.T;c=np.cos(phi);s=np.sin(phi)
    return np.c_[R*c,R*s,z],np.c_[vr*c-vp*s,vr*s+vp*c,vz]

def invariants(field,y):
    x,p=cartesian(y);potential=[];acc=[]
    for start in range(0,len(y),128):
        a,b=field.evaluate(x[start:start+128]);potential.extend(a);acc.extend(b)
    potential=np.array(potential);acc=np.array(acc)
    R,z,phi,vr,vp,vz=y.T;c=np.cos(phi);s=np.sin(phi)
    ar=acc[:,0]*c+acc[:,1]*s;ap=-acc[:,0]*s+acc[:,1]*c
    J=.5*np.sum(p*p,axis=1)+potential-OMEGA*R*vp
    gradient=np.c_[-ar-OMEGA*vp,-acc[:,2],-R*ap,vr,vp-OMEGA*R,vz]
    return J,gradient

def fit_weights(y,sigma,centres,broadening=0):
    variance=sigma*sigma+broadening*broadening
    logK=-.5*((y[:,None]-centres[None,:])**2/variance[:,None]+np.log(2*np.pi*variance[:,None]))
    weights=np.ones(len(centres))/len(centres);old=-np.inf
    history=[]
    for iteration in range(500):
        logits=logK+np.log(np.maximum(weights,1e-300))
        norm=logsumexp(logits,axis=1)
        objective=float(np.mean(norm));history.append(objective)
        assert objective>=old-1e-10
        response=np.exp(logits-norm[:,None])
        updated=response.mean(axis=0);updated/=updated.sum()
        if iteration and abs(objective-old)<1e-9:
            weights=updated;break
        weights=updated;old=objective
    nll=-float(np.mean(logsumexp(logK+np.log(np.maximum(weights,1e-300)),axis=1)))
    return dict(weights=weights.tolist(),training_mean_negative_log_density=nll,
                iterations=iteration+1,converged=iteration<499,
                fixed_J_broadening_kms2=broadening,
                interpretation='In-sample library diagnostic; not a gravitational-model score')

def support(y,sigma,centres):
    distance=np.min(abs(y[:,None]-centres[None,:]),axis=1)
    return dict(fraction_beyond_3_approximate_sigma=float(np.mean(distance>3*sigma)),
                median_nearest_gap_kms2=float(np.median(distance)),
                p90_nearest_gap_kms2=float(np.quantile(distance,.9)),
                fraction_outside_centre_range=float(np.mean((y<centres.min())|(y>centres.max()))))

def anchors(y,identifiers):
    R,z,phi,vr,vp,vz=y.T
    # Components in physical bar Cartesian space, plus cylindrical velocities.
    x,_=cartesian(y);features=np.c_[x,vr,vp,vz]
    indices=[];strata=[]
    for ir,(lo,hi) in enumerate(zip([.5,3.5,5],[3.5,5,9])):
        for iz,(zl,zh) in enumerate(zip([0,.2,.5],[.2,.5,1.5])):
            mask=(R>=lo)&((R<hi) if hi<9 else (R<=hi))&(abs(z)>=zl)&((abs(z)<zh) if zh<1.5 else (abs(z)<=zh))
            loc=np.where(mask)[0]
            if not len(loc):continue
            f=features[loc];median=np.median(f,axis=0)
            scale=np.maximum(np.quantile(f,.84,axis=0)-np.quantile(f,.16,axis=0),[.2,.2,.1,10,10,10])
            norm=(f-median)/scale
            chosen=[int(np.argmin(np.sum(norm*norm,axis=1)))]
            nearest=np.sum((norm-norm[chosen[0]])**2,axis=1)
            for _ in range(min(8,len(loc))-1):
                j=int(np.argmax(nearest));chosen.append(j)
                nearest=np.minimum(nearest,np.sum((norm-norm[j])**2,axis=1));nearest[chosen]=-np.inf
            indices.extend(loc[chosen]);strata.extend([(ir,iz)]*len(chosen))
    indices=np.array(indices)
    assert len(np.unique(indices))==len(indices)
    return indices,strata

def main():
    src=CACHE/'stellar-errors-conditional_minus_0.017.parquet'
    d=pd.read_parquet(src,filters=[('holdout_role','==','training')])
    assert (d.holdout_role=='training').all()
    chemistry=pd.read_parquet(CACHE/'matched-with-gaia-covariance.parquet',columns=['source_id','FE_H','ALPHA_M'])
    d=d.merge(chemistry,on='source_id',validate='one_to_one').sort_values('source_id')
    d=d[d.mean_R_kpc.between(.5,9)&(d.mean_z_kpc.abs()<=1.5)&(d.FE_H>=-.5)&(d.ALPHA_M<.15)].reset_index(drop=True)
    y=np.stack([d['mean_'+k].to_numpy(float) for k in LABELS],axis=1)
    y[:,2]-=np.radians(P['bar_angle_degrees'])
    covariance=np.zeros((len(d),6,6))
    for i,k in enumerate(LABELS):
        for j in range(i+1):
            covariance[:,i,j]=covariance[:,j,i]=d['cov_'+k+'__'+LABELS[j]].to_numpy()
    eig=np.linalg.eigvalsh(covariance)
    assert eig.min()>-1e-8
    chosen,strata=anchors(y,d.source_id.to_numpy())
    x,p=cartesian(y[chosen]);launch=d.loc[chosen,['source_id','sky_pixel','FIELD','FE_H','ALPHA_M']].copy()
    for i,name in enumerate(['x_kpc','y_kpc','z_kpc','vx_kms','vy_kms','vz_kms']):
        launch[name]=np.c_[x,p][:,i]
    launch['R_stratum']=[s[0] for s in strata];launch['height_stratum']=[s[1] for s in strata]
    launch['holdout_role']='training';launch.to_parquet(OUT/'expanded-launches.parquet',index=False)
    summary=dict(training_stars=len(d),bulge_stars=int((y[:,0]<3.5).sum()),
                 expanded_training_launches=len(chosen),uncertainty_source_sha256=digest(src),
                 launches_sha256=digest(OUT/'expanded-launches.parquet'),models={})
    diagnostics=d[['source_id','sky_pixel']].copy()
    disk=DiskGrid(256,129,None)
    old=json.loads((ORBIT/'results.json').read_text())
    initial=np.array(old['initial_conditions']);R=np.hypot(initial[:,0],initial[:,1]);phi=np.arctan2(initial[:,1],initial[:,0]);c=np.cos(phi);s=np.sin(phi)
    demo=np.c_[R,initial[:,2],phi,initial[:,3]*c+initial[:,4]*s,-initial[:,3]*s+initial[:,4]*c,initial[:,5]]
    for model in P['models']:
        field=Field(disk,model)
        J,g=invariants(field,y)
        sigma=np.sqrt(np.maximum(np.einsum('ni,nij,nj->n',g,covariance,g),1e-12))
        centres,_=invariants(field,demo)
        expanded=J[chosen]
        diagnostics[model+'_J_kms2']=J;diagnostics[model+'_J_sigma_kms2']=sigma
        masks={'all':np.ones(len(d),bool),'bulge':y[:,0]<3.5,'disk_control':y[:,0]>=5}
        sample={name:dict(stars=int(mask.sum()),demo_support=support(J[mask],sigma[mask],centres),
                         expanded_launch_support=support(J[mask],sigma[mask],expanded)) for name,mask in masks.items()}
        summary['models'][model]=dict(regions=sample,
            demo_mixture_fits=[fit_weights(J,sigma,centres,b) for b in [0.,1000.]],
            expanded_mixture_fits=[fit_weights(J,sigma,expanded,b) for b in [0.,1000.]],
            median_approximate_J_sigma_kms2=float(np.median(sigma)),demo_J_centres_kms2=centres.tolist())
        save('results.json',summary)
        print(model,sample['all'],flush=True)
    diagnostics.to_parquet(OUT/'training-J-diagnostic.parquet',index=False)
    summary['diagnostics_sha256']=digest(OUT/'training-J-diagnostic.parquet')
    summary['heldout_scores_evaluated']=False
    summary['full_position_velocity_likelihood_fitted']=False
    summary['expanded_orbits_integrated']=False
    save('results.json',summary)
    print('Training orbit-support and positive-weight mixture diagnostics complete.',flush=True)

if __name__=='__main__':main()

"""Spatial holdout allocation and Gaia/StarHorse uncertainty propagation.

Uses known Gaussian conditioning and coordinate transformations, not a new
gravity law. StarHorse joint-posterior reconstruction is explicitly approximate.
"""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy.special import ndtri
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential
import astropy.units as u

H=Path(__file__).resolve().parent
ROOT=H.parents[2]
CACHE=ROOT/'research_work/data-cache/stellar-catalogs'
P=json.loads((H/'protocol.json').read_text())
save=lambda name,obj:(H/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def split(pixel):
    v=int.from_bytes(hashlib.sha256((P['split']['hash_salt']+':'+str(pixel)).encode()).digest()[:8],'big')%10
    return 'training' if v<6 else ('validation' if v<8 else 'test')

def coordinates(ra,dec,distance,pmra,pmdec,rv):
    gc=Galactocentric(galcen_distance=8.2*u.kpc,z_sun=.0208*u.kpc,
        galcen_v_sun=CartesianDifferential([11.1,248.,7.25]*u.km/u.s))
    c=SkyCoord(ra=ra*u.deg,dec=dec*u.deg,distance=distance*u.kpc,
        pm_ra_cosdec=pmra*u.mas/u.yr,pm_dec=pmdec*u.mas/u.yr,
        radial_velocity=rv*u.km/u.s).transform_to(gc)
    x=-c.x.to_value(u.kpc);y=c.y.to_value(u.kpc)
    vx=-c.v_x.to_value(u.km/u.s);vy=c.v_y.to_value(u.km/u.s)
    R=np.hypot(x,y)
    return np.c_[R,c.z.to_value(u.kpc),np.arctan2(y,x),
                 (x*vx+y*vy)/R,(x*vy-y*vx)/R,c.v_z.to_value(u.km/u.s)]

def region(R,z):
    a=np.abs(z);code=np.zeros(len(R),dtype=np.int8)
    code[(R>=.5)&(R<=3.5)&(a<.2)]=1
    code[(R>=.5)&(R<=3.5)&(a>=.5)&(a<=1.5)]=2
    code[(R>=5)&(R<=9)&(a<.2)]=3
    code[(R>=5)&(R<=9)&(a>=.5)&(a<=1.5)]=4
    return code

def main():
    src=CACHE/P['parent_catalog'];assert digest(src)==P['parent_sha256']
    d=pd.read_parquet(src);d=d[d.gaia_quality_candidate].reset_index(drop=True)
    assert d.source_id.is_unique and d.source_id.dtype.kind in 'iu'
    pixels=d.source_id.to_numpy()//P['split']['pixel_integer_divisor']
    roles=np.array([split(int(k)) for k in pixels]);d['holdout_role']=roles;d['sky_pixel']=pixels
    groups=d.groupby('sky_pixel').holdout_role.nunique();assert groups.max()==1
    allocation=d[['source_id','sky_pixel','holdout_role','FIELD']].copy()
    allocation.to_parquet(CACHE/'stellar-spatial-holdouts.parquet',index=False)
    save('split-groups.json',[{'sky_pixel':int(k),'role':split(int(k)),'stars':int(v)} for k,v in d.groupby('sky_pixel').size().items()])
    # All five astrometric uncertainties and all ten pairwise correlations.
    names=['ra','dec','parallax','pmra','pmdec'];n=len(d)
    err=np.stack([d[k+'_error'].to_numpy(float) for k in names],axis=1)
    corr=np.repeat(np.eye(5)[None,:,:],n,axis=0)
    for i in range(5):
        for j in range(i+1,5):
            corr[:,i,j]=corr[:,j,i]=d[names[i]+'_'+names[j]+'_corr'].to_numpy(float)
    assert np.isfinite(corr).all() and np.isfinite(err).all() and (err>0).all()
    eig=np.linalg.eigvalsh(corr)
    assert eig.min()>0, 'Do not silently repair invalid catalog covariance'
    cov=corr*err[:,:,None]*err[:,None,:]
    idx=[0,1,3,4]
    marginal=cov[:,idx,:][:,:,idx]
    cross=cov[:,idx,2]
    regression=cross/cov[:,2,2,None]
    conditional=marginal-np.einsum('ni,nj->nij',cross,cross)/cov[:,2,2,None,None]
    chol_ind=np.linalg.cholesky(marginal);chol_cond=np.linalg.cholesky(conditional)
    base=np.stack([d[k].to_numpy(float) for k in ['ra','dec','pmra','pmdec']],axis=1)
    dist=d.dist50.to_numpy(float);parallax=d.parallax.to_numpy(float)
    below=np.log(dist/d.dist16.to_numpy(float))/ndtri(.84)
    above=np.log(d.dist84.to_numpy(float)/dist)/ndtri(.84)
    assert (below>=0).all() and (above>=0).all()
    central=coordinates(base[:,0],base[:,1],dist,base[:,2],base[:,3],d.VHELIO_AVG.to_numpy(float))
    nominal_region=region(central[:,0],central[:,1])
    counts={role:{str(k):int(((roles==role)&(nominal_region==k)).sum()) for k in range(5)} for role in ['training','validation','test']}
    shared_fields=int((d.groupby('FIELD').holdout_role.nunique()>1).sum())
    summary={'sources':n,'roles':pd.Series(roles).value_counts().to_dict(),
       'sky_pixels':int(len(groups)),'fields_spanning_roles':shared_fields,
       'coverage_by_role_region':counts,'minimum_correlation_eigenvalue':float(eig.min()),
       'role_allocation_sha256':digest(CACHE/'stellar-spatial-holdouts.parquet'),
       'region_codes':{'0':'outside targets','1':'plane beneath bulge','2':'off-plane bulge','3':'disk plane control','4':'off-plane disk control'},
       'modes':{}}
    alpha=P['photon_prediction']['alpha_per_Mpc'];c=P['photon_prediction']['speed_of_light_kms']
    pred=d[['source_id','holdout_role']].copy()
    pred['D_kpc']=dist;pred['z_conversion']=np.expm1(alpha*dist/1000)
    pred['equivalent_los_kms']=c*pred.z_conversion
    pred['photon_energy_fraction_transferred']=pred.z_conversion/(1+pred.z_conversion)
    pred.to_parquet(CACHE/'stellar-frozen-redshift-predictions.parquet',index=False)
    summary['frozen_conversion_equivalent_los_kms_quantiles']=np.quantile(pred.equivalent_los_kms,[0,.16,.5,.84,1]).tolist()
    for mode in P['uncertainty']['modes']:
        rng=np.random.default_rng(P['uncertainty']['seed'])
        independent=mode=='independent'
        bias=0 if independent else (-float(mode.split('_')[-1]) if 'minus' in mode else float(mode.split('_')[-1]))
        chol=chol_ind if independent else chol_cond
        mean=np.zeros((n,6)); M2=np.zeros((n,6,6));migrated=np.zeros(n)
        draws=P['uncertainty']['draws_per_mode']
        for it in range(draws):
            q=rng.standard_normal(n)
            dd=dist*np.exp(q*np.where(q>=0,above,below))
            mu=np.zeros((n,4)) if independent else regression*(1/dd+bias-parallax)[:,None]
            delta=mu+np.einsum('nij,nj->ni',chol,rng.standard_normal((n,4)))
            ra=base[:,0]+delta[:,0]/(3600000*np.cos(np.radians(base[:,1])))
            dec=base[:,1]+delta[:,1]/3600000
            rv=d.VHELIO_AVG.to_numpy(float)+d.VERR.to_numpy(float)*rng.standard_normal(n)
            a=coordinates(ra,dec,dd,base[:,2]+delta[:,2],base[:,3]+delta[:,3],rv)
            # Unwrap azimuth locally to avoid an artificial error at +/-pi.
            a[:,2]=central[:,2]+(a[:,2]-central[:,2]+np.pi)%(2*np.pi)-np.pi
            old=a-mean;mean+=old/(it+1);M2+=np.einsum('ni,nj->nij',old,a-mean)
            migrated+=(region(a[:,0],a[:,1])!=nominal_region)
            if (it+1)%16==0:print(mode,it+1,'/',draws,flush=True)
        covariance=M2/(draws-1)
        labels=['R_kpc','z_kpc','phi_rad','vR_kms','vphi_kms','vz_kms']
        out=allocation.copy()
        for i,k in enumerate(labels):
            out['mean_'+k]=mean[:,i]
            for j in range(i+1):out['cov_'+k+'__'+labels[j]]=covariance[:,i,j]
        out['target_region_change_probability']=migrated/draws
        filename='stellar-errors-'+mode+'.parquet';out.to_parquet(CACHE/filename,index=False)
        sd=np.sqrt(np.maximum(0,np.diagonal(covariance,axis1=1,axis2=2)))
        ms={'file':filename,'sha256':digest(CACHE/filename),
           'median_uncertainty':dict(zip(labels,np.median(sd,axis=0).tolist())),
           'q84_uncertainty':dict(zip(labels,np.quantile(sd,.84,axis=0).tolist())),
           'median_velocity_shift_from_catalog_kms':np.median(np.abs(mean[:,3:]-central[:,3:]),axis=0).tolist(),
           'target_stars_region_change_fraction':float(np.mean(migrated[nominal_region>0]/draws))}
        summary['modes'][mode]=ms
        if mode=='conditional_minus_0.017':
            # Training summaries only; no held-out velocity scores are inspected.
            stats=[]
            for k in range(1,5):
                w=(roles=='training')&(nominal_region==k)
                for i,label in enumerate(labels[3:],start=3):
                    y=central[w,i]
                    stats.append({'region':k,'observable':label,'stars':int(w.sum()),
                         'catalog_mean':float(y.mean()),'catalog_std':float(y.std(ddof=1)),
                         'median_propagated_sd':float(np.median(sd[w,i])),
                         'warning':'Selected sample moments; not intrinsic dispersion or a force measurement'})
            save('training-moments.json',stats)
    save('results.json',summary)
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()

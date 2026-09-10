"""Training-only decomposition of the observed vertical-velocity proxy."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord,Galactocentric,CartesianDifferential

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache'; OUT=CACHE/'bulge-velocity-components'; OUT.mkdir(exist_ok=True)

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def moments(q):
    names=['rv_vertical','pm_vertical','vz_median_distance','VHELIO_AVG','mean_vz_kms']
    x=q[names].to_numpy(float)
    cov=np.cov(x,rowvar=False,ddof=0)
    identity=float(cov[0,0]+cov[1,1]+2*cov[0,1])
    assert abs(identity-cov[2,2])<1e-8
    return {'stars':len(q),'means_kms':dict(zip(names,x.mean(axis=0).tolist())),
            'stds_kms':dict(zip(names,np.sqrt(np.diag(cov)).tolist())),
            'vertical_variance_decomposition_kms2':{'rv':float(cov[0,0]),'proper_motion_distance':float(cov[1,1]),'twice_cross_covariance':float(2*cov[0,1]),'total':float(cov[2,2])},
            'median_abs_b_degrees':float(q.abs_b.median()),'survey_fields':int(q.FIELD.nunique())}

def main():
    paths={'parent':CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet',
       'moments':CACHE/'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet',
       'association_screen':CACHE/'stellar-association-window/training-screen.parquet'}
    cols=['source_id','ra','dec','pmra','pmdec','dist50','VHELIO_AVG','VERR','FE_H','ALPHA_M','FIELD']
    screen=pd.read_parquet(paths['association_screen'],columns=['source_id','holdout_role','window_residual_arcsec'])
    assert screen.holdout_role.eq('training').all()
    d=pd.read_parquet(paths['parent'],columns=cols).merge(screen,on='source_id',validate='one_to_one')
    m=pd.read_parquet(paths['moments'],columns=['source_id','mean_R_kpc','mean_z_kpc','mean_vz_kms'],filters=[('holdout_role','==','training')])
    d=d.merge(m,on='source_id',validate='one_to_one')
    selected=d.mean_R_kpc.between(.5,9)&(d.mean_z_kpc.abs()<=1.5)&(d.FE_H>=-.5)&(d.ALPHA_M<.15)
    d=d[selected].copy(); assert len(d)==27884
    frame=Galactocentric(galcen_distance=8.2*u.kpc,z_sun=.0208*u.kpc,galcen_v_sun=CartesianDifferential([11.1,248,7.25]*u.km/u.s))
    def transform(pm_on,rv_on):
        return SkyCoord(ra=d.ra.to_numpy()*u.deg,dec=d.dec.to_numpy()*u.deg,distance=d.dist50.to_numpy()*u.kpc,
          pm_ra_cosdec=d.pmra.to_numpy()*pm_on*u.mas/u.yr,pm_dec=d.pmdec.to_numpy()*pm_on*u.mas/u.yr,
          radial_velocity=d.VHELIO_AVG.to_numpy()*rv_on*u.km/u.s).transform_to(frame).v_z.to_value(u.km/u.s)
    zero=transform(0,0);rv=transform(0,1)-zero;pm=transform(1,0)-zero;total=transform(1,1)
    np.testing.assert_allclose(zero,7.25,rtol=0,atol=1e-10)
    error=float(np.max(abs(total-(zero+rv+pm))));assert error<1e-9
    d['rv_vertical']=rv;d['pm_vertical']=pm;d['vz_median_distance']=total
    gal=SkyCoord(d.ra.to_numpy()*u.deg,d.dec.to_numpy()*u.deg).galactic
    d['longitude']=gal.l.wrap_at(180*u.deg).degree;d['abs_b']=abs(gal.b.degree)
    masks={'bulge_plane':(d.mean_R_kpc<3.5)&(d.mean_z_kpc.abs()<.2),
           'bulge_offplane':(d.mean_R_kpc<3.5)&d.mean_z_kpc.abs().between(.5,1.5),
           'disk_plane':(d.mean_R_kpc>=5)&(d.mean_z_kpc.abs()<.2),
           'disk_offplane':(d.mean_R_kpc>=5)&d.mean_z_kpc.abs().between(.5,1.5)}
    results={}
    for name,mask in masks.items():
        part=d[mask];flag=part.window_residual_arcsec>.5
        variants={'all':moments(part),'flagged':moments(part[flag]) if flag.sum()>1 else None,'unflagged_or_unassessed':moments(part[~flag])}
        fields=[]
        for field,g in part.groupby('FIELD',observed=True):
            label=field.decode().strip() if isinstance(field,bytes) else str(field)
            fields.append({'field':label,'stars':len(g),'flagged':int((g.window_residual_arcsec>.5).sum())})
        variants['fields_with_flags']=sorted([x for x in fields if x['flagged']>0],key=lambda x:(-x['flagged'],-x['stars']))
        results[name]=variants
    # Common sky/metallicity cells compare flagged and other bulge-plane stars.
    # This is conditional descriptive balancing, not a survey-selection model.
    part=d[masks['bulge_plane']].copy()
    part['lbin']=np.floor(part.longitude/5).astype(int)
    part['bbin']=np.floor(part.abs_b/2).astype(int)
    part['metalbin']=np.floor(part.FE_H/.25).astype(int)
    cell_results=[]
    for key,g in part.groupby(['lbin','bbin','metalbin'],observed=True):
        f=g.window_residual_arcsec>.5
        if f.sum()>=3 and (~f).sum()>=3:
            cell_results.append({'cell':[int(x) for x in key],'flagged':moments(g[f]),'other':moments(g[~f])})
    balanced={}
    weights=np.array([min(c['flagged']['stars'],c['other']['stars']) for c in cell_results],float)
    weights/=weights.sum()
    for side in ['flagged','other']:
        values={}
        for key in ['rv_vertical','pm_vertical','vz_median_distance','VHELIO_AVG']:
            means=np.array([c[side]['means_kms'][key] for c in cell_results])
            variance=np.array([c[side]['stds_kms'][key]**2 for c in cell_results])
            mean=float(weights@means)
            values[key]={'weighted_mean_kms':mean,'within_cell_spread_kms':float(np.sqrt(weights@variance)),
                         'mixture_spread_kms':float(np.sqrt(weights@(variance+means**2)-mean**2))}
        balanced[side]={'stars':sum(c[side]['stars'] for c in cell_results),'values':values}
    d.to_parquet(OUT/'training-components.parquet',index=False)
    result={'selected_training_stars':len(d),'regions':results,'common_bulge_cells':cell_results,
       'common_cell_balanced_summary':balanced,'coordinate_sum_max_error_kms':error,'source_hashes':{k:sha(p) for k,p in paths.items()},
       'formula_status':'Known linear velocity-coordinate transformation and variance identity. Not a companion-gravity formula.',
       'scope':'Descriptive training sensitivity, no intrinsic-dispersion likelihood, no association correction, no gravity fit, no held-out outcomes.'}
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    for name,v in results.items():
        print(name,'all',v['all']['stds_kms'],'without flags',v['unflagged_or_unassessed']['stds_kms'])
    print('common cells',len(cell_results),'coordinate error',error)

if __name__=='__main__': main()

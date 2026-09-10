"""Training-only positional screen with a physically bounded observation epoch.

No source replacement, distance fit, gravity tuning or holdout scoring.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache'
OUT=CACHE/'stellar-association-window'
OUT.mkdir(exist_ok=True)

def digest(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def bounded_fit(offset,pm,lo,hi):
    denom=np.sum(pm*pm,axis=1)
    dt=np.divide(np.sum(offset*pm,axis=1),denom,out=np.zeros(len(pm)),where=denom>0)
    clipped=np.clip(dt,lo,hi)
    return dt,clipped,np.linalg.norm(offset-clipped[:,None]*pm,axis=1)

def main():
    paths={
      'parent':CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet',
      'roles':CACHE/'stellar-catalogs/stellar-spatial-holdouts.parquet',
      'moments':CACHE/'stellar-catalogs/stellar-errors-conditional_minus_0.017.parquet',
      'provenance':CACHE/'stellar-orbit-support/stellar-input-provenance.parquet'}
    roles=pd.read_parquet(paths['roles'],columns=['source_id','holdout_role'],filters=[('holdout_role','==','training')])
    cols=['source_id','APOGEE_ID','RA','DEC','ra','dec','pmra','pmdec','FE_H','ALPHA_M']
    d=pd.read_parquet(paths['parent'],columns=cols).merge(roles,on='source_id',validate='one_to_one')
    assert len(d)==77927 and d.source_id.is_unique
    g=SkyCoord(d.ra.to_numpy()*u.deg,d.dec.to_numpy()*u.deg)
    a=SkyCoord(d.RA.to_numpy()*u.deg,d.DEC.to_numpy()*u.deg)
    dx,dy=g.spherical_offsets_to(a)
    offset=np.c_[dx.arcsec,dy.arcsec]
    pm=d[['pmra','pmdec']].to_numpy()/1000
    assert np.isfinite(offset).all() and np.isfinite(pm).all()
    lo,hi=Time(['1997-06-07','2001-02-16']).jyear-2016
    dt,clip,res=bounded_fit(offset,pm,lo,hi)
    unlimited=np.linalg.norm(offset-dt[:,None]*pm,axis=1)
    eligible=d.APOGEE_ID.str.fullmatch(r'2M\d{8}[+-]\d{7}').to_numpy()
    d['window_residual_arcsec']=res
    d['unrestricted_residual_arcsec']=unlimited
    d['best_unrestricted_jyear']=2016+dt
    d['best_window_jyear']=2016+clip
    d['two_mass_designation']=d.APOGEE_ID.str.removeprefix('2M')
    epochs=g.spherical_offsets_by(pm[:,0]*clip*u.arcsec,pm[:,1]*clip*u.arcsec)
    exact=epochs.separation(a).arcsec
    assert np.max(abs(exact-res))<1e-4
    assert np.all(res>=unlimited-1e-10)
    d.loc[~eligible,'window_residual_arcsec']=np.nan
    # Independent brute grid must never improve the analytic bounded minimum.
    rng=np.random.default_rng(247); ii=rng.choice(len(d),1000,replace=False)
    grid=np.linspace(lo,hi,257)
    brute=np.min(np.linalg.norm(offset[ii,None,:]-grid[None,:,None]*pm[ii,None,:],axis=2),axis=1)
    assert np.all(brute>=res[ii]-1e-10)
    # Fixtures exercise lower/upper endpoint, interior, and zero proper motion.
    _,t,e=bounded_fit(np.array([[2.,1.],[-2.,0.],[.5,0.],[1.,2.]]),np.array([[1.,0.],[1.,0.],[1.,0.],[0.,0.]]),-1,1)
    np.testing.assert_allclose(t,[1,-1,.5,0]);np.testing.assert_allclose(e,[np.sqrt(2),1,0,np.sqrt(5)])
    d.to_parquet(OUT/'training-screen.parquet',index=False)
    # The crossmatch screen is defined without consulting velocities.
    moments=pd.read_parquet(paths['moments'],filters=[('holdout_role','==','training')])
    needed=['source_id']+[f'mean_{v}' for v in ['R_kpc','z_kpc','vR_kms','vphi_kms','vz_kms']]
    d=d.merge(moments[needed],on='source_id',validate='one_to_one')
    prov=pd.read_parquet(paths['provenance'],columns=['source_id','starhorse_used_parallax','parallax_interval_disjoint_5'])
    d=d.merge(prov,on='source_id',validate='one_to_one')
    selected=d.mean_R_kpc.between(.5,9)&(d.mean_z_kpc.abs()<=1.5)&(d.FE_H>=-.5)&(d.ALPHA_M<.15)
    masks={'bulge_plane':selected&(d.mean_R_kpc<3.5)&(d.mean_z_kpc.abs()<.2),
           'bulge_offplane':selected&(d.mean_R_kpc<3.5)&d.mean_z_kpc.abs().between(.5,1.5),
           'disk_plane':selected&(d.mean_R_kpc>=5)&(d.mean_z_kpc.abs()<.2),
           'disk_offplane':selected&(d.mean_R_kpc>=5)&d.mean_z_kpc.abs().between(.5,1.5)}
    thresholds=[.3,.5,1.,1.5]
    counts={str(x):{'window_flagged':int(((res>x)&eligible).sum()),'unrestricted_flagged':int(((unlimited>x)&eligible).sum()),'newly_flagged':int(((res>x)&(unlimited<=x)&eligible).sum())} for x in thresholds}
    regions={}
    vs=['mean_vR_kms','mean_vphi_kms','mean_vz_kms']
    for name,mask in masks.items():
        part=d[mask]; variants=[]
        for threshold in [None]+thresholds:
            q=part if threshold is None else part[(part.window_residual_arcsec<=threshold)|part.window_residual_arcsec.isna()]
            variants.append({'threshold_arcsec':threshold,'stars':len(q),'means_kms':q[vs].mean().tolist(),'stds_kms':q[vs].std(ddof=0).tolist()})
        regions[name]={'stars':len(part),'parallax_disagreement_stars':int(part.parallax_interval_disjoint_5.sum()),'flagged_0.5_arcsec':int((part.window_residual_arcsec>.5).sum()),'flagged_and_parallax_disagreement':int(((part.window_residual_arcsec>.5)&part.parallax_interval_disjoint_5).sum()),'sensitivity':variants}
    flagged=d[d.window_residual_arcsec>.5].copy()
    flagged[['source_id','APOGEE_ID','two_mass_designation','window_residual_arcsec','unrestricted_residual_arcsec','best_unrestricted_jyear']].assign(source_id=lambda x:x.source_id.astype(str)).to_csv(HERE/'review-candidates.csv',index=False,lineterminator='\n')
    result={'training_stars':len(d),'selected_training_stars':int(selected.sum()),'epoch_window_utc':['1997-06-07','2001-02-16'],'threshold_counts':counts,'regions':regions,'velocity_order':['vR','vphi','vz'],'max_spherical_vs_tangent_arcsec':float(np.max(abs(exact-res))),'source_hashes':{k:digest(p) for k,p in paths.items()},'screen_sha256':digest(OUT/'training-screen.parquet'),'verification':{'analytic_vs_grid':True,'synthetic_boundary_fixtures':True,'training_only':True,'unchanged_source_files':True},'interpretation':'Training positional diagnostic and descriptive sensitivity only. No adopted deletion, replacement, association probability or gravity ranking. Approximate distance-derived positions/moments retain previously documented provenance limitations.'}
    result['epoch_screen_applicable_training_stars']=int(eligible.sum())
    result['non_2mass_identifiers_not_epoch_screened']=int((~eligible).sum())
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'thresholds':counts,'regions':{k:{x:v[x] for x in ['stars','flagged_0.5_arcsec','flagged_and_parallax_disagreement']} for k,v in regions.items()}},indent=2))

if __name__=='__main__':main()

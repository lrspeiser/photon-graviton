"""Audit recovered observations and prepare an exploratory APOGEE/StarHorse sample.

No dynamical fit, distance inversion, or companion prediction is performed.
Run with Python, numpy, pandas, astropy, pyarrow and requests installed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential
import astropy.units as u

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache/stellar-catalogs'
URLS = {
    'APOGEE_DR17_EDR3_StarHorse_v2.fits': 'https://s3.data.aip.de:9000/shaqueiroz2023/APOGEE_DR17_EDR3_StarHorse_v2.fits',
    'allStarLite-dr17-synspec_rev1.fits': 'https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/allStarLite-dr17-synspec_rev1.fits',
}

def save(name, obj):
    (HERE/name).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8*1024*1024), b''): h.update(chunk)
    return h.hexdigest()

def frame(path, cols):
    t = Table.read(path, hdu=1)
    return t[cols].to_pandas()

def acquire():
    import requests
    CACHE.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        p = CACHE/name
        if p.exists(): continue
        with requests.get(url, stream=True, timeout=(30, 90)) as r:
            r.raise_for_status()
            with p.with_suffix('.part').open('wb') as f:
                for chunk in r.iter_content(4*1024*1024): f.write(chunk)
        # Only promote complete, parseable FITS files.
        Table.read(p.with_suffix('.part'), format='fits', hdu=1)
        p.with_suffix('.part').replace(p)

def main(local):
    records = []
    candidates = list((local/'data/gaia/slices').glob('raw_*.csv'))
    candidates += list((local/'data/gaia').glob('eilers_apogee*.csv'))
    candidates += [local/'data/gaia/6d_brava_full.fits',
                   local/'data/bulge_kinematics/BRAVA/brava_catalog.tbl',
                   local/'data/bulge_kinematics/GIBS/gibs_catalog_1.fits']
    for p in candidates:
        if not p.exists(): continue
        if p.suffix == '.csv':
            d = pd.read_csv(p, dtype={'source_id':'Int64'})
            cols = list(d.columns)
        else:
            d = Table.read(p, format='ascii.ipac' if p.suffix=='.tbl' else 'fits')
            cols = d.colnames
        records.append(dict(path=str(p), bytes=p.stat().st_size, sha256=digest(p),
                            rows=len(d), columns=cols,
                            status='Recovered; not adopted as a validated dynamical reduction'))
    save('local-inventory.json', records)
    files = []
    for name,url in URLS.items():
        p=CACHE/name
        files.append(dict(file=name, url=url, bytes=p.stat().st_size,sha256=digest(p)))
    save('download-manifest.json',files)
    sh=frame(CACHE/'APOGEE_DR17_EDR3_StarHorse_v2.fits',
             ['APOGEE_ID','DR3_source_id','dist16','dist50','dist84','StarHorse_OUTPUTFLAGS'])
    cols=['APOGEE_ID','GAIAEDR3_SOURCE_ID','RA','DEC','GLON','GLAT','FIELD','TELESCOPE',
          'SNR','NVISITS','VHELIO_AVG','VSCATTER','VERR','STARFLAG','ASPCAPFLAG',
          'FE_H','FE_H_ERR','ALPHA_M','ALPHA_M_ERR','TEFF','LOGG',
          'GAIAEDR3_PMRA','GAIAEDR3_PMDEC','GAIAEDR3_PMRA_ERROR','GAIAEDR3_PMDEC_ERROR']
    ap=frame(CACHE/'allStarLite-dr17-synspec_rev1.fits',cols)
    for d in [sh,ap]:
        d['APOGEE_ID']=d['APOGEE_ID'].map(lambda x:x.decode().strip() if isinstance(x,bytes) else str(x).strip())
    sh=sh.rename(columns={'DR3_source_id':'source_id'})
    ap=ap.rename(columns={'GAIAEDR3_SOURCE_ID':'source_id'})
    for d in [sh,ap]:
        assert d.source_id.dtype.kind in 'iu', 'Do not convert Gaia identifiers through floating point'
    # Reject ambiguous duplicate StarHorse pairs rather than arbitrarily selecting distances.
    keys=['APOGEE_ID','source_id']
    dup=sh.duplicated(keys,keep=False)
    counts={'starhorse_rows':len(sh),'apogee_rows':len(ap),
            'ambiguous_starhorse_rows_excluded':int(dup.sum())}
    sh=sh.loc[~dup & (sh.source_id>0)]
    ap=ap.loc[ap.source_id>0].sort_values(['SNR','FIELD','TELESCOPE'],ascending=[False,True,True])
    ap=ap.drop_duplicates(keys)
    j=sh.merge(ap,on=keys,how='inner',validate='one_to_one')
    counts['matched_pairs']=len(j)
    # A Gaia source must contribute at most once to exploratory moments.
    counts['multiple_apogee_ids_same_gaia_rows']=int(j.duplicated('source_id',keep=False).sum())
    j=j.sort_values(['SNR','APOGEE_ID'],ascending=[False,True]).drop_duplicates('source_id')
    counts['unique_matched_sources']=len(j)
    j['fractional_distance_halfwidth']=(j.dist84-j.dist16)/(2*j.dist50)
    # Astropy masks blank FITS strings; these become NaN in pandas.
    flags=j.StarHorse_OUTPUTFLAGS.fillna('').map(lambda x:x.decode().strip() if isinstance(x,bytes) else str(x).strip())
    cuts={
        'finite_positive_ordered_distance':np.isfinite(j.dist50)&(j.dist16>0)&(j.dist16<=j.dist50)&(j.dist50<=j.dist84),
        'distance_halfwidth_le_20pct':j.fractional_distance_halfwidth<=.2,
        'snr_ge_70':j.SNR>=70,
        'apogee_no_star_or_aspcap_flags':(j.STARFLAG==0)&(j.ASPCAPFLAG==0),
        'starhorse_no_output_warnings':flags.eq(''),
        'giants_logg_0_to_3_5':j.LOGG.between(0,3.5),
        'valid_rv':np.isfinite(j.VHELIO_AVG)&(j.VHELIO_AVG>-9990)&j.VERR.between(0,5,inclusive='right'),
        'stable_rv':j.VSCATTER.between(0,2),
        'valid_pm':np.isfinite(j.GAIAEDR3_PMRA)&np.isfinite(j.GAIAEDR3_PMDEC)&(j.GAIAEDR3_PMRA.abs()<9990)&(j.GAIAEDR3_PMDEC.abs()<9990),
        'valid_pm_errors':j.GAIAEDR3_PMRA_ERROR.between(0,1,inclusive='right')&j.GAIAEDR3_PMDEC_ERROR.between(0,1,inclusive='right'),
        'valid_chemistry':j.FE_H.between(-3,1)&j.ALPHA_M.between(-1,1),
    }
    keep=np.ones(len(j),dtype=bool);cutflow=[]
    for name,mask in cuts.items():
        keep &= mask.to_numpy();cutflow.append(dict(cut=name,remaining=int(keep.sum())))
    save('cutflow.json',cutflow)
    d=j.loc[keep].copy()
    # Declared frame convention, not a fit to the companion hypothesis.
    gc=Galactocentric(galcen_distance=8.2*u.kpc,z_sun=.0208*u.kpc,
                     galcen_v_sun=CartesianDifferential([11.1,248.,7.25]*u.km/u.s))
    c=SkyCoord(ra=d.RA.to_numpy(dtype=float)*u.deg,dec=d.DEC.to_numpy(dtype=float)*u.deg,
        distance=d.dist50.to_numpy(dtype=float)*u.kpc,pm_ra_cosdec=d.GAIAEDR3_PMRA.to_numpy(dtype=float)*u.mas/u.yr,
        pm_dec=d.GAIAEDR3_PMDEC.to_numpy(dtype=float)*u.mas/u.yr,radial_velocity=d.VHELIO_AVG.to_numpy(dtype=float)*u.km/u.s).transform_to(gc)
    # X points from Galactic centre toward Sun; Y follows disk rotation.
    x=-c.x.to_value(u.kpc);y=c.y.to_value(u.kpc);z=c.z.to_value(u.kpc)
    vx=-c.v_x.to_value(u.km/u.s);vy=c.v_y.to_value(u.km/u.s)
    R=np.hypot(x,y)
    d['R_kpc']=R;d['z_kpc']=z;d['phi_deg']=np.degrees(np.arctan2(y,x))
    d['bar_relative_phi_deg']=(d.phi_deg-27+180)%360-180
    d['vR_kms']=(x*vx+y*vy)/R;d['vphi_kms']=(x*vy-y*vx)/R
    d['vz_kms']=c.v_z.to_value(u.km/u.s)
    d['radius_region']=pd.cut(R,[.5,1.5,2.5,3.5,5,7,9,15,25]).astype(str)
    d['height_region']=pd.cut(np.abs(z),[0,.2,.5,1,1.5,3],include_lowest=True).astype(str)
    d['hemisphere']=np.where(z>=0,'above','below')
    d['metallicity_bin']=pd.cut(d.FE_H,[-3,-.5,0,1],include_lowest=True).astype(str)
    d['alpha_bin']=np.where(d.ALPHA_M>=.15,'alpha_ge_0.15','alpha_lt_0.15')
    d['bar_azimuth_bin']=pd.cut(d.bar_relative_phi_deg,[-180,-90,-45,0,45,90,180],include_lowest=True).astype(str)
    CACHE.mkdir(parents=True,exist_ok=True)
    d.to_parquet(CACHE/'matched-exploratory.parquet',index=False)
    groups=['radius_region','height_region','hemisphere','metallicity_bin','alpha_bin','bar_azimuth_bin']
    moments=[]
    for key,g in d.groupby(groups,observed=True):
        if len(g)<30 or 'nan' in key:continue
        row=dict(zip(groups,key));row['n']=len(g)
        for col in ['vR_kms','vphi_kms','vz_kms']:
            row['mean_'+col]=float(g[col].mean());row['std_'+col]=float(g[col].std())
        moments.append(row)
    save('exploratory-moments.json',moments)
    counts['exploratory_selected_sources']=len(d);counts['cells_with_at_least_30']=len(moments)
    regions={
       'plane_under_bulge':d.R_kpc.between(.5,3.5)&(d.z_kpc.abs()<.2),
       'above_below_bulge':d.R_kpc.between(.5,3.5)&d.z_kpc.abs().between(.5,1.5),
       'outer_disk_plane_control':d.R_kpc.between(5,9)&(d.z_kpc.abs()<.2),
       'outer_disk_off_plane_control':d.R_kpc.between(5,9)&d.z_kpc.abs().between(.5,1.5),
    }
    counts['target_region_counts']={k:int(v.sum()) for k,v in regions.items()}
    counts['status']='Exploratory point-estimate kinematics; no error deconvolution, selection correction or gravity-model prediction'
    counts['frame']={'galcen_distance_kpc':8.2,'z_sun_kpc':.0208,'solar_velocity_astropy_kms':[11.1,248.,7.25],'bar_angle_deg':27}
    counts['prepared_parquet_sha256']=digest(CACHE/'matched-exploratory.parquet')
    assert d.source_id.is_unique
    assert np.isfinite(d[['R_kpc','z_kpc','vR_kms','vphi_kms','vz_kms']]).all().all()
    save('readiness.json',counts)
    print(json.dumps(counts,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--local-project',type=Path,required=True)
    p.add_argument('--download',action='store_true')
    a=p.parse_args()
    if a.download:acquire()
    main(a.local_project)

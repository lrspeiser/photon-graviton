"""Training-only positional association diagnostic; never replaces identifiers."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

ROOT=Path(__file__).resolve().parents[3]
CACHE=ROOT/'research_work/data-cache'
OUT=CACHE/'stellar-orbit-support'
HERE=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
parent=CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet'
d=pd.read_parquet(parent)
roles=pd.read_parquet(CACHE/'stellar-catalogs/stellar-spatial-holdouts.parquet',columns=['source_id','holdout_role'])
d=d.merge(roles,on='source_id',validate='one_to_one')
d=d[d.holdout_role.eq('training')].copy()
g=SkyCoord(d.ra.to_numpy()*u.deg,d.dec.to_numpy()*u.deg)
a=SkyCoord(d.RA.to_numpy()*u.deg,d.DEC.to_numpy()*u.deg)
dx,dy=g.spherical_offsets_to(a)
offset=np.column_stack([dx.arcsec,dy.arcsec]); pm=d[['pmra','pmdec']].to_numpy()/1000
# Minimize separation over ALL possible epochs under constant proper motion.
# This avoids silently interpreting coordinate equinox J2000 as observation date.
dt=(offset*pm).sum(axis=1)/(pm*pm).sum(axis=1)
residual=np.linalg.norm(offset-dt[:,None]*pm,axis=1)
d['minimum_any_epoch_linear_residual_arcsec']=residual
d['best_linear_epoch_year']=2016+dt
d[['source_id','minimum_any_epoch_linear_residual_arcsec','best_linear_epoch_year']].to_parquet(OUT/'training-association-screen.parquet',index=False)
epochs=pd.read_csv(OUT/'four-2mass-epochs.csv',dtype={'designation':str})
neighbors=pd.read_csv(OUT/'four-gaia-neighbors.csv',dtype={'source_id':'int64'})
examples=[]
for e in epochs.itertuples():
    row=d[d.APOGEE_ID.eq('2M'+e.designation)].iloc[0]
    old=SkyCoord(e.ra*u.deg,e.dec*u.deg)
    elapsed=Time(e.jdate,format='jd').jyear-2016
    candidates=[]
    position_only=[]
    for n in neighbors.itertuples():
        c=SkyCoord(n.ra*u.deg,n.dec*u.deg)
        if c.separation(old).arcsec>5.0001: continue
        if not np.isfinite(n.pmra+n.pmdec):
            position_only.append(dict(source_id=str(n.source_id),uncorrected_offset_arcsec=float(c.separation(old).arcsec)))
            continue
        # Spherical tangent displacement, known constant-proper-motion approximation.
        earlier=c.spherical_offsets_by(n.pmra*elapsed*u.mas,n.pmdec*elapsed*u.mas)
        candidates.append(dict(source_id=str(n.source_id),epoch_corrected_offset_arcsec=float(earlier.separation(old).arcsec),
                               parallax_mas=float(n.parallax),parallax_error_mas=float(n.parallax_error)))
    candidates.sort(key=lambda x:x['epoch_corrected_offset_arcsec'])
    original=next(x for x in candidates if x['source_id']==str(row.source_id))
    examples.append(dict(APOGEE_ID=row.APOGEE_ID,original_source_id=str(row.source_id),
        observation_jyear=float(2016+elapsed),original_epoch_corrected_offset_arcsec=original['epoch_corrected_offset_arcsec'],
        minimum_any_epoch_offset_arcsec=float(row.minimum_any_epoch_linear_residual_arcsec),
        two_mass_position_error_major_arcsec=e.err_maj,starhorse_distance_kpc=float(row.dist50),candidates=candidates,
        candidates_without_proper_motion=position_only))
result=dict(training_stars=len(d),parent_sha256=sha(parent),
    training_minimum_any_epoch_residual_counts={str(t):int((residual>t).sum()) for t in [.3,.5,1,1.5]},
    examples=examples,source_files={p.name:sha(p) for p in [OUT/'four-2mass-epochs.csv',OUT/'four-gaia-neighbors.csv']},
    caveats=['Positional diagnostics, not formal association probabilities.',
    'Constant proper motion ignores perspective acceleration, orbital motion and annual parallax; no individual epoch covariance fit.',
    'Nearby candidates are not automatically adopted; crowding, blended photometry and spectra need review.',
    'No distances, identifiers, holdout roles or gravity parameters changed.'])
assert len(examples)==4 and len(d)==77927
assert all(e['original_epoch_corrected_offset_arcsec']>1 for e in examples)
(HERE/'association-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='examples'},indent=2))
for e in examples: print(e['APOGEE_ID'],e['original_epoch_corrected_offset_arcsec'],e['minimum_any_epoch_offset_arcsec'],e['candidates'][0])

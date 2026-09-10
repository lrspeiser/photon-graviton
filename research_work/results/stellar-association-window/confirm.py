"""Confirm the window-screen candidates against actual 2MASS epochs."""
from run import HERE,CACHE,OUT,digest
import json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

c=pd.read_csv(HERE/'review-candidates.csv',dtype={'source_id':'int64','two_mass_designation':str})
epochs=pd.read_csv(OUT/'2mass-review-epochs.csv',dtype={'designation':str})
c=c.merge(epochs,left_on='two_mass_designation',right_on='designation',validate='one_to_one')
parent=pd.read_parquet(CACHE/'stellar-catalogs/matched-with-gaia-covariance.parquet',columns=['source_id','ra','dec','RA','DEC','pmra','pmdec'])
c=c.rename(columns={'ra':'tm_ra','dec':'tm_dec'}).merge(parent,on='source_id',validate='one_to_one')
g=SkyCoord(c.ra.to_numpy()*u.deg,c.dec.to_numpy()*u.deg)
a=SkyCoord(c.RA.to_numpy()*u.deg,c.DEC.to_numpy()*u.deg)
tm=SkyCoord(c.tm_ra.to_numpy()*u.deg,c.tm_dec.to_numpy()*u.deg)
dt=Time(c.jdate.to_numpy(),format='jd').jyear-2016
prior=g.spherical_offsets_by(c.pmra.to_numpy()*dt*u.mas,c.pmdec.to_numpy()*dt*u.mas)
c['actual_epoch_gaia_to_2mass_arcsec']=prior.separation(tm).arcsec
c['apogee_to_2mass_arcsec']=a.separation(tm).arcsec
c['actual_observation_jyear']=dt+2016
assert len(c)==140 and c.source_id.is_unique
assert np.isfinite(c.actual_epoch_gaia_to_2mass_arcsec).all()
cols=['source_id','APOGEE_ID','actual_epoch_gaia_to_2mass_arcsec','apogee_to_2mass_arcsec','actual_observation_jyear','err_maj','err_min','window_residual_arcsec']
c[cols].assign(source_id=lambda x:x.source_id.astype(str)).to_csv(HERE/'confirmed-review.csv',index=False,lineterminator='\n')
q={'candidates':len(c),'actual_epoch_residual_above_0.5_arcsec':int((c.actual_epoch_gaia_to_2mass_arcsec>.5).sum()),'min_actual_epoch_residual_arcsec':float(c.actual_epoch_gaia_to_2mass_arcsec.min()),'max_apogee_to_2mass_arcsec':float(c.apogee_to_2mass_arcsec.max()),'max_2mass_error_major_arcsec':float(c.err_maj.max()),'actual_epoch_range_jyear':[float(min(dt)+2016),float(max(dt)+2016)],'epoch_input_sha256':digest(OUT/'2mass-review-epochs.csv'),'interpretation':'Positional consistency check, not mismatch probabilities or proof of a correct replacement. Constant proper motion omits annual parallax, perspective acceleration and binary/crowding effects.'}
(HERE/'confirmation.json').write_text(json.dumps(q,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(q,indent=2))

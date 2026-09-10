"""Independent coordinate fixtures and Monte Carlo error-propagation check."""
import json
import hashlib
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential
from run import HERE, ROOT, CACHE, P, coords, variances

solar=P['frame']
frame=Galactocentric(galcen_distance=solar['galcen_distance_kpc']*u.kpc,
    z_sun=solar['z_sun_kpc']*u.kpc,
    galcen_v_sun=CartesianDifferential(np.array(solar['v_sun_astropy_xyz_kms'])*u.km/u.s))
R=np.array([6.4,8.7,12.1,17.5]);phi=np.radians([-25,5,15,28]);z=np.array([.1,-.2,.3,-.4])
vr=np.array([-15,0,20,-5]);vp=np.array([220,235,240,225]);vz=np.array([5,-2,10,3])
x=-R*np.cos(phi);y=R*np.sin(phi)
vx=(x*vr+y*vp)/R;vy=(y*vr-x*vp)/R
fixture=SkyCoord(x=x*u.kpc,y=y*u.kpc,z=z*u.kpc,v_x=vx*u.km/u.s,v_y=vy*u.km/u.s,v_z=vz*u.km/u.s,frame=frame).icrs
d=pd.DataFrame(dict(ra=fixture.ra.deg,dec=fixture.dec.deg,distance_kpc=fixture.distance.to_value(u.kpc),
    pmra=fixture.pm_ra_cosdec.to_value(u.mas/u.yr),pmdec=fixture.pm_dec.to_value(u.mas/u.yr),radial_velocity=fixture.radial_velocity.to_value(u.km/u.s)))
expected=np.c_[R,np.degrees(phi),z,vr,vp,vz]
err=np.max(abs(coords(d)-expected),axis=0)
assert np.max(err)<1e-9

t=pd.read_parquet(CACHE/'cepheid-common-frame-train-selected.parquet').sort_values('source_id')
sample=t.iloc[np.linspace(0,len(t)-1,20,dtype=int)].reset_index(drop=True)
analytic=variances(sample)
rng=np.random.default_rng(2030910);n=4000
draws=sample.loc[sample.index.repeat(n)].reset_index(drop=True)
e=rng.normal(size=(len(draws),4));corr=draws.pmra_pmdec_corr.to_numpy()
delta=dict(pmra=e[:,0]*draws.pmra_error.to_numpy(),
    pmdec=(corr*e[:,0]+np.sqrt(1-corr**2)*e[:,1])*draws.pmdec_error.to_numpy(),
    rv=e[:,2]*draws.radial_velocity_error.to_numpy())
vel=coords(draws,scale=1+.07*e[:,3],delta=delta)[:,3:].reshape(20,n,3)
empirical=vel.var(axis=1,ddof=1)
relative=abs(empirical/analytic-1)
assert np.max(relative)<.15, relative
split=pd.read_parquet(CACHE/'cepheid-common-frame-split.parquet')
assert split.source_id.is_unique and set(t.source_id)<=set(split.loc[split.role.eq('train'),'source_id'])
bins=json.loads((HERE/'training-bins.json').read_text())
assert sum(row['n'] for row in bins)==len(t)
for name,expected_hash in P['frozen_files'].items():
    assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected_hash
out=dict(coordinate_fixture_max_abs_errors=err.tolist(),coordinate_columns=['R_kpc','phi_deg','z_kpc','vR_kms','vphi_kms','vz_kms'],
    monte_carlo_stars=20,draws_per_star=n,max_relative_variance_difference=float(relative.max()),
    median_relative_variance_difference=float(np.median(relative)),
    training_only_identifier_check=True,bin_count_check=True,frozen_hashes_unchanged=True,
    qualification='Monte Carlo checks the declared independent Gaussian 7-percent error scenario, not its astrophysical adequacy.')
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))

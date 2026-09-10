"""Check source matching and coordinate transformations independently of fits."""
import json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, Galactocentric, CartesianDifferential
from astropy.table import Table
import astropy.units as u
from prepare import CACHE, HERE, save, digest

d=pd.read_parquet(CACHE/'matched-exploratory.parquet')
assert d.source_id.is_unique and d.source_id.dtype.kind in 'iu'
phi=np.radians(d.phi_deg.to_numpy()); R=d.R_kpc.to_numpy()
x=R*np.cos(phi); y=R*np.sin(phi)
vx=d.vR_kms*np.cos(phi)-d.vphi_kms*np.sin(phi)
vy=d.vR_kms*np.sin(phi)+d.vphi_kms*np.cos(phi)
gc=Galactocentric(galcen_distance=8.2*u.kpc,z_sun=.0208*u.kpc,
    galcen_v_sun=CartesianDifferential([11.1,248.,7.25]*u.km/u.s))
c=SkyCoord(x=-x*u.kpc,y=y*u.kpc,z=d.z_kpc.to_numpy()*u.kpc,
    v_x=-vx.to_numpy()*u.km/u.s,v_y=vy.to_numpy()*u.km/u.s,
    v_z=d.vz_kms.to_numpy()*u.km/u.s,frame=gc).icrs
rv_error=float(np.max(np.abs(c.radial_velocity.to_value(u.km/u.s)-d.VHELIO_AVG)))
pm_error=float(np.max(np.abs(c.pm_ra_cosdec.to_value(u.mas/u.yr)-d.GAIAEDR3_PMRA)))
dist_error=float(np.max(np.abs(c.distance.to_value(u.kpc)-d.dist50)))
assert rv_error<1e-7 and pm_error<1e-7 and dist_error<1e-7
save('verification.json',{'unique_integer_ids':True,'rows':len(d),
    'roundtrip_max_rv_error_kms':rv_error,'roundtrip_max_pmra_error_masyr':pm_error,
    'roundtrip_max_distance_error_kpc':dist_error,'dynamical_prediction_tested':False})

p=CACHE/'gaia-dr3-quality-covariance.fits'
if p.exists():
    g=Table.read(p).to_pandas()
    assert g.source_id.dtype.kind in 'iu' and g.source_id.is_unique
    e=d.merge(g,on='source_id',validate='one_to_one')
    assert len(e)==len(d)
    delta=np.abs(e.pmra-e.GAIAEDR3_PMRA)
    # Catalog packaging may round single-precision values; report rather than assume equality.
    mask=(e.ruwe<1.4)&(e.visibility_periods_used>=9)&e.astrometric_params_solved.isin([31,95])&(~e.duplicated_source)&(e.non_single_star==0)
    e['gaia_quality_candidate']=mask
    e.to_parquet(CACHE/'matched-with-gaia-covariance.parquet',index=False)
    report={'matched_rows':len(e),'gaia_quality_candidates':int(mask.sum()),
        'maximum_pmra_packaging_difference_masyr':float(delta.max()),
        'quality_rule':'ruwe<1.4; visibility_periods_used>=9; astrometric_params_solved in (31,95); duplicated_source false; non_single_star=0',
        'caveat':'An undetected binary may pass; Gaia quality cuts are selection, not completeness corrections.',
        'enriched_file_sha256':digest(CACHE/'matched-with-gaia-covariance.parquet')}
    selected=e.loc[mask]
    report['target_region_counts']={
       'plane_under_bulge':int((selected.R_kpc.between(.5,3.5)&(selected.z_kpc.abs()<.2)).sum()),
       'above_below_bulge':int((selected.R_kpc.between(.5,3.5)&selected.z_kpc.abs().between(.5,1.5)).sum()),
       'outer_disk_plane_control':int((selected.R_kpc.between(5,9)&(selected.z_kpc.abs()<.2)).sum()),
       'outer_disk_off_plane_control':int((selected.R_kpc.between(5,9)&selected.z_kpc.abs().between(.5,1.5)).sum())}
    save('gaia-quality.json',report)
    print(json.dumps(report,indent=2))
print('Coordinate and identifier checks passed')

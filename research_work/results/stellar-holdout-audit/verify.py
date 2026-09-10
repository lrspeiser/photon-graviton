"""Verify spatial separation, uncertainty-file integrity and a conditioning identity."""
import numpy as np
import pandas as pd
from run import CACHE,H,P,save,digest,split

roles=pd.read_parquet(CACHE/'stellar-spatial-holdouts.parquet')
assert roles.source_id.is_unique
assert roles.groupby('sky_pixel').holdout_role.nunique().max()==1
assert all(role==split(int(pixel)) for pixel,role in zip(roles.sky_pixel,roles.holdout_role))
for mode in P['uncertainty']['modes']:
    d=pd.read_parquet(CACHE/('stellar-errors-'+mode+'.parquet'))
    assert np.array_equal(d.source_id,roles.source_id)
    for c in d.columns:
        if c.startswith(('mean_','cov_')):assert np.isfinite(d[c]).all()
    assert d.target_region_change_probability.between(0,1).all()
    for c in d.columns:
        if c.startswith('cov_') and len(set(c[4:].split('__')))==1:assert (d[c]>=0).all()
photons=pd.read_parquet(CACHE/'stellar-frozen-redshift-predictions.parquet')
z=photons.z_conversion.to_numpy()
loss=photons.photon_energy_fraction_transferred.to_numpy()
assert np.allclose(1/(1+z)+loss,1,rtol=0,atol=3e-16)
assert np.allclose(np.log1p(z),P['photon_prediction']['alpha_per_Mpc']*photons.D_kpc/1000,rtol=1e-13)
# Independent two-variable Gaussian experiment: posterior parallax mixed through
# a conditional proper motion must have the analytic mixture mean/variance.
rng=np.random.default_rng(921)
parallax=rng.normal(.7,np.sqrt(.3),200000)
pm=2+1.2*(parallax-.5)+rng.normal(0,np.sqrt(4-1.2**2),len(parallax))
expected_mean=2+1.2*(.7-.5);expected_var=4-1.2**2+1.2**2*.3
assert abs(pm.mean()-expected_mean)<.02
assert abs(pm.var()-expected_var)/expected_var<.02
save('verification.json',{'sources':len(roles),'sky_pixel_role_leakage':False,
     'conditioning_sample_mean':float(pm.mean()),'conditioning_expected_mean':expected_mean,
     'conditioning_sample_variance':float(pm.var()),'conditioning_expected_variance':expected_var,
     'photon_reservoir_fraction_identity':True,'full_astrophysical_energy_budget_proven':False,
     'held_out_orbital_scores_evaluated':False,'uncertainty_modes_checked':P['uncertainty']['modes']})
print('Spatial-allocation, covariance-file, conditioning, and per-photon bookkeeping checks passed.')

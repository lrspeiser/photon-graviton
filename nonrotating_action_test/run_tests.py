"""Restricted candidate checks. Existing data reused; no new blinded validation.
Run from an extracted archive with the three sibling data directories intact.
"""
from pathlib import Path
import contextlib, io, sys, json, hashlib
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar, brentq
from scipy.linalg import cho_solve
from scipy.stats import chi2

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'shared_interaction_test'))
with contextlib.redirect_stdout(io.StringIO()):
    import brightness as b
gamma=0.000077315/1e6
# Directly reported frequency-ratio drifts, not inferred alpha limits.
clocks=[]
for name,mu,sig in [('Yb_E3_over_E2',-6.8e-18,7.5e-18),('Yb_E3_over_Cs',-3.1e-17,3.4e-17)]:
    clocks.append(dict(ratio=name,observed_drift_per_year=mu,sigma_per_year=sig,
        zero_prediction_residual_in_sigma=-mu/sig,
        residual_response_coefficient_95_interval=[(mu-1.96*sig)/gamma,(mu+1.96*sig)/gamma]))
old_prediction=-2*gamma
old_ratio_to_sigma=abs(old_prediction)/clocks[1]['sigma_per_year']
# If a compensating response multiplies the hyperfine law by n^b,
# optical/hyperfine drift=-(2+b)*gamma in the leading fixed-optical approximation.
coeff=clocks[1]['residual_response_coefficient_95_interval']
comp_b=[-2-coeff[1],-2-coeff[0]]
# A pure linear secular signal lies exactly in the removed drift subspace.
t=np.linspace(0,2826942/31557600,1001)
X=np.column_stack([np.ones(len(t)),t]); signal=-gamma*t
detrended=signal-X@np.linalg.lstsq(X,signal,rcond=None)[0]

df=pd.read_csv(P.parent/'redshift_paper/all_164_groups.csv')
test=df[df['split']=='test']
if len(test)!=25:
    test=df[df['split']=='final_test']
assert len(test)==25,df['split'].value_counts()
null_rmse=float(np.sqrt(np.mean(test.observed_cmb_cz_kms.to_numpy()**2)))
original_rmse=float(np.sqrt(np.mean(test.time_residual_kms.to_numpy()**2)))

x=np.log1p(b.z)
def shape(a,sign):
    if a==0:return np.ones_like(x)
    y=a*x
    return np.sinc(y/np.pi) if sign==1 else np.sinh(y)/y
def objective(a,sign):
    r=b.r-5*np.log10(shape(a,sign))
    return float(r@cho_solve(b.vf,r))
fits=[]
for sign,upper in [(1,.999*np.pi/max(x)),(-1,10.)]:
    sol=minimize_scalar(lambda a:objective(a,sign),bounds=(0,upper),method='bounded',options={'xatol':1e-10})
    a,ch=min([(0.,objective(0.,sign)),(float(sol.x),float(sol.fun))],key=lambda v:v[1])
    fits.append(dict(sign=sign,a=a,chi_squared=ch,radius_Mpc=None if a==0 else 1/(b.kap*a)))
neg=fits[1];Rc=neg['radius_Mpc']
lo=brentq(lambda a:objective(a,-1)-neg['chi_squared']-3.841458820694124,0,neg['a'])
hi=brentq(lambda a:objective(a,-1)-neg['chi_squared']-3.841458820694124,neg['a'],2)
radius_interval=[1/(b.kap*hi),1/(b.kap*lo)]
zgrid=np.array([.1,.5,1.,2.])
R=np.log1p(zgrid)/b.kap
DA=Rc*np.sinh(R/Rc)
distances=[dict(z=float(z),R_Mpc=float(r),DA_Mpc=float(da),DL_Mpc=float((1+z)*da)) for z,r,da in zip(zgrid,R,DA)]
c=299792458.;G=6.67430e-11;pc=149597870700.*648000/np.pi
# Einstein-form static condition rho_energy+p=-c^4/(4*pi*G*Rc^2) for negative curvature.
enthalpy=-c**4/(4*np.pi*G*(Rc*1e6*pc)**2)
result=dict(status='Restricted analytic checks and published-summary/catalog reanalysis; no complete successful new action',
    gamma_per_year=gamma,clock_data=clocks,
    old_fixed_optical_prediction_per_year=old_prediction,
    old_prediction_to_reported_sigma=old_ratio_to_sigma,
    clock_caveat='Leading constitutive model is not an exact Yb/Cs atomic calculation; huge ratio is a scale comparison, not an exact exclusion significance.',
    compensating_hyperfine_exponent_b_95=comp_b,
    compensation_caveat='Fitting a response exponent to a slope is not deriving a protective symmetry or validating it independently.',
    cavity=dict(predicted_fixed_atom_fixed_length_drift_per_year=-gamma,
        expected_fractional_change_over_2826942_seconds=-gamma*2826942/31557600,
        ideal_linear_detrending_max_remaining=float(max(abs(detrended))),
        status='Published Kennedy analysis removes overall linear drift; no raw undetrended-data secular fit performed.',
        source='https://arxiv.org/abs/2008.08773'),
    universal_metric=dict(action='Einstein-Hilbert[q] + minimally coupled Standard Model[q] + canonical scalar[q]',
        homogeneous_metric='ds_q^2=-c0^2 dt^2/n(t)^2+a0^2 dSigma_k^2',
        result='dTau=dt/n removes lapse; fixed a0 gives zero cosmological redshift for comoving material observers.',
        original_25_test_objects_reused=True,null_zero_cosmological_redshift_RMSE_km_s=null_rmse,
        original_temporal_RMSE_km_s=original_rmse,
        caveat='Descriptive zero-total-redshift benchmark with no peculiar-velocity fit; not a likelihood for a model of galaxy velocities.'),
    brightness=dict(N=len(b.z),M_fixed=b.M,kappa_per_Mpc=b.kap,flat_chi_squared=b.s['chi_squared'],
        fixed_reference_chi_squared=b.sr['chi_squared'],fits=fits,
        negative_radius_conditional_profile_95_Mpc=radius_interval,
        flat_four_bin_trend_chi_squared=b.bchi,flat_four_bin_trend_dof=3,
        flat_four_bin_trend_p=float(chi2.sf(b.bchi,3)),
        distances=distances,
        caveat='Exploratory reuse of same Pantheon+ data and conventional released covariance/calibration. Radius interval conditional on fixed law; angular distances are predictions, not tested angular data.'),
    background=dict(required_negative_curvature_rho_plus_p_J_m3=float(enthalpy),
        required_effective_mass_density_kg_m3=float(enthalpy/c**2),
        conclusion='Einstein gravity with ordinary positive-enthalpy matter and a canonical homogeneous scalar cannot support this negative-curvature static background, for any scalar potential or cosmological constant.',
        scope='Modified gravity or additional interactions change the effective equations; they are not excluded by this restricted result.'),
    untested=['Derivative interaction without specified action','Nonlinear tensor-sharing completion with fixed atoms',
        'Species-specific exact atomic and solid-state calculations','Stable modified-gravity nonexpanding solution',
        'New angular-distance data fit','CMB temperature or polarization likelihood'])
(P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
pd.DataFrame(distances).to_csv(P/'conditional_distances.csv',index=False)
print(json.dumps(result,indent=2))

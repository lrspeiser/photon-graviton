"""Shared-action consequences in a restricted low-energy completion.
The clock comparison is leading-order, not a full many-electron calculation.
"""
from pathlib import Path
import json,math
P=Path(__file__).resolve().parent
gamma=0.000077315/1e6
mu=-3.1e-17;sigma=3.4e-17
c=299792458.;yr=31557600.;pc=149597870700.*648000/math.pi
dist=[]
for Rmpc in [30,40,50]:
    lightyears=Rmpc*1e6*pc/(c*yr); x=gamma*lightyears
    delay=(x+math.expm1(-x))/gamma
    dist.append(dict(R_Mpc=Rmpc,affine_photon_vs_metric_GW_delay_years=delay,
        convention='Simultaneous emission; metric GW speed c0; photon reception today; positive means GW arrives later'))
out=dict(gamma_per_year=gamma,
    constitutive=dict(photon_metric='diag(-1/n^2,1,1,1) in the matter rest frame',
        epsilon_relative='Z*n',mu_relative='n/Z',photon_speed_relative='1/n',
        fixed_Coulomb_choice='Z=1/n',fixed_Coulomb_epsilon=1,magnetic_mu='n^2',
        leading_hydrogen_optical_frequency='constant',
        leading_hyperfine_frequency='proportional to n^2 for fixed moments and wavefunctions'),
    clock=dict(source='https://arxiv.org/pdf/2010.06620',
        directly_reported_E3_over_Cs_drift_per_year=mu,sigma_per_year=sigma,
        leading_fixed_optical_model_drift_per_year=-2*gamma,
        measured_95_interval_for_net_coefficient_multiplying_gamma=[(mu-1.96*sigma)/gamma,(mu+1.96*sigma)/gamma],
        caveat='Uses directly measured ratio slope, not the published alpha interpretation. Prediction assumes leading hyperfine scaling and negligible optical magnetic corrections; no exact confidence exclusion without species-specific calculation.'),
    GW=dist,
    gravity=dict(flat_static_condition='rho_total=0 and rho_total+p_total=0',
        rho_plus_p_today='rho_b + 4*rho_gamma/3 + K*gamma^2 > 0 in the stated ordinary-matter, homogeneous rest-frame branch',
        caveat='At n=1 the photon metric coincides with g. Homogeneous vector U=0 and lambda=0 at this epoch. A constant potential or cosmological constant cannot cancel rho+p. Curvature or new gravity changes the equations and requires a new background solution.'))
(P/'interaction_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

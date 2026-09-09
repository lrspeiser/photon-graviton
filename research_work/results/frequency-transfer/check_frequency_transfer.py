"""Energy-conserving spectral counterexamples for a common mean loss law.

Kernels are stipulated, not derived from a microphysical action. Geometry,
momentum, clocks and capture must be supplied separately.
"""
from pathlib import Path
import os
import json,shutil
import numpy as np
from scipy.stats import poisson

OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'frequency-transfer'
OUT.mkdir(parents=True,exist_ok=True)
checks=[];cases=[]
for S in [1.1,2.,5.]:
    for eps in [.5,.1,.01]:
        mu=np.log(S)/eps
        n=np.arange(int(poisson.ppf(1-1e-14,mu))+1)
        w=poisson.pmf(n,mu)
        energy=np.exp(n*np.log1p(-eps))
        mean=float(w@energy)
        second=float(w@(energy**2))
        expected_mean=1/S
        expected_second=np.exp(mu*((1-eps)**2-1))
        error=max(abs(mean-expected_mean),abs(second-expected_second))
        assert error<1e-11 and abs(w.sum()-1)<1e-11
        # Energy carried into no-loss companions equals the telescoping photon loss.
        companion=float(w@(1-energy))
        assert abs(mean+companion-1)<1e-11
        cv=np.sqrt(np.expm1(eps*np.log(S)))
        direct=np.sqrt(second/mean**2-1)
        assert abs(cv-direct)<1e-10
        cases.append(dict(centroid_shift_factor=S,event_fraction=eps,mean_events=mu,photon_energy=mean,companion_energy=companion,photon_survival=1,relative_line_rms=float(cv),moment_check_error=error))
checks.append(dict(check='Poisson_distribution_moments_and_energy_partition',cases=len(cases),max_moment_error=max(c['moment_check_error'] for c in cases)))
for S in [1.1,2.,5.]:
    # Whole-photon removal has the same first energy moment but leaves survivors
    # at their original energy. The zero-energy mass is no longer a photon.
    survival=np.exp(-np.log(S))
    assert np.isclose(survival,1/S,rtol=1e-14)
checks.append(dict(check='whole_photon_removal_same_mean_energy_surviving_line_unshifted',cases=3))
rng=np.random.default_rng(20260909)
residual=[]
for _ in range(100):
    energy=1.;companion=0.
    for fraction in rng.uniform(0,.2,size=100):
        transferred=energy*fraction
        energy-=transferred;companion+=transferred
    residual.append(abs(energy+companion-1))
assert max(residual)<1e-13
checks.append(dict(check='individual_random_transfer_histories_conserve_energy',histories=100,max_residual=max(residual)))
targets=[]
for S in [1.1,2.,5.]:
    for velocity_rms_km_s in [1.,10.,100.]:
        delta=velocity_rms_km_s/299792.458
        epsilon_max=np.log1p(delta**2)/np.log(S)
        mean_events_min=np.log(S)/epsilon_max
        targets.append(dict(centroid_shift_factor=S,illustrative_added_rms_km_s=velocity_rms_km_s,maximum_event_energy_fraction=float(epsilon_max),minimum_mean_events=float(mean_events_min)))
result=dict(scope='Stipulated forward energy kernels, fixed path and speed, unchanged clocks; no microscopic rate, inverse process, angular transfer or observed line-width bound established',checks=checks,poisson_cases=cases,illustrative_width_requirements=targets,identifiability='Same mean photon/companion energy partition does not determine photon number, line profile or transient time dilation',time_mapping='For common fixed path R and speed v, t_arrive=t_emit+R/v; dt_arrive/dt_emit=1 in the stipulated clock frame')
(OUT/'frequency-transfer-checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(__file__,OUT/Path(__file__).name)
print(json.dumps(dict(check_groups=len(checks),poisson_cases=len(cases),z1_10kms=[t for t in targets if t['centroid_shift_factor']==2 and t['illustrative_added_rms_km_s']==10][0]),indent=2))

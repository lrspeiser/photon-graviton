from pathlib import Path
import json,math,hashlib
import numpy as np
from scipy.linalg import expm
from scipy.stats import binom
from scipy.constants import h,e
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
obs_path=ROOT/'research_work/results/electromagnetic-audit/observations.json'
rate_path=ROOT/'research_work/results/electromagnetic-audit/protocol.json'
obs=json.loads(obs_path.read_text())['radio_methanol'];alpha=json.loads(rate_path.read_text())['alpha_per_mpc']
c_kms=299792.458;factor=2*math.sqrt(2*math.log(2))
S=1+obs['independent_centroid_groups'][0]['z'];q=1/S;A=math.log(S)
nu_emit=obs['independent_centroid_groups'][0]['frequency_GHz']*1e9;nu_obs=nu_emit*q
maxerr=0;maxmoment=0
for n in [32,64,128]:
    Q=np.zeros((n+1,n+1));k=np.arange(n+1)
    for m in range(1,n+1):Q[m,m]=-m;Q[m-1,m]=m
    for depth in [.01,.1,A]:
        probability=expm(depth*Q)[:,n];survival=math.exp(-depth)
        analytic=binom.pmf(k,n,survival)
        maxerr=max(maxerr,float(np.max(abs(probability-analytic))))
        mean=float(k@probability);var=float((k-mean)**2@probability)
        maxmoment=max(maxmoment,abs(mean-n*survival),abs(var-n*survival*(1-survival)))
        assert abs(probability.sum()-1)<1e-11
        assert np.max(abs(probability-analytic))<1e-11
        assert abs(mean-n*survival)<1e-10 and abs(var-n*survival*(1-survival))<1e-10
        assert np.all(k+(n-k)==n)
variance=(obs['FWHM_km_s']/factor/c_kms)**2
limit=variance*nu_obs/(1-q)
cases=[]
for step_hz in [1,10,100,1000]:
    cases.append(dict(step_energy_over_h_hz=step_hz,step_energy_ev=h*step_hz/e,
        lower_line_conversion_fwhm_kms=factor*c_kms*math.sqrt(step_hz/nu_obs*(1-q)),
        higher_line_conversion_fwhm_kms=factor*c_kms*math.sqrt(step_hz/(60.531489e9*q)*(1-q)),
        mean_companion_steps=(nu_emit-nu_obs)/step_hz))
synthetic=[]
q100=math.exp(-alpha*30.660139)
for nu in [48.3724558e9,6e14]:
    synthetic.append(dict(emitted_frequency_hz=nu,distance_mpc=30.660139,mean_redshift=1/q100-1,
        step_energy_over_h_hz=10,mean_companion_steps=nu*(1-q100)/10,
        conversion_fwhm_kms=factor*c_kms*math.sqrt(10/(nu*q100)*(1-q100)),
        event_interval_ratio=1,required_shared_stretch=1/q100,energy_flux_ratio_unchanged_time=q100))
result=dict(status='Conditional energy-process and exposed linewidth test, not a completed field theory.',
    alpha_per_mpc=alpha,observed_factor=S,observed_line_frequency_hz=nu_obs,
    variance_budget_scenarios=[dict(fraction=f,max_step_energy_over_h_hz=f*limit,max_step_energy_ev=f*limit*h/e) for f in [1,.1]],
    step_cases=cases,synthetic=synthetic,
    checks=dict(max_distribution_error=maxerr,max_moment_absolute_error=maxmoment),
    limitations=['No source distance used; not independent rate validation','Finite-line profiles and covariance not fitted',
        'Energy probabilities do not establish optical phase coherence','No event stretching at unchanged stationary transport',
        'Photon zero-energy absorption boundary and continuum completion unresolved','Momentum, angular redistribution and clocks not derived'],
    sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [obs_path,rate_path,HERE/'protocol.md']})
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))

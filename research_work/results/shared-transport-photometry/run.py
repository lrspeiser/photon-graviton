from pathlib import Path
import json, math, csv, hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import lambertw
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=ROOT/'research_work/results/electromagnetic-audit/protocol.json'
alpha=json.loads(source.read_text())['alpha_per_mpc']
def planck(x): return x**3/np.expm1(x)
baseline=quad(planck,1e-5,80,epsabs=1e-10,epsrel=1e-11)[0]
assert abs(baseline-math.pi**4/15)<1e-10
photon0=quad(lambda x:planck(x)/x,1e-5,80,epsabs=1e-10)[0]
rows=[]
for D in [0,30.660139,40.7,100,1000]:
    S=math.exp(alpha*D)
    # y is observed frequency h*nu_obs/(k*T_emit).
    f=lambda y:planck(S*y)/S
    energy=quad(f,1e-5/S,80/S,epsabs=1e-10,epsrel=1e-11)[0]/baseline
    photon=quad(lambda y:f(y)/y,1e-5/S,80/S,epsabs=1e-10)[0]/photon0
    grid=np.geomspace(.001,40,300)
    # B_nu(nu,T/S), in fixed emitted-temperature units.
    lower_temperature=grid**3/np.expm1(S*grid)
    relative=np.max(abs(np.array([f(y) for y in grid])/(S*S*lower_temperature)-1))
    assert abs(energy-S**-2)<1e-10
    assert abs(photon-S**-1)<1e-10
    assert relative<1e-12
    rows.append(dict(distance_mpc=D,stretch=S,redshift=S-1,
        bolometric_flux_ratio=energy,photon_rate_ratio=photon,
        event_fluence_ratio=energy*S,temperature_inferred_kelvin=6000/S,
        thermal_shape_amplitude=S*S,thermal_shape_relative_error=float(relative),
        bolometric_magnitude_extra=5*math.log10(S),
        naive_candle_distance_mpc=D*S,
        line_fractional_width_ratio=1.0,extra_frequency_dependent_delay_seconds=0.0))
with (HERE/'predictions.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
Dapp=40.7
Dtrue=float(lambertw(alpha*Dapp).real/alpha)
root=brentq(lambda d:d*math.exp(alpha*d)-Dapp,0,Dapp,xtol=1e-12)
assert abs(root-Dtrue)<1e-10
result=dict(alpha_per_mpc=alpha,predictions=rows,
    sbf_bolometric_sensitivity=dict(apparent_distance_mpc=Dapp,
        inferred_static_path_mpc=Dtrue,naive_distance_bias_percent=100*(Dapp/Dtrue-1),
        source='https://arxiv.org/abs/1801.06080',status='Idealized calibration sensitivity, not an adopted correction'),
    checks=dict(planck_integral_error=abs(baseline-math.pi**4/15),
        lambert_vs_bracket_mpc=abs(root-Dtrue),energy_photon_rate_and_shape='pass'),
    missing=['physical endpoint clock coupling','local rods and clocks','causal interaction and momentum closure',
        'derived spatial coupling','passband and source calibration likelihood','fresh observed spectra/flux test'],
    sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,HERE/'protocol.md']})
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))

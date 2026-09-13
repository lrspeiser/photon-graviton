"""Spectral discrimination of the previously bolometrically degenerate histories."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad

OUT=Path(__file__).resolve().parent;path=OUT/'history-results.json'
histories=json.loads(path.read_text())
# Exact SI constants; 6000 K is an explicit illustrative emission temperature.
h=6.62607015e-34;c=299792458.;kb=1.380649e-23;temp=6000.
normalization=np.pi**4/15
bands=[('below_0.4_um',0.,.4),('0.4_to_0.7_um',.4,.7),('0.7_to_5_um',.7,5.),('5_to_1000_um',5.,1000.),('above_1000_um',1000.,np.inf)]

def planck_fraction(lo,hi):
    # Wavelength bounds become reversed dimensionless frequency bounds.
    ymin=0. if np.isinf(hi) else h*c/(hi*1e-6*kb*temp)
    ymax=np.inf if lo==0 else h*c/(lo*1e-6*kb*temp)
    if ymin>750:return 0.
    def kernel(y):
        if y==0:return 0.
        if y>700:return y**3*np.exp(-y)
        return y**3/np.expm1(y)
    return quad(kernel,ymin,min(ymax,750),epsabs=1e-13,epsrel=1e-10,points=[p for p in [1.,3.,10.,30.] if ymin<p<min(ymax,750)])[0]/normalization

def emitted_band_at_age(age,lo,hi):
    # lambda_obs=lambda_emit*exp(age), with surviving photon energy exp(-age).
    scale=np.exp(-age)
    return scale*planck_fraction(lo*scale,hi*scale)

cases=[]
for case in histories['cases']:
    ages=case['burst_ages_alpha_c_T'];width=case['burst_width_alpha_c_T']
    matrix=np.array([[quad(lambda age:emitted_band_at_age(age,lo,hi),t-width/2,t+width/2,epsabs=1e-14)[0]/width for t in ages] for _,lo,hi in bands])
    original=matrix@np.array(case['baseline_emitted_energies'])
    altered=matrix@np.array(case['alternative_emitted_energies'])
    total=case['baseline_current_photon_companion_deposit'][0]
    assert abs(sum(original)/total-1)<1e-9
    assert abs(sum(altered)/total-1)<1e-9
    cases.append(dict(old_burst_age_alpha_c_T=case['old_burst_age_alpha_c_T'],
        bands=[dict(name=name,reference_energy=float(p),alternative_energy=float(q),difference_over_reference_total=float((q-p)/total)) for (name,_,_),p,q in zip(bands,original,altered)],
        spectral_total_variation_fraction=float(np.sum(abs(altered-original))/(2*total))))
result=dict(scope='Synthetic identical 6000 K blackbody emitters and unchanged positive histories; no dust, stellar evolution or observed background fit. Finite-band spectral differences do not establish uniqueness of arbitrary histories.',
    history_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),emission_temperature_K=temp,
    transport='Photon number conserved; each photon frequency and energy multiplied by exp(-alpha*c*age). This is a homogeneous radiation snapshot, not an event-arrival-time or luminosity-distance calculation.',cases=cases)
(OUT/'spectrum-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(result,indent=2))

"""Solar photon energy transfer and maximum locally retained mass equivalent."""
import json,hashlib
from pathlib import Path
import numpy as np
from scipy.special import zeta
HERE=Path(__file__).resolve().parent
p=HERE.parent/'brightness-distance-consistency/regular-area-results.json'
alpha=json.loads(p.read_text())['alpha_per_mpc']
L=3.828e26;T=5772.;h=6.62607015e-34;k=1.380649e-23;c=299792458.
AU=149597870700.;Mpc=3.085677581491367e22;Rsun=6.957e8;year=365.25*86400
mean=np.pi**4/(30*zeta(3))*k*T
rows=[]
for radius in [1.,30.,100.,1000.]:
    length=radius*AU-Rsun
    f=-np.expm1(-alpha*length/Mpc)
    power=L*f;rate=power/c**2
    rows.append(dict(radius_AU=radius,path_from_photosphere_m=length,
                     converted_energy_fraction=f,converted_power_W=power,
                     full_photon_energy_equivalents_per_s=power/mean,
                     max_retained_mass_equivalent_kg_per_s=rate,
                     max_retained_mass_equivalent_kg_per_year=rate*year))
out=dict(status='Sun-dominated outgoing luminosity; blackbody photon-count approximation; conditional conversion and perfect-capture upper bounds',
         alpha_input_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),alpha_per_Mpc=alpha,
         nominal_solar_luminosity_W=L,nominal_effective_temperature_K=T,
         mean_photon_energy_J=mean,mean_photon_energy_eV=mean/1.602176634e-19,
         emitted_photons_per_s=L/mean,
         full_conversion_mass_equivalent_kg_per_s=L/c**2,
         full_conversion_mass_equivalent_kg_per_year=L/c**2*year,
         average_photon_mass_equivalent_kg=mean/c**2,
         example_graviton_frequency_Hz=100.,
         example_gravitons_per_fully_converted_average_photon=mean/(h*100),
         rows=rows,
         sources=['https://arxiv.org/abs/1510.07674','https://www.iau.org/common/Uploaded%20files/IAUGA2015-Resolution-B3-recommended-nominal-conversion.pdf'])
assert all(0<r['converted_power_W']<L for r in rows)
assert abs(mean*(L/mean)/L-1)<1e-12
(HERE/'solar-photon-budget-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))

"""Physical capacity of a thermal linear-dispersion bosonic bath."""
from pathlib import Path
import json,math
from scipy.constants import Boltzmann as kb,hbar,c,electron_volt,Stefan_Boltzmann
P=Path(__file__).resolve().parent;MSUN=1.98847e30;KPC=3.085677581491367e19;YEAR=365.25*86400
old=json.loads((P/'coupled-torque-results.json').read_text())
coef=math.pi**2*2*kb**4/(30*hbar**3*c**3)
assert abs(coef/(4*Stefan_Boltzmann/c)-1)<1e-9
out=dict(scope='Two thermal bosonic polarizations with speed c, zero chemical potential; local storage and thermal outlet benchmarks, no universal bath exclusion',rows=[])
for b in ['I','II']:
 source=next(r for r in old['rows'] if r['baryons']==b and r['receiver_region_kpc']==[30,60])
 mass=source['mass_Msun'];rest=mass*MSUN*c*c
 for Ea in [1e-8,2e-10]:
  for eps in [.1,.5,.9]:
   Ti=.01*Ea*electron_volt/kb;Tmax=eps*Ea*electron_volt/(kb*math.log(1000))
   du=coef*(Tmax**4-Ti**4);assert du>0
   for name,Q in [('previous_orbital_cooling',source['released_energy_J']),('all_deposit_energy_from_protected_excitations',rest*eps/(1-eps))]:
    V=Q/du;R=(3*V/(4*math.pi))**(1/3)
    row=dict(baryons=b,bright_gap_eV=Ea,relaxation_fraction=eps,initial_temperature_K=Ti,max_temperature_K=Tmax,energy_case=name,heat_J=Q,inventory_mass_Msun=mass,required_bath_radius_kpc=R/KPC,benchmarks=[])
    for rk in [30.,120.,1000.]:
     radius=rk*KPC;volume=4*math.pi*radius**3/3;capacity=du*volume;ratio=Q/capacity
     luminosity=math.pi*radius**2*c*du
     row['benchmarks'].append(dict(radius_kpc=rk,thermal_capacity_J=capacity,heat_over_capacity=ratio,thermal_outlet_net_power_W=luminosity,minimum_thermal_export_years=Q/luminosity/YEAR,required_mode_speed_m_s=c/ratio**(1/3),required_polarization_count=2*ratio))
     assert abs(capacity*ratio/Q-1)<1e-12
    out['rows'].append(row)
(P/'bath-budget-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([r for r in out['rows'] if r['baryons']=='I' and r['relaxation_fraction']==.5],indent=2))

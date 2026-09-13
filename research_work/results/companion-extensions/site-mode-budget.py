"""Necessary steady-storage lifetime under a finite-site bosonic escape closure."""
from pathlib import Path
import json, math
from scipy.constants import hbar,c,electron_volt
from scipy.integrate import quad
P=Path(__file__).resolve().parent
MSUN=1.98847e30; KPC=3.085677581491367e19; YEAR=365.25*86400
source=json.loads((P/'coupled-torque-results.json').read_text())
out=dict(scope='Conditional entire fitted deposit energy in finite excitations; necessary steady-state bound, not a predicted lifetime',rows=[])
for b in ['I','II']:
 row=next(r for r in source['rows'] if r['baryons']==b and r['receiver_region_kpc']==[30,60])
 energy=row['mass_Msun']*MSUN*c*c
 for Ea in [1e-8,2e-10]:
  for radius in [30.,120.,1000.]:
   V=4*math.pi*(radius*KPC)**3/3; tau=3*radius*KPC/(4*c)
   for bandwidth in [.01,1.]:
    e0=.5*Ea*electron_volt; lo=e0*math.exp(-bandwidth/2);hi=e0*math.exp(bandwidth/2)
    modes=2*V*(hi**3-lo**3)/(6*math.pi**2*hbar**3*c**3)
    integral=quad(lambda x:x*x,math.exp(-bandwidth/2),math.exp(bandwidth/2))[0]
    numerical=2*V*e0**3*integral/(2*math.pi**2*hbar**3*c**3)
    assert abs(numerical/modes-1)<1e-12
    for target in [.75,.9,.99]:
     # Max energy per site with Pb>=target and Pa<=1-Pb.
     sites_min=energy/(Ea*electron_volt*(1-.5*target))
     mu_min=sites_min/modes
     n_cap=(1-target)/(2*target-1)
     lifetime=mu_min*tau*target/n_cap
     # Above the bound on gamma, no Pa<=1-Pb can provide positive replenishment.
     gamma=1.01/lifetime; n=mu_min*tau*gamma*target
     max_flow=(1+n)*(1-target)-n*target
     assert max_flow<0
     # Scaling by radius and excitation energy independently.
     out['rows'].append(dict(baryons=b,mass_Msun=row['mass_Msun'],bright_gap_eV=Ea,radius_kpc=radius,log_bandwidth=bandwidth,target_protected_fraction=target,sites_lower_bound=sites_min,coupled_modes=modes,sites_per_mode_lower_bound=mu_min,escape_years=tau/YEAR,occupation_upper_bound=n_cap,protected_lifetime_lower_bound_years=lifetime/YEAR))
for row in out['rows']:
 match=next(r for r in out['rows'] if all(r[k]==row[k] for k in ['baryons','bright_gap_eV','log_bandwidth','target_protected_fraction']) and r['radius_kpc']==120)
 ratio=row['protected_lifetime_lower_bound_years']/match['protected_lifetime_lower_bound_years']
 assert abs(ratio/(120/row['radius_kpc'])**2-1)<1e-12
out['checks']={'cases':len(out['rows']),'mode_quadrature_relative_tolerance':1e-12,'radius_scaling_relative_tolerance':1e-12,'above_bound_maximum_flow_negative':True}
(P/'site-mode-budget-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
for r in out['rows']:
 if r['baryons']=='I' and r['radius_kpc']==120 and r['target_protected_fraction']==.9:
  print(r)

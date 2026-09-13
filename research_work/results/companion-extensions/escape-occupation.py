"""Optically thin source occupancy for escaping, nonthermal relaxation waves."""
from pathlib import Path
import json,math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import hbar,c,electron_volt
P=Path(__file__).resolve().parent;KPC=3.085677581491367e19;YEAR=31557600.
q=json.loads((P/'bath-budget-results.json').read_text())
# Independent uniform-volume, isotropic-direction integration of path to a unit sphere.
x,w=leggauss(256);r=(x+1)/2;wr=w/2;mu=x
distance=-r[:,None]*mu+np.sqrt(1-r[:,None]**2*(1-mu**2))
mean=float(np.sum(3*r*r*wr*(distance@w/2)))
assert abs(mean-.75)<1e-8
out=dict(scope='Uniform isotropic emission, no absorption, constant occupation across a stipulated log-energy band; diagnostic before feedback',mean_escape_path_over_radius=mean,target_up_down_ratio=.001,rows=[])
limit=.001/.999
for source in q['rows']:
 if source['relaxation_fraction']!=.5:continue
 energy=.5*source['bright_gap_eV']*electron_volt
 for rk in [30.,120.,1000.]:
  R=rk*KPC;V=4*math.pi*R**3/3;tesc=.75*R/c
  for band in [.01,1.]:
   lo=energy*math.exp(-band/2);hi=energy*math.exp(band/2)
   u_unit=2*(hi**4-lo**4)/(8*math.pi**2*hbar**3*c**3)
   time_min=source['heat_J']*tesc/(V*u_unit*limit)
   examples=[]
   for years in [1e6,1e9,1e12]:
    duration=years*YEAR;L=source['heat_J']/duration
    dwell=(R/c)*float(np.sum(3*r*r*wr*(np.minimum(distance,duration*c/R)@w/2)))
    u=L*dwell/V;n=u/u_unit
    examples.append(dict(emission_duration_years=years,power_W=L,steady_field_reached=bool(duration>=2*R/c),effective_residence_years=dwell/YEAR,mean_mode_occupation=n,up_over_down=n/(1+n)))
   out['rows'].append(dict(baryons=source['baryons'],gap_eV=source['bright_gap_eV'],energy_case=source['energy_case'],heat_J=source['heat_J'],radius_kpc=rk,log_bandwidth=band,escape_time_years=tesc/YEAR,required_emission_duration_years=time_min/YEAR,examples=examples))
(P/'escape-occupation-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([r for r in out['rows'] if r['baryons']=='I' and r['gap_eV']==1e-8 and r['radius_kpc']==120 and r['energy_case']=='previous_orbital_cooling'],indent=2))

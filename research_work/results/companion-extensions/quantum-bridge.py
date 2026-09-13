"""Conditional bridge between discrete redshift losses and fitted bound mass."""
from pathlib import Path
import json,math
from scipy.constants import h,electron_volt
P=Path(__file__).resolve().parent
mass=json.loads((P/'settling-results.json').read_text())['selected']['mass_eV']
conversion=json.loads((P/'conversion-results.json').read_text())
alpha=conversion['alpha_per_Mpc'];D=30.660139;A=alpha*D
energies=[('1 GHz radio',h*1e9/electron_volt),('100 GHz microwave',h*1e11/electron_volt),('0.1 eV infrared',.1),('2 eV optical',2.),('10 eV ultraviolet',10.),('1 keV X-ray',1000.),('1 MeV gamma',1e6)]
out=dict(scope='Illustrative linewidth design constraints under independent fixed-fraction Poisson losses; effective mass is not an established creation gap',distance_Mpc=D,alpha_per_Mpc=alpha,mean_energy_loss_fraction=-math.expm1(-A),effective_bound_mass_eV=mass,rows=[])
for width in [1e-3,1e-5,1e-6]:
 delta=math.log1p(width**2)/A;events=A/delta
 assert 0<delta<1 and abs(math.expm1(events*delta**2)-width**2)<1e-15
 for label,E in energies:
  step=E*delta;loss=E*(-math.expm1(-A))
  row=dict(band=label,photon_energy_eV=E,illustrative_width=width,max_first_step_loss_eV=step,max_fraction_per_step=delta,minimum_mean_events=events,whole_path_mean_loss_eV=loss)
  row['hypothetical_gaps']=[]
  for factor in [1,2]:
   gap=factor*mass
   row['hypothetical_gaps'].append(dict(gap_eV=gap,minimum_first_step_contributions=math.ceil(gap/step),whole_path_photon_equivalent=gap/loss,one_step_can_supply_gap=bool(step>=gap),binding_fraction_if_gap_reduced_to_step=max(0,1-step/gap)))
  out['rows'].append(row)
# Exact Poisson moments from the probability generating function.
# A fixed absolute gap instead has fractional mean-loss rate Gamma*gap/E.
out['fixed_gap_rate_requirement']='Gamma(E)*gap/E=alpha, hence Gamma(E)=alpha*E/gap for fixed positive gap; this tunes the mean only, not fluctuations, thresholds or duration stretch.'
out['massive_travel_examples']=[]
for ratio in [2,10,1000,1e7]:
 beta=math.sqrt(1-ratio**-2)
 deficit=ratio**-2/(1+beta)
 out['massive_travel_examples'].append(dict(total_energy_over_rest=ratio,total_energy_eV=mass*ratio,speed_deficit_fraction=deficit))
(P/'quantum-bridge-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([r for r in out['rows'] if r['illustrative_width']==1e-5],indent=2))

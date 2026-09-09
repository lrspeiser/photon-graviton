from pathlib import Path
import numpy as np,json
from scipy.integrate import solve_ivp
p=.26390611655880836;gamma=7.7315e-11
rows=[]
for beta in [.1,1.]:
 for du,dv in [(0,0),(.1,0),(-.1,0),(.5,0),(-.5,0),(0,.1),(0,-.1)]:
  ni=.3;end=(1.5**(1-p)-ni**(1-p))/(1-p)
  def rhs(t,y):
   n,v,u=y;return [v,p*n**(2*p-1)-beta/n**2+u/n,-v*u/n]
  def boundary(t,y):return y[0]-.05
  boundary.terminal=True;boundary.direction=-1
  s=solve_ivp(rhs,[0,end],[ni,ni**p*(1+dv),beta*(1+du)/ni],events=boundary,rtol=1e-10,atol=1e-12)
  n,v,u=s.y;energy=.5*v*v-.5*n**(2*p)-beta/n+u
  rows.append(dict(beta=beta,radiation_fraction_change=du,initial_velocity_fraction_change=dv,hit_lower_boundary=bool(len(s.t_events[0])),n_final=float(n[-1]),target_n_at_same_time=float((ni**(1-p)+(1-p)*s.t[-1])**(1/(1-p))),final_speed_relative_to_target_law=float(v[-1]/n[-1]**p-1),energy_initial=float(energy[0]),energy_final=float(energy[-1]),max_energy_drift=float(max(abs(energy-energy[0])))))
opt=[]
for Om in [3.,10.,100.]:
 F=3*(Om*Om-.25);rr=[]
 for kin in [.5,1.,2.]:
  A=Om*Om+F+kin*kin;y=2*kin*kin*Om*Om/(A+np.sqrt(A*A-4*kin*kin*Om*Om));w=np.sqrt(y)
  rr.append(dict(initial_frequency=kin,final_frequency=float(w),frequency_ratio=float(w/kin)))
 epsdc=1+F/Om**2
 opt.append(dict(resonance=Om,oscillator_strength=F,mode_scan=rr,dc_permittivity=epsdc,leading_fixed_mass_charge_optical_energy_ratio=epsdc**-2))
sec=31557600.;pred=-gamma/sec;nu=299792458/(1542e-9)
rates=[-2.6e-19,-7.7e-20,-1e-19,5.2e-20]
cav=[dict(name=n,reported_fractional_drift_per_second=r,prediction_over_observed_magnitude=abs(pred/r),required_additional_drift_per_second=r-pred) for n,r in zip(['Si2','Si3','Si5','Si6'],rates)]
# Common offset is unidentifiable if all four instrument drift offsets are free.
A=np.column_stack([np.ones(4),np.eye(4)])
out=dict(status='Mixed published-summary comparison and synthetic tests; no raw cavity fit',feedback=rows,
 feedback_note='Potential fixed at original beta. No damping or retuning. Conserved energy and photon invariant prevent a unique full-state attractor; finite-time target-law deviations tested.',
 optical_medium=opt,optical_note='Lossless Lorentz dispersion endpoint test with externally prescribed oscillator-strength change, conserved k, adiabatic following of lower branch; not a full driven-medium energy calculation or pulse test.',
 cavity_prediction_per_second=pred,cavity_prediction_microHz_per_second_at_1542nm=pred*nu*1e6,cavity_summary=cav,
 cavity_source='Lee et al. arXiv:2509.13503v1, long-term drift section, https://arxiv.org/html/2509.13503v1',
 cavity_nuisance_design_rank=int(np.linalg.matrix_rank(A)),cavity_nuisance_parameters=5,
 caveats=['Reported recent drift values lack uncertainties/covariance here; no exclusion sigma assigned.', 'A common predicted cavity term could be masked by instrument drift; required offsets calculated, not justified.', 'No independent brightness or angular-distance data added.', 'No complete atomic protection, driver, or gravitational solution constructed.'])
Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

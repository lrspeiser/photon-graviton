"""Gaussian receiving overlap: color response and conditional angular spread."""
from pathlib import Path
import json,math
import numpy as np
from scipy.special import gammainc,gamma
from scipy.integrate import quad
from scipy.constants import hbar,c,electron_volt
P=Path(__file__).resolve().parent
gap=1e-8

def angular(E,a):
 outgoing=E-gap;x=2*a*a*E*outgoing
 if x<.1:
  I=quad(lambda u:(2-2*u+u*u)*math.exp(-x*u),0,2,epsabs=1e-13)[0]
  J=quad(lambda u:u*(2-2*u+u*u)*math.exp(-x*u),0,2,epsabs=1e-13)[0]
 else:
  moments=[gamma(k+1)*gammainc(k+1,2*x)/x**(k+1) for k in range(4)]
  I=2*moments[0]-2*moments[1]+moments[2]
  J=2*moments[1]-2*moments[2]+moments[3]
 # Independent integral with the boundary layer rescaled to finite width.
 if x>=.1:
  upper=min(2*x,80.)
  independent=quad(lambda y:(2-2*y/x+(y/x)**2)*math.exp(-y),0,upper,epsabs=1e-12)[0]/x
  assert abs(independent/I-1)<1e-10
 return math.exp(-(a*gap)**2)*I/(8/3),J/I

def rate(E,a):
 A,u=angular(E,a)
 return gap*(E-gap)**3*A,u,A

if __name__=='__main__':
 out=dict(scope='Hypothesized Gaussian spatial overlap in prior heavy-store forward scalar interaction; no absolute coupling or full path simulation',gap_eV=gap,rows=[],optical_angular_examples=[])
 for a in [0.,1.,1e3,1e6,1e8]:
  ref=rate(1.,a)[0]
  for E in [.01,.1,1.,2.,10.,100.]:
   loss,u,A=rate(E,a);h=1e-4
   slope=(math.log(rate(E*math.exp(h),a)[0])-math.log(rate(E*math.exp(-h),a)[0]))/(2*h)
   out['rows'].append(dict(a_inverse_eV=a,spatial_sigma_m=a*hbar*c/electron_volt,photon_energy_eV=E,relative_fractional_loss=loss/ref,logarithmic_slope=slope,angular_suppression=A,mean_one_minus_cos=u))
  # Small-optical-depth, fixed-gap diagnostic only at a 2 eV normalization.
  optical=2.;alpha=.0002488993286382367;D=30.660139;tau=alpha*D;events=tau*optical/gap
  _,u,A=rate(optical,a);disp=2*events*u
  out['optical_angular_examples'].append(dict(a_inverse_eV=a,spatial_sigma_m=a*hbar*c/electron_volt,mean_events_leading_order=events,orientation_decorrelation_exponent=events*u,small_angle_approximation_applicable=bool(disp<.01),conditional_rms_arcsec=(math.sqrt(disp)*180/math.pi*3600 if disp<.01 else None),relative_energy_width_leading_order=math.sqrt(events)*gap/optical,angular_suppression=A))
 assert abs(angular(2.,0)[0]-1)<1e-12 and abs(angular(2.,0)[1]-1)<1e-12
 (P/'spatial-response-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps(out['optical_angular_examples'],indent=2))

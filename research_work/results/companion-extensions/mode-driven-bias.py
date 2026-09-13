"""Occupation-driven bias of a reciprocal two-mode exchange interaction."""
from pathlib import Path
import math,json
import numpy as np
P=Path(__file__).resolve().parent
kb=8.617333262145e-5
out=dict(scope='Selected equal-weight mode pairs, flat ladder matrix elements, weak incoherent rate limit; no integrated astronomical spectrum',occupation_designs=[],thermal_checks=[],matrix_checks=[])
for nh in [1e-12,1e-6,1.,1000.]:
 for target in [.1,.9,.99]:
  nl=target*nh/(1+(1-target)*nh)
  r=nl*(1+nh)/(nh*(1+nl))
  assert abs(r/target-1)<1e-12
  # Equal mode counts, closed reservoirs: maximum net transfer before bias changes sign.
  fraction=(nh-nl)/(2*nh)
  middle_h=nh*(1-fraction);middle_l=nl+nh*fraction
  assert abs(middle_h-middle_l)<=1e-12*nh
  out['occupation_designs'].append(dict(high_occupation=nh,desired_reverse_forward=target,required_low_occupation=nl,high_photon_fraction_transferred_at_equal_occupation=fraction))
for T in [3000.,6000.,10000.]:
 for delta in [1e-8,2e-10]:
  for dilution in [1.,1e-9]:
   high=2.;low=high-delta
   nh=dilution/math.expm1(high/(kb*T));nl=dilution/math.expm1(low/(kb*T))
   r=nl*(1+nh)/(nh*(1+nl));assert r>1
   if dilution==1:assert abs(math.log(r)-delta/(kb*T))<1e-14
   out['thermal_checks'].append(dict(temperature_K=T,transfer_eV=delta,dilution=dilution,reverse_forward_minus_one=r-1))
# Exact number-state ladder operators independently check stimulated factors and hermiticity.
N=5;ann=np.diag(np.sqrt(np.arange(1,N)),1);raise_site=np.array([[0.,0.],[1.,0.]])
forward=np.kron(raise_site,np.kron(ann,ann.T));H=forward+forward.T
assert np.array_equal(H,H.T)
for nh in [1,2,3]:
 for nl in [0,1,2]:
  initial=(0*N+nh)*N+nl;final=(1*N+nh-1)*N+nl+1
  amplitude=H[final,initial]
  assert abs(amplitude**2-nh*(nl+1))<1e-12
  assert H[initial,final]==amplitude
  out['matrix_checks'].append(dict(high_number=nh,low_number=nl,squared_amplitude=float(amplitude**2)))
(P/'mode-driven-bias-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(out['occupation_designs'][1]);print(out['thermal_checks'][4:8])

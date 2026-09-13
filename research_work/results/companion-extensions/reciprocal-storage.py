"""Reciprocal inverse scattering for an ideally flat forward fractional rate."""
from pathlib import Path
import json,math
import numpy as np
P=Path(__file__).resolve().parent
out=dict(scope='Equal-degeneracy time-reversal-invariant heavy-store channels; ideal flat forward coefficient; inherited bosonic pair occupation model',rows=[])
for x in [.001,.01,.1,.5]:
 for loss in [0.,1e-6,.1,1.]:
  A=1.;B=(1+x)**3;C=loss;ratio=A/(B+C);mean=ratio/(1-ratio)
  cutoff=math.ceil(math.log(1e-15)/math.log(ratio));n=np.arange(cutoff+1);prob=(1-ratio)*ratio**n
  normalization=prob.sum();nummean=float(np.dot(n,prob));m2=float(np.dot(n*n,prob));plus=float(np.dot((n+1)**2,prob))
  birth=A*plus;inverse=B*m2;decay=C*m2
  assert abs(normalization-1)<2e-14 and abs(nummean/mean-1)<1e-10
  assert abs(birth-inverse-decay)/birth<1e-10
  out['rows'].append(dict(gap_over_photon_energy=x,inverse_over_forward=B,spontaneous_over_forward=C,mean_pairs=mean,stored_energy_over_incident_photon=x*mean,stationary_inverse_fraction_of_birth=inverse/birth,stationary_decay_fraction_of_birth=decay/birth,net_photon_loss_equals_decay_normalized_residual=(birth-inverse-decay)/birth,occupation_cutoff=cutoff))
out['analytic_small_gap_zero_decay_limit']=1/3
out['old_point_response_small_gap_limit']=1/6
(P/'reciprocal-storage-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([q for q in out['rows'] if q['spontaneous_over_forward']==0],indent=2))

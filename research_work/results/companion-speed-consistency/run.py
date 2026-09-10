from pathlib import Path
import json,hashlib
import numpy as np
from kernel import solve
HERE=Path(__file__).resolve().parent
prior=HERE.parent/'continuous-signal-feedback/results.json'
old=json.loads(prior.read_text())
runs=[];refinements=[]
for energy in [.01,.1]:
 for speed in [.5,.75,1.]:
  pair=[]
  for modes in [32,64]:
   x,_=solve(energy,32,modes=modes,wave_speed=speed);runs.append(x);pair.append(x)
  keys=['mean_carrier_stretch','energy_weighted_duration_ratio','field_energy_gain']
  change={key:abs(pair[0][key]-pair[1][key]) for key in keys}
  assert max(change.values())<1e-6
  refinements.append({'energy':energy,'wave_speed':speed,'mode_changes':change})
  if speed==.5:
   reference=next(x for x in old['runs'] if x['signal_energy']==energy and x['quadrature_nodes']==32)
   for key in keys:np.testing.assert_allclose(pair[0][key],reference[key],rtol=1e-10,atol=1e-12)
  print(energy,speed,pair[-1]['mean_carrier_stretch'],pair[-1]['sampled_companion_to_photon_speed_ratio'],flush=True)
result={'runs':runs,'refinements':refinements,'scope':'Fixed-speed companions in the existing optical Hamiltonian; no claim that fixed speed one is locally luminal throughout the evolved field.',
 'source_hashes':{str(p.relative_to(HERE.parents[2])):hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'run.py',HERE/'kernel.py',prior]},
 'known_speed_relation':'v_photon=c0/n_bar; v_companion=v_constant; ratio=v_constant*n_bar/c0, independent of a shared clock/ruler conversion.',
 'matched_speed_full_signal_tested':False}
(HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

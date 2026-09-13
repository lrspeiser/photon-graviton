"""Local special-relativistic kinematics; energies in units of M0*c^2."""
import json,math
from pathlib import Path
rows=[]
for E in [1e-9,.01,1.,10.]:
 for beta in [0.,.5,1.]:
  M=1.;Mf=math.sqrt(M*M+2*M*E+(1-beta*beta)*E*E)
  gain=E*(2*M+(1-beta*beta)*E)/(Mf+M)
  recoil=beta*beta*E*E/(M+E+Mf)
  assert abs((gain+recoil)/E-1)<1e-12
  rows.append(dict(incident_energy_over_initial_rest_energy=E,directional_imbalance=beta,stored_rest_energy_fraction=gain/E,recoil_energy_fraction=recoil/E,final_speed_over_c=beta*E/(M+E)))
Path(__file__).with_name('capture-ledger-results.json').write_text(json.dumps({'scope':'Absorption kinematics only; no interaction probability or cooling law','rows':rows},indent=2)+'\n',encoding='utf-8',newline='\n')

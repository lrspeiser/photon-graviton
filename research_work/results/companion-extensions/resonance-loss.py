"""Competing photon removal and small-gap transfer, with a full drift ledger."""
from pathlib import Path
import json,math
from scipy.integrate import solve_ivp
P=Path(__file__).resolve().parent

def ledger(E0,gap,z,r):
 tau=math.log1p(z);Ef=E0/(1+z);lost=E0*z/(1+z);optical=r*lost/gap
 T=math.exp(-optical)
 C=(gap/r*(-math.expm1(-optical)) if r else lost)
 beam=T*Ef;removed=E0-beam-C
 return dict(transmission=T,transmitted_energy_eV=beam,companion_transfer_eV=C,removed_sector_energy_eV=removed,removal_optical_depth=optical)

if __name__=='__main__':
 A=.0002488993286382367*30.660139
 out=dict(scope='Constant removal-to-inelastic rate ratio, deterministic small-gap energy drift; no inferred microscopic branching ratio',cases=[],checks=[])
 for E in [.01,2.,100.]:
  for gap in [1e-8,1e-6,1e-4]:
   for z in [math.expm1(A),.1,1.]:
    rmax=-math.log(.9)*gap/(E*z/(1+z))
    rows=[]
    for factor in [0.,.1,1.,10.]:
     r=factor*rmax;v=ledger(E,gap,z,r);v.update(removal_ratio=r,multiple_of_90_percent_limit=factor);rows.append(v)
     assert abs(sum(v[k] for k in ['transmitted_energy_eV','companion_transfer_eV','removed_sector_energy_eV'])/E-1)<1e-12
    out['cases'].append(dict(initial_energy_eV=E,gap_eV=gap,z=z,ratio_limit_for_90_percent_transmission=rmax,rows=rows))
 # Independent integration for a representative energy, all three paths and ratios.
 for z in [math.expm1(A),.1,1.]:
  E=2.;gap=1e-8;rmax=-math.log(.9)*gap/(E*z/(1+z))
  for factor in [.1,1.,10.]:
   r=rmax*factor
   def rhs(t,y):
    photon=E*math.exp(-t);removal=r*photon/gap
    return [-removal*y[0],y[0]*photon,y[0]*removal*photon]
   sol=solve_ivp(rhs,[0,math.log1p(z)],[1.,0.,0.],rtol=1e-11,atol=1e-13)
   exact=ledger(E,gap,z,r);T,C,H=sol.y[:,-1]
   error=max(abs(T-exact['transmission']),abs(C-exact['companion_transfer_eV'])/E,abs(H-exact['removed_sector_energy_eV'])/E)
   assert sol.success and error<1e-9
   out['checks'].append(dict(z=z,ratio=r,maximum_normalized_ODE_error=error))
 (P/'resonance-loss-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps([q for q in out['cases'] if q['initial_energy_eV']==2 and q['gap_eV']==1e-8],indent=2))

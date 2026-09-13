"""Conditional aggregate release scaling and reversible assembly first passage."""
from pathlib import Path
import math,json
import numpy as np
P=Path(__file__).resolve().parent
old=json.loads((P/'site-mode-budget-results.json').read_text())
def logsuccess(K,r):
 if r==1:return -math.log(K)
 x=math.log(r)
 if r>1:
  return math.log(math.expm1(x))-K*x-math.log1p(-math.exp(-K*x))
 return math.log(-math.expm1(x))-math.log(-math.expm1(K*x))
ratios=[.9,1.,1.000001,1.01**3]
out=dict(scope='Conditional enlarged excitation and single collective relaxation quantum; assembly is a separate constant-rate ladder hypothesis',rows=[],checks=[])
for row in old['rows']:
 if row['radius_kpc']!=120 or row['target_protected_fraction']!=.9:continue
 for benchmark in [1e9,1e12,1e15]:
  L=row['protected_lifetime_lower_bound_years'];K=math.ceil((L/benchmark)**.25)
  assert L/K**4<=benchmark and L/(K-1)**4>benchmark
  out['rows'].append(dict(baryons=row['baryons'],small_transfer_eV=row['bright_gap_eV'],log_bandwidth=row['log_bandwidth'],lifetime_benchmark_years=benchmark,required_contributions=K,aggregate_bright_energy_eV=K*row['bright_gap_eV'],collective_release_eV=.5*K*row['bright_gap_eV'],conditional_lifetime_bound_years=L/K**4,assembly=[dict(reverse_forward_ratio=r,log10_success_per_seed=logsuccess(K,r)/math.log(10)) for r in ratios]))
for K in [5,20,100]:
 for r in ratios:
  up=1/(1+r);down=r/(1+r)
  mat=np.eye(K-1)-np.diag(np.full(K-2,up),1)-np.diag(np.full(K-2,down),-1)
  rhs=np.zeros(K-1);rhs[-1]=up
  sol=np.linalg.solve(mat,rhs)
  exact=math.exp(logsuccess(K,r));error=abs(sol[0]/exact-1)
  assert error<1e-10
  out['checks'].append(dict(K=K,r=r,relative_error=error))
# Effective aggregate energy is independent of input quantum size except rounding.
for row in out['rows']:
 match=next(r for r in out['rows'] if r['baryons']==row['baryons'] and r['log_bandwidth']==row['log_bandwidth'] and r['lifetime_benchmark_years']==row['lifetime_benchmark_years'] and r['small_transfer_eV']==1e-8)
 assert abs(row['aggregate_bright_energy_eV']-match['aggregate_bright_energy_eV'])<=1.01e-8
(P/'collective-assembly-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
for row in out['rows']:
 if row['baryons']=='I' and row['lifetime_benchmark_years']==1e12:print(row)

"""Restarted assembly mean time and gross versus net energy traffic."""
from pathlib import Path
import math,json
import numpy as np
P=Path(__file__).resolve().parent
YEAR=365.25*86400;MSUN=1.98847e30;c=299792458.
old=json.loads((P/'collective-assembly-results.json').read_text())
inventory=json.loads((P/'coupled-torque-results.json').read_text())
def logtime(K,r):
 if r==1:return math.log(K*(K+1)/2)
 x=math.log(r)
 if r>1:
  y=(K+1)*x
  if y>50:lnum=y+math.log1p(-((K+1)*r-K)*math.exp(-y))
  else:lnum=math.log(math.expm1(y)-(K+1)*math.expm1(x))
 else:lnum=math.log(K*(1-r)-r+math.exp((K+1)*x))
 return lnum-2*math.log(abs(1-r))
out=dict(scope='Constant-rate ladder, restart at zero, ideal recovery of downward energy, immediate terminal protection; mean first-fill benchmark only',rows=[],checks=[])
for row in old['rows']:
 mass=next(r['mass_Msun'] for r in inventory['rows'] if r['baryons']==row['baryons'] and r['receiver_region_kpc']==[30,60])
 K=row['required_contributions']
 # Formation duration is an independent illustrative choice, not the lifetime benchmark.
 duration=1e12*YEAR
 netpower=2*mass*MSUN*c*c/duration
 for item in row['assembly']:
  r=item['reverse_forward_ratio'];lt=logtime(K,r);lover=lt-math.log(K)
  out['rows'].append(dict(baryons=row['baryons'],small_transfer_eV=row['small_transfer_eV'],log_bandwidth=row['log_bandwidth'],lifetime_design_years=row['lifetime_benchmark_years'],K=K,reverse_forward_ratio=r,mean_fill_benchmark_years=1e12,log10_dimensionless_mean_time=lt/math.log(10),log10_required_forward_rate_per_second=(lt-math.log(duration))/math.log(10),log10_gross_net_energy_ratio=lover/math.log(10),net_assembly_power_W=netpower,log10_gross_forward_power_W=math.log10(netpower)+lover/math.log(10),log10_recycling_loss_fraction_for_extra_input_below_net=-math.log(math.expm1(lover))/math.log(10) if lover<700 else -lover/math.log(10)))
for K in [5,20,100]:
 for r in [.9,1.,1.000001,1.01**3]:
  # Backward mean hitting-time equation, zero reflecting, K absorbing.
  mat=np.zeros((K,K));mat[0,0]=1
  if K>1:mat[0,1]=-1
  for j in range(1,K):
   mat[j,j]=1+r;mat[j,j-1]=-r
   if j+1<K:mat[j,j+1]=-1
  value=np.linalg.solve(mat,np.ones(K))[0]
  exact=math.exp(logtime(K,r));error=abs(value/exact-1)
  assert error<1e-8
  # Finite positive sum independently checks cancellation-sensitive expressions.
  summed=sum((K-j)*r**j for j in range(K))
  assert abs(summed/exact-1)<1e-8
  out['checks'].append(dict(K=K,r=r,relative_error=error))
(P/'assembly-recycling-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
for row in out['rows']:
 if row['baryons']=='I' and row['small_transfer_eV']==1e-8 and row['log_bandwidth']==.01 and row['lifetime_design_years']==1e12:print(row)

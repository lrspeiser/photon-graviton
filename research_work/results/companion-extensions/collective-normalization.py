"""Normalized collective excitation and forward/reverse Dicke factors."""
from pathlib import Path
import math,json
import numpy as np
P=Path(__file__).resolve().parent
out=dict(scope='Ideal phase-matched symmetric two-level domains, no dephasing, no protection sink; matrix normalization audit',matrix_checks=[],states=[],stationary=[])
for N in [2,4,6,8]:
 dim=2**N;S=np.zeros((dim,dim))
 for state in range(dim):
  for bit in range(N):
   if not state&(1<<bit):S[state|(1<<bit),state]=1
 for m in range(N+1):
  vec=np.array([1 if j.bit_count()==m else 0 for j in range(dim)],dtype=float)/math.sqrt(math.comb(N,m))
  up=float(np.dot(S@vec,S@vec));down=float(np.dot(S.T@vec,S.T@vec))
  assert abs(up-(N-m)*(m+1))<1e-10 and abs(down-m*(N-m+1))<1e-10
  out['matrix_checks'].append(dict(N=N,m=m,up_factor=up,down_factor=down))
for N in [10,100,10000]:
 for m in sorted(set([0,N//2,N-1,N])):
  out['states'].append(dict(N=N,m=m,up_factor=(N-m)*(m+1),down_factor=m*(N-m+1),up_per_constituent=(N-m)*(m+1)/N,preparation_energy_eV=m*1e-8))
 for r in [.9,1.,1.01**3]:
  m=np.arange(N+1,dtype=float);logp=-m*math.log(r);p=np.exp(logp-logp.max());p/=p.sum()
  up=(N-m)*(m+1);down=m*(N-m+1)
  meanup=float(p@up);meandown=float(p@down)
  assert abs(meanup-r*meandown)<1e-8*max(1,meanup)
  assert np.max(abs(p[:-1]*up[:-1]-r*p[1:]*down[1:]))<1e-8*max(1,meanup)
  out['stationary'].append(dict(N=N,r=r,mean_excitation=float(p@m),mean_up_per_constituent=meanup/N,net_current_over_forward_coefficient=meanup-r*meandown))
# A normalized single bright state has amplitude sqrt(N), not N.
for N in [10,100,10000]:
 bright=np.ones(N)/math.sqrt(N);created=np.ones(N)
 assert abs(abs(np.vdot(bright,created))**2/N-1)<1e-12
(P/'collective-normalization-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(out['stationary'][-3:])

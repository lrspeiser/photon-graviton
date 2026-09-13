"""Self-consistent steady source/outlet populations for a recycled loading ladder."""
from pathlib import Path
import math,json
import numpy as np
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
def logF(K,r):
 if r==0:return math.log(K)
 u=r-1
 if abs(K*u)<1e-4:
  return math.log(K*(K+1)/2)+math.log1p((K-1)*u/3+(K-1)*(K-2)*u*u/12)
 if r<1:num=K*(1-r)-r+math.exp((K+1)*math.log(r));return math.log(num)-2*math.log(1-r)
 y=(K+1)*math.log(r)
 numlog=y+math.log1p(-((K+1)*r-K)*math.exp(-y)) if y>50 else math.log(math.expm1(y)-(K+1)*u)
 return numlog-2*math.log(u)
def solve(K,ns,A,B):
 def current(j):
  nh=ns*max(0,1-A*j);nl=ns*B*j
  if nh==0:return 0.
  a=nh*(1+nl);b=nl*(1+nh);r=b/a
  return math.exp(math.log(K)+math.log(a/ns)-logF(K,r))
 j=brentq(lambda j:current(j)-j,0,1/A,xtol=1e-14)
 nh=ns*(1-A*j);nl=ns*B*j;a=nh*(1+nl);b=nl*(1+nh);r=b/a
 assert abs(current(j)-j)<1e-8*max(j,1e-10)
 result=dict(K=K,source_occupation=ns,source_load=A,outlet_load=B,high_occupation=nh,low_occupation=nl,reverse_forward_ratio=r,net_rate_over_Gamma_source_occupation=j,source_depletion_fraction=A*j,packet_rate_over_Gamma=ns*j/K,interior_approx_net_rate=1/(1+A+B))
 if K==20:
  Q=np.zeros((K,K))
  for k in range(K):
   Q[(k+1)%K,k]+=a;Q[k,k]-=a
   if k>0:Q[k-1,k]+=b;Q[k,k]-=b
  mat=Q.copy();rhs=np.zeros(K);mat[-1,:]=1;rhs[-1]=1
  prob=np.linalg.solve(mat,rhs)
  exact=K*a*prob[-1]/ns
  assert prob.min()>-1e-10 and abs(exact/j-1)<1e-8
  result['matrix_relative_error']=abs(exact/j-1)
 # E_H=2 eV, E_L=2-delta: terminal half/half output.
 delta=1e-8;rate=ns*j;incoming=2*rate;outgoing=(2-delta)*rate
 protected=.5*delta*rate;release=protected
 assert abs(incoming-outgoing-protected-release)<1e-14*incoming
 return result
out=dict(scope='Steady externally maintained high field, low escape, instantaneous terminal packet export/reset; terminal sink is not derived',rows=[])
for K in [20,399951]:
 for ns in [1e-12,1.,1000.]:
  for A in [.1,10.]:
   for B in [.1,10.,1000.]:out['rows'].append(solve(K,ns,A,B))
(P/'source-outlet-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
for r in out['rows']:
 if r['K']==399951 and r['source_occupation']==1e-12:print(r)

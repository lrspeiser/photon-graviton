"""Reciprocal collective transfer into one finite protected level."""
from pathlib import Path
import json,math
import numpy as np
P=Path(__file__).resolve().parent

def solve(N,r,nb,gamma):
 kappa=.01;size=2*(N+1);Q=np.zeros((size,size));pin=np.zeros(size);J=np.zeros(size);pb=np.zeros(size);mvals=np.zeros(size)
 def add(i,j,rate):Q[j,i]+=rate;Q[i,i]-=rate
 for m in range(N+1):
  up=(N-m)*(m+1);down=m*(N-m+1)
  for p in [0,1]:
   i=2*m+p;pb[i]=p;mvals[i]=m;pin[i]=up-r*down
   if m<N:add(i,2*(m+1)+p,up)
   if m>0:add(i,2*(m-1)+p,r*down)
   if p==0 and m>0:
    rate=kappa*(nb+1)*down;add(i,2*(m-1)+1,rate);J[i]=rate
   if p==1 and m<N:
    rate=kappa*nb*up;add(i,2*(m+1),rate);J[i]=-rate
   if p==1:add(i,2*m,gamma)
 mat=Q.copy();rhs=np.zeros(size);mat[-1,:]=1;rhs[-1]=1
 prob=np.linalg.solve(mat,rhs)
 assert prob.min()>-1e-9 and abs(prob.sum()-1)<1e-10
 residual=float(np.max(abs(Q@prob)));assert residual<1e-8
 protected=float(prob@pb);flow=float(prob@J);power=float(prob@pin)
 assert abs(flow-gamma*protected)<1e-8 and abs(power-flow)<1e-8
 # Energy per step=1, half relaxes; remaining half exits if protected level decays.
 assert abs(power-.5*flow-.5*gamma*protected)<1e-8
 exact=1 if nb==0 else (nb+1)/(nb+1+r*nb)
 if gamma==0:
  assert abs(protected-exact)<1e-8
  weights=np.exp(-np.arange(N+1)*math.log(r));weights/=weights.sum()
  expected=np.array([w*(1-exact) if p==0 else w*exact for w in weights for p in [0,1]])
  assert np.max(abs(prob-expected))<1e-8
 return dict(N=N,reverse_forward=r,bath_occupation=nb,protected_decay=gamma,protected_probability=protected,mean_collective_excitation=float(prob@mvals),net_transfer_rate=flow,net_input_power=power,relaxation_power=.5*flow,protected_decay_power=.5*gamma*protected,no_decay_analytic_probability=exact,stationary_residual=residual)
out=dict(scope='Finite symmetric ladder plus one protected two-state level; fixed pump and relaxation bath; no spatial transport or derived protection lifetime',rows=[])
for N in [10,100]:
 for r in [.9,1.,1.01**3]:
  for nb in [0.,.1,10.]:
   for gamma in [0.,.001,.1]:out['rows'].append(solve(N,r,nb,gamma))
(P/'collective-reservoir-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
for row in out['rows']:
 if row['reverse_forward']==.9 and row['bath_occupation'] in [0,10] and row['protected_decay'] in [0,.1]:print(row)

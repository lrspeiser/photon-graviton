"""Interior finite rays through a cube: tetrahedral field vs polar volume integral."""
from pathlib import Path
from itertools import permutations
import sys,json
import numpy as np
from scipy.integrate import quad_vec
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H.parent/'bar-volume-gravity'))
from tetra import gravity
cells=[]
for order in permutations(range(3)):
 v=[np.full(3,-1.)]
 for i in order:
  a=v[-1].copy();a[i]=1.;v.append(a)
 cells.append(v)
def reference(b,L,tol):
 corners=[np.arctan2(y-b[1],x-b[0])%(2*np.pi) for x in (-1,1) for y in (-1,1)]
 cuts=sorted([0.,2*np.pi]+corners)
 def primitive(r,a):return .5*(r*np.sqrt(r*r+a*a)+a*a*np.arcsinh(r/a))
 def angular(t):
  e=np.array([np.cos(t),np.sin(t)])
  bounds=[((1 if e[j]>0 else -1)-b[j])/e[j] for j in (0,1) if abs(e[j])>1e-15]
  r=min(bounds)
  return .5*e*(primitive(r,L+1)-primitive(r,L-1))
 return sum((quad_vec(angular,a,c,epsabs=tol,epsrel=tol)[0] for a,c in zip(cuts[:-1],cuts[1:])),np.zeros(2))
rows=[]
for b in ([.137,.219],[.8,.7],[.99,.97]):
 for L in (5.,20.):
  def f(z):return 2*gravity(cells,np.ones(6)/6,[[b[0],b[1],z]])[1][0,:2]
  calc,est=quad_vec(f,-L,L,points=[-1,1],epsabs=1e-9,epsrel=1e-9)
  r1=reference(b,L,1e-8);r2=reference(b,L,1e-11)
  err=float(np.linalg.norm(calc-r2)/np.linalg.norm(r2));refine=float(np.linalg.norm(r2-r1)/np.linalg.norm(r2))
  rows.append(dict(impact=b,L=L,ray_integral=calc.tolist(),polar_reference=r2.tolist(),relative_difference=err,reference_refinement=refine,integrator_estimate=float(est),passes=max(err,refine)<1e-7))
  print(b,L,err,refine,flush=True)
(H/'interior-results.json').write_text(json.dumps(dict(records=rows,all_pass=all(r['passes'] for r in rows)),indent=2)+'\n',encoding='utf8',newline='\n')
assert all(r['passes'] for r in rows)

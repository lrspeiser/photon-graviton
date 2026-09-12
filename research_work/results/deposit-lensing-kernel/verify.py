"""Conditional equal-potential lensing: tetrahedron field vs volume quadrature.
G=M=c=1 for numerical coefficients; finite ray segment, exterior projected rays.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import quad_vec
from scipy.special import roots_legendre
H=Path(__file__).resolve().parent
sys.path.insert(0,str(H.parent/'bar-volume-gravity'))
from tetra import gravity
v=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1.]])
rows=[]
for b in ([3.,2.],[1.5,.2]):
 for L in (5.,20.):
  def f(z):return 2*gravity([v],[1],[[b[0],b[1],z]])[1][0,:2]
  integrated,estimated=quad_vec(f,-L,L,epsabs=1e-10,epsrel=1e-10)
  refs=[]
  for n in (16,32):
   nodes,w=roots_legendre(n);nodes=(nodes+1)/2;w/=2
   u,t,q=np.meshgrid(nodes,nodes,nodes,indexing='ij');wu,wt,wq=np.meshgrid(w,w,w,indexing='ij')
   bary=np.stack([(1-u)*(1-t)*(1-q),u,(1-u)*t,(1-u)*(1-t)*q],axis=-1)
   x=bary@v;weights=6*wu*wt*wq*(1-u)**2*(1-t)
   assert abs(weights.sum()-1)<1e-12
   d=x[...,:2]-np.array(b);r2=np.sum(d*d,axis=-1)
   hi=L-x[...,2];lo=-L-x[...,2]
   end=hi/np.sqrt(r2+hi*hi)-lo/np.sqrt(r2+lo*lo)
   ref=2*np.sum(weights[...,None]*d/r2[...,None]*end[...,None],axis=(0,1,2))
   refs.append(ref)
  err=float(np.linalg.norm(integrated-refs[-1])/np.linalg.norm(refs[-1]))
  refine=float(np.linalg.norm(refs[-1]-refs[0])/np.linalg.norm(refs[-1]))
  rows.append(dict(impact=b,L=L,deflection_coefficient=integrated.tolist(),volume_reference=refs[-1].tolist(),relative_difference=err,volume_refinement=refine,integrator_estimate=float(estimated),passes=max(err,refine)<1e-7))
  assert rows[-1]['passes']
# Factor and sign check, independent analytic finite-segment point lens.
point=[]
for b in (.5,2.,10.):
 L=20.;numeric=quad_vec(lambda z:np.array([-2*b/(b*b+z*z)**1.5]),-L,L,epsabs=1e-12)[0][0]
 exact=-4*L/(b*np.sqrt(b*b+L*L))
 error=float(abs(numeric-exact)/abs(exact));assert error<1e-10
 point.append(dict(b=b,L=L,numeric=float(numeric),exact=float(exact),relative_error=error))
(H/'results.json').write_text(json.dumps(dict(tetrahedron=rows,point_lens=point,all_pass=True,scope='conditional static equal-potential weak-field exterior-ray benchmark; not a galaxy observation'),indent=2)+'\n',encoding='utf8',newline='\n')
print('Four tetrahedron ray checks:',max(r['relative_difference'] for r in rows),'reference refinement:',max(r['volume_refinement'] for r in rows))
print('Three finite point-lens checks:',max(r['relative_error'] for r in point))

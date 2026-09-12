"""New trajectory checks inside four previously unprobed symmetry classes per source.
Select largest mass-weighted final cell diameter. This is targeted diagnosis, not a holdout.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
sys.path[:0]=[str(H.parent/'bar-field-foundation'),str(H.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole
from axis import evaluate
from mesh import canonical,key
R=int(sys.argv[1]);assert R in (1,3)
M=json.loads((H/'mesh.json').read_text());V=np.array(M['vertices']);F=np.array(M['faces'])
L=np.array(M['lookup']);S=np.array(M['signs']);times=np.array([0,.005,.01,.025,.05,.1,.25])
D=json.loads((H/'local-refinement-diagnostic.json').read_text())[0 if R==1 else 1]
occupied={r['face'] for r in D['records'][0]['ranked_cells']}
rows=json.loads((H/f'prepared-R{R}.json').read_text())['records'];splines=[]
for row in rows:
 p=ROOT/row['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==row['cache_sha256']
 c=np.load(p);splines.append(CubicHermiteSpline(c['ages'],c['position'],c['velocity']))
C=np.stack([s(times) for s in splines],axis=1)[:,L]*S[None,:,:]
a,b,c=np.moveaxis(V[F],1,0)
omega=2*np.arctan2(abs(np.einsum('ij,ij->i',a,np.cross(b,c))),1+np.sum(a*b+b*c+c*a,axis=1))
cap=np.array([r['capture'] for r in rows])[L];mass=omega*np.mean(cap[F],axis=1);mass/=mass.sum()
diam=np.max(np.stack([np.linalg.norm(C[-1,F[:,i]]-C[-1,F[:,j]],axis=1) for i,j in ((0,1),(0,2),(1,2))]),axis=0)
chosen=[];seen=set()
for i in np.argsort(-mass*diam**2):
 if int(i) in occupied:continue
 center=V[F[i]].sum(axis=0);center/=np.linalg.norm(center)
 canon,sign=canonical(center);k=key(canon)
 if k in seen:continue
 seen.add(k);chosen.append((int(i),center,canon,sign))
 if len(chosen)==4:break
field=ROOT/'research_work/data-cache/bar-field/bar-L40.npz';bar=FastMultipole.load(field)
def cross(v):return np.array([-v[1],v[0],0.])
def rhs(t,y):
 _,a=evaluate(bar,y[:3]);return np.r_[y[3:]-37.5*cross(y[:3]),a[0]-37.5*cross(y[3:])]
result=[]
for face,center,direction,sign in chosen:
 sol=solve_ivp(rhs,(0,.25),np.r_[R*direction,-50*direction],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.0002,dense_output=True)
 assert sol.success
 y=sol.sol(times).T;actual=y[:,:3]*sign
 predicted=C[:,F[face]].mean(axis=1)
 error=np.linalg.norm(predicted-actual,axis=1)/R
 probe=sol.sol(np.linspace(0,.25,301)).T;pot=evaluate(bar,probe[:,:3])[0]
 J=.5*np.sum(probe[:,3:]**2,axis=1)+pot-37.5*(probe[:,0]*probe[:,4]-probe[:,1]*probe[:,3])
 drift=float(max(abs(J-J[0]))/220**2)
 result.append(dict(face=face,direction=center.tolist(),source_mass=float(mass[face]),
  times=times.tolist(),error_over_R=error.tolist(),actual_positions=actual.tolist(),
  predicted_positions=predicted.tolist(),jacobi_error=drift,orbit_check_passes=drift<1e-5))
 (H/f'unsampled-probes-R{R}.json').write_text(json.dumps(dict(R=R,expected=4,field_hash=hashlib.sha256(field.read_bytes()).hexdigest(),mesh_hash=hashlib.sha256((H/'mesh.json').read_bytes()).hexdigest(),selection='four unsampled symmetry classes ranked by source mass times final cell diameter squared',records=result),indent=2)+'\n',encoding='utf8',newline='\n')
 print(R,face,'final relative interpolation error',float(error[-1]),'orbit check',drift,flush=True)

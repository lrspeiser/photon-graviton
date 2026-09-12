"""Tighter-tolerance checks of selected new trajectories, including the axis."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole
from axis import evaluate
bar=FastMultipole.load(ROOT/'research_work/data-cache/bar-field/bar-L40.npz')
mesh=json.loads((HERE/'mesh.json').read_text(encoding='utf8'))
chosen=sorted(set([18,40,70,97]+[i for i,v in enumerate(mesh['directions']) if np.hypot(v[0],v[1])==0]))
def cross(v):return np.array([-v[1],v[0],0.])
def rhs(t,y):
    _,a=evaluate(bar,y[:3])
    return np.r_[y[3:]-37.5*cross(y[:3]),a[0]-37.5*cross(y[3:])]
outpath=HERE/'nested-orbit-checks.json'
rows=[]
for R in (1.,3.):
    d=json.loads((HERE/f'prepared-R{R:g}.json').read_text(encoding='utf8'));assert len(d['records'])==d['expected']==98
    for i in chosen:
        r=d['records'][i];direction=np.array(r['direction'])
        p=ROOT/r['cache'];assert hashlib.sha256(p.read_bytes()).hexdigest()==r['cache_sha256']
        c=np.load(p);s=CubicHermiteSpline(c['ages'],c['position'],c['velocity'])
        try:
            sol=solve_ivp(rhs,(0,.25),np.r_[R*direction,-50*direction],method='DOP853',rtol=2e-13,atol=2e-15,max_step=.0002,dense_output=True)
        except ValueError as error:
            rows.append(dict(R=R,index=i,passes=False,error=str(error),status='outside_original_field_domain'))
            outpath.write_text(json.dumps(dict(expected=2*len(chosen),records=rows),indent=2)+'\n',encoding='utf8',newline='\n')
            print(rows[-1],flush=True)
            continue
        assert sol.success
        t=np.linspace(0,.25,301);y=sol.sol(t).T
        velocity=y[:,3:]-37.5*np.c_[-y[:,1],y[:,0],np.zeros(len(y))]
        dx=float(np.max(np.linalg.norm(y[:,:3]-s(t),axis=1)))
        dv=float(np.max(np.linalg.norm(velocity-s(t,1),axis=1)))
        pot=evaluate(bar,y[:,:3])[0];J=.5*np.sum(y[:,3:]**2,axis=1)+pot-37.5*(y[:,0]*y[:,4]-y[:,1]*y[:,3])
        error=float(max(abs(J-J[0]))/220**2)
        rows.append(dict(R=R,index=i,position_change=dx,velocity_change=dv,invariant_error=error,
                         passes=dx<1e-5 and dv<.05 and error<1e-5))
        (HERE/'nested-orbit-checks.json').write_text(json.dumps(dict(expected=2*len(chosen),records=rows),indent=2)+'\n',encoding='utf8',newline='\n')
        print(rows[-1],flush=True)

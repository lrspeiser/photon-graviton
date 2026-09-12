"""Compute genuine new midpoint trajectories; reuse matching cached directions."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole
from axis import evaluate
R=float(sys.argv[1]);assert R==3.
mesh=json.loads((HERE/'mesh-adaptive.json').read_text(encoding='utf8'))
previous=json.loads((HERE/f'prepared2-R{R:g}.json').read_text(encoding='utf8'))
previous=[r for r in previous['records'] if r['R']==R]
path=ROOT/'research_work/data-cache/bar-field/bar-L40.npz';bar=FastMultipole.load(path)
field_hash=hashlib.sha256(path.read_bytes()).hexdigest()
cache=ROOT/'research_work/data-cache/bar-volumes';cache.mkdir(exist_ok=True,parents=True)
ages=np.unique(np.r_[np.linspace(0,.25,4097),.05,.1])
def cross(v):return np.array([-v[1],v[0],0.])
def rhs(t,y):
    _,a=evaluate(bar,y[:3])
    return np.r_[y[3:]-37.5*cross(y[:3]),a[0]-37.5*cross(y[3:])]
records=[]
for i,direction in enumerate(np.array(mesh['directions'])):
    distances=np.array([np.linalg.norm(direction-r['direction']) for r in previous])
    closest=int(np.argmin(distances));reused=distances[closest]<1e-12
    if reused:
        old=previous[closest];p=ROOT/old['cache']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==old['cache_sha256']
        row=dict(old,index=i,reused=True,direction_match_error=float(distances[closest]))
    else:
        sol=solve_ivp(rhs,(0,.25),np.r_[R*direction,-50*direction],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.0002,dense_output=True)
        assert sol.success
        y=sol.sol(ages).T;velocity=y[:,3:]-37.5*np.c_[-y[:,1],y[:,0],np.zeros(len(y))]
        probe=sol.sol(np.linspace(0,.25,301)).T;pot=evaluate(bar,probe[:,:3])[0]
        J=.5*np.sum(probe[:,3:]**2,axis=1)+pot-37.5*(probe[:,0]*probe[:,4]-probe[:,1]*probe[:,3])
        error=float(max(abs(J-J[0]))/220**2);W=-float(pot[0])
        target=cache/f'adaptive-R{R:g}-orbit{i:03d}.npz'
        np.savez_compressed(target,ages=ages,position=y[:,:3],velocity=velocity)
        row=dict(R=R,index=i,direction=direction.tolist(),capture=(W/(W+40000))**6,
                 invariant_error=error,passes=error<1e-5,reused=False,cache=str(target.relative_to(ROOT)),
                 cache_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    records.append(row)
    (HERE/f'prepared-adaptive-R{R:g}.json').write_text(json.dumps(dict(R=R,expected=len(mesh['directions']),field_hash=field_hash,
        mesh_hash=hashlib.sha256((HERE/'mesh-adaptive.json').read_bytes()).hexdigest(),records=records),indent=2)+'\n',encoding='utf8',newline='\n')
    print(R,i+1,len(mesh['directions']),row['reused'],row['passes'],flush=True)

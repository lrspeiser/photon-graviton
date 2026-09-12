"""Cache resolved bar trajectories for reusable mass-volume reconstructions."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole
n=int(sys.argv[1]);assert n in (6,8)
source=HERE.parent/('bar-source-persistence/audit-fine-L40.json' if n==6 else 'bar-angular-refinement/grid8-L40.json')
original=json.loads(source.read_text(encoding='utf8'))
path=ROOT/'research_work/data-cache/bar-field/bar-L40.npz'
source_hash=hashlib.sha256(path.read_bytes()).hexdigest();assert source_hash==original['source_hash']
bar=FastMultipole.load(path)
cache=ROOT/'research_work/data-cache/bar-volumes';cache.mkdir(exist_ok=True,parents=True)
ages=np.unique(np.r_[np.linspace(0,.25,4097),.05,.1])
def cross(v):return np.array([-v[1],v[0],0.])
def rhs(t,y):
    _,a=bar.evaluate(y[:3])
    return np.r_[y[3:]-37.5*cross(y[:3]),a[0]-37.5*cross(y[3:])]
records=[]
for i,r in enumerate(original['records']):
    R,u,phi=r['R'],r['mu'],r['phi']
    direction=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
    sol=solve_ivp(rhs,(0,.25),np.r_[R*direction,-50*direction],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.0002,dense_output=True)
    assert sol.success
    y=sol.sol(ages).T;velocity=y[:,3:]-37.5*np.c_[-y[:,1],y[:,0],np.zeros(len(y))]
    probe=sol.sol(np.linspace(0,.25,301)).T;pot=bar.evaluate(probe[:,:3])[0]
    J=.5*np.sum(probe[:,3:]**2,axis=1)+pot-37.5*(probe[:,0]*probe[:,4]-probe[:,1]*probe[:,3])
    error=float(max(abs(J-J[0]))/220**2)
    target=cache/f'source{n}-orbit{i:03d}.npz'
    np.savez_compressed(target,ages=ages,position=y[:,:3],velocity=velocity)
    W=-float(pot[0]);row={k:r[k] for k in ('R','mu','phi','raw_weight')}
    row.update(direction=direction.tolist(),capture=(W/(W+40000))**6,invariant_error=error,passes=error<1e-5,
               cache=str(target.relative_to(ROOT)),cache_sha256=hashlib.sha256(target.read_bytes()).hexdigest())
    records.append(row)
    (HERE/f'prepared{n}.json').write_text(json.dumps(dict(n=n,expected=n*n,source_hash=source_hash,records=records),indent=2)+'\n',encoding='utf8',newline='\n')
    print(n,i+1,len(original['records']),row['passes'],flush=True)

"""Whole-trajectory audit with turning-point-bracketed core crossings."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole

from kernel import integrate_kernel,compare,POINTS
grid='fine';order=40
launch_R=float(sys.argv[1]);assert launch_R in (1.,3.)
original=json.loads((HERE.parent/'bar-angular-refinement/grid12-L40.json').read_text())
path=ROOT/'research_work/data-cache/bar-field/bar-L40.npz'
assert hashlib.sha256(path.read_bytes()).hexdigest()==original['source_hash']
bar=FastMultipole.load(path)

def cross_z(v):return np.array([-v[1],v[0],0.])

def integrate(y0,tol):
    def rhs(t,y):
        _,a=bar.evaluate(y[:3])
        return np.r_[y[3:]-37.5*cross_z(y[:3]),a[0]-37.5*cross_z(y[3:])]
    def turn(t,y):return y[:3]@y[3:]/np.linalg.norm(y[:3])
    def boundary(t,y):return np.linalg.norm(y[:3])-.1
    sol=solve_ivp(rhs,(0,.25),y0,method='DOP853',rtol=tol,atol=tol*.01,max_step=.0002,
                  events=(turn,boundary),dense_output=True)
    assert sol.success
    brackets=np.r_[0.,sol.t_events[0],.25]
    crossings=[]
    def f(t):return np.linalg.norm(sol.sol(t)[:3])-.1
    for a,b in zip(brackets[:-1],brackets[1:]):
        if f(a)*f(b)<0:crossings.append(float(brentq(f,a,b,xtol=1e-13)))
    durations=[]
    for T in (.05,.1,.25):
        inside=False;previous=0.;duration=0.
        for t in crossings+[.25]:
            if inside:duration+=max(0.,min(t,T)-previous)
            previous=t;inside=not inside
        durations.append(duration)
    times=np.linspace(0,.25,301);y=sol.sol(times).T
    E=.5*np.sum(y[:,3:]**2,axis=1)+bar.evaluate(y[:,:3])[0]
    Lz=y[:,0]*y[:,4]-y[:,1]*y[:,3];J=E-37.5*Lz
    entry=crossings[0] if crossings else None
    turns_before=[float(t) for t in sol.t_events[0] if entry is not None and t<entry]
    kernels=[]
    for T in (.05,.1,.25):
        position=lambda t:sol.sol(t)[:3].T
        a=integrate_kernel(position,sol.t,T,nodes=4)
        b=integrate_kernel(position,sol.t,T,nodes=8)
        initial=compare(a,b);extra=None
        if not initial['passes']:extra=integrate_kernel(position,sol.t,T,nodes=16)
        final=compare(b,extra) if extra else initial
        kernels.append(dict(T=T,coarse=a,refined=b,extra_refined=extra,initial_check=initial,final_check=final))
    return dict(kernels=kernels,entry_time=entry,inside_durations=durations,crossings=crossings,
                naive_crossing_count=len(sol.t_events[1]),entry_after_turn=bool(turns_before),
                endpoint=sol.y[:,-1].tolist(),invariant_error=float(max(abs(J-J[0]))/220**2),evaluations=sol.nfev)

records=[]
launches=[r for r in original['records'] if r['R']==launch_R]
assert len(launches)==72
for i,r in enumerate(launches):
    R,u,phi=r['R'],r['mu'],r['phi']
    direction=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
    result=integrate(np.r_[R*direction,-50*direction],2e-11)
    row={k:r[k] for k in ('R','mu','phi','raw_weight')}
    row.update(result=result,passes=result['invariant_error']<1e-5 and all(k['final_check']['passes'] for k in result['kernels']))
    records.append(row)
    print(launch_R,i+1,72,row['passes'],flush=True)
    (HERE/f'orbits12-R{launch_R:g}.json').write_text(json.dumps(dict(records=records,expected=72,source_hash=original['source_hash'],
        points=POINTS.tolist()),indent=2)+'\n',encoding='utf8',newline='\n')

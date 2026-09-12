"""Whole-trajectory audit with turning-point-bracketed core crossings."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole

from scipy.special import roots_legendre
n=int(sys.argv[1]);assert n in (8,12)
order=40
path=ROOT/f'research_work/data-cache/bar-field/bar-L{order}.npz'
bar=FastMultipole.load(path)
source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
mu,weights=roots_legendre(n)
launches=[]
for R in (1.,3.):
    for u,w in zip(mu[mu>0],weights[mu>0]):
        for j in range(n):
            phi=(j+.5)*np.pi/n
            direction=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
            W=-float(bar.evaluate(R*direction)[0][0])
            launches.append(dict(R=R,mu=float(u),phi=float(phi),raw_weight=float(w/n*(W/(W+40000))**6)))

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
    return dict(entry_time=entry,inside_durations=durations,crossings=crossings,
                naive_crossing_count=len(sol.t_events[1]),entry_after_turn=bool(turns_before),
                endpoint=sol.y[:,-1].tolist(),invariant_error=float(max(abs(J-J[0]))/220**2),evaluations=sol.nfev)

def check(a,b):
    agreement=(a['entry_time'] is None)==(b['entry_time'] is None)
    dt=0. if a['entry_time'] is None or b['entry_time'] is None else abs(a['entry_time']-b['entry_time'])
    duration=max(abs(x-y) for x,y in zip(a['inside_durations'],b['inside_durations']))
    dv=float(np.linalg.norm(np.array(a['endpoint'][3:])-b['endpoint'][3:]))
    return dict(entry_agreement=agreement,entry_time_difference=dt,duration_difference=duration,velocity_difference=dv,
                passes=agreement and dt<1e-5 and duration<1e-5 and dv<.05 and b['invariant_error']<1e-5)

records=[]
for i,r in enumerate(launches):
    R,u,phi=r['R'],r['mu'],r['phi']
    direction=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
    y0=np.r_[R*direction,-50*direction]
    a=integrate(y0,2e-11)
    b=integrate(y0,2e-13) if i%8==0 or a['invariant_error']>=1e-5 else None
    comparison=check(a,b) if b else None
    row=dict(r,primary=a,refined=b,comparison=comparison,
             passes=(b or a)['invariant_error']<1e-5 and (comparison is None or comparison['passes']))
    records.append(row)
    print(n,i+1,len(launches),row['passes'],flush=True)
    (HERE/f'grid{n}-L40.json').write_text(json.dumps(dict(nmu=n,nphi=2*n,order=40,expected=len(launches),
        records=records,source_hash=source_hash),indent=2)+'\n',encoding='utf8',newline='\n')

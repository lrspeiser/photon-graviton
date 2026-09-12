"""Whole-trajectory audit with turning-point-bracketed core crossings."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole

grid=sys.argv[1];order=int(sys.argv[2])
original=json.loads((HERE/f'{grid}-L{order}.json').read_text())
assert len(original['records'])==(16 if grid=='coarse' else 36)
path=ROOT/f'research_work/data-cache/bar-field/bar-L{order}.npz'
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
    return dict(entry_time=entry,inside_durations=durations,crossings=crossings,
                naive_crossing_count=len(sol.t_events[1]),entry_after_turn=bool(turns_before),
                endpoint=sol.y[:,-1].tolist(),invariant_error=float(max(abs(J-J[0]))/220**2),evaluations=sol.nfev)

records=[]
for i,r in enumerate(original['records']):
    R,u,phi=r['R'],r['mu'],r['phi']
    direction=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
    y0=np.r_[R*direction,-50*direction]
    a=integrate(y0,2e-9);b=integrate(y0,2e-11)
    def check(a,b):
        agreement=(a['entry_time'] is None)==(b['entry_time'] is None)
        dt=0. if a['entry_time'] is None or b['entry_time'] is None else abs(a['entry_time']-b['entry_time'])
        duration=max(abs(x-y) for x,y in zip(a['inside_durations'],b['inside_durations']))
        dv=float(np.linalg.norm(np.array(a['endpoint'][3:])-b['endpoint'][3:]))
        return dict(entry_agreement=agreement,entry_time_difference=dt,duration_difference=duration,velocity_difference=dv,
                    passes=agreement and dt<1e-5 and duration<1e-5 and dv<.05 and b['invariant_error']<1e-5)
    initial=check(a,b);extra=None;final=initial
    if not initial['passes']:
        extra=integrate(y0,2e-13);final=check(b,extra)
    row={k:r[k] for k in ('R','mu','phi','raw_weight')}
    row.update(coarse=a,refined=b,extra_refined=extra,initial_check=initial,final_check=final)
    records.append(row)
    print(grid,order,i+1,(extra or b)['entry_time'],(extra or b)['inside_durations'],final['passes'],flush=True)
    (HERE/f'audit-{grid}-L{order}.json').write_text(json.dumps(dict(grid=grid,order=order,records=records,source_hash=original['source_hash']),indent=2)+'\n',encoding='utf8',newline='\n')

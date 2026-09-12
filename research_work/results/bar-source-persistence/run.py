from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.special import roots_legendre
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole

grid=sys.argv[1];order=int(sys.argv[2]);assert grid in ('coarse','fine') and order in (40,64)
nmu,nphi=(4,8) if grid=='coarse' else (6,12)
path=ROOT/f'research_work/data-cache/bar-field/bar-L{order}.npz'
bar=FastMultipole.load(path);omega=37.5
OUT=HERE/f'{grid}-L{order}.json'

def cross_z(v):return np.array([-v[1],v[0],0.])

def integrate(y0,tol):
    # Coordinates rotate with the bar; p is inertial velocity expressed in these axes.
    def rhs(t,y):
        _,a=bar.evaluate(y[:3])
        return np.r_[y[3:]-omega*cross_z(y[:3]),a[0]-omega*cross_z(y[3:])]
    def core(t,y):return np.linalg.norm(y[:3])-.1
    def turn(t,y):return y[:3]@y[3:]/np.linalg.norm(y[:3])
    core.terminal=True;core.direction=-1;turn.terminal=False;turn.direction=1
    sol=solve_ivp(rhs,(0,.25),y0,method='DOP853',rtol=tol,atol=tol*.01,max_step=.0002,events=(core,turn),dense_output=True)
    assert sol.success
    times=np.linspace(0,sol.t[-1],201);y=sol.sol(times).T
    pot=bar.evaluate(y[:,:3])[0]
    Lz=y[:,0]*y[:,4]-y[:,1]*y[:,3]
    E=.5*np.sum(y[:,3:]**2,axis=1)+pot;J=E-omega*Lz
    hit=float(sol.t_events[0][0]) if len(sol.t_events[0]) else None
    return dict(entry_time=hit,terminal_time=float(sol.t[-1]),endpoint=sol.y[:,-1].tolist(),
                turns_before_entry=sol.t_events[1].tolist(),entry_after_turn=hit is not None and len(sol.t_events[1])>0,
                invariant_error=float(max(abs(J-J[0]))/220**2),evaluations=sol.nfev)

mu,weights=roots_legendre(nmu)
records=[]
for R in (1.,3.):
    for u,w in zip(mu[mu>0],weights[mu>0]):
        for j in range(nphi//2):
            phi=(j+.5)*2*np.pi/nphi
            direction=np.array([np.sqrt(1-u*u)*np.cos(phi),np.sqrt(1-u*u)*np.sin(phi),u])
            pos=R*direction;W=-float(bar.evaluate(pos)[0][0])
            weight=float(w/(nphi//2)*(W/(W+40000))**6)
            y0=np.r_[pos,-50*direction]
            a=integrate(y0,2e-9);b=integrate(y0,2e-11)
            def check(a,b):
                event=(a['entry_time'] is None)==(b['entry_time'] is None)
                time=abs(a['terminal_time']-b['terminal_time'])
                dv=float(np.linalg.norm(np.array(a['endpoint'][3:])-b['endpoint'][3:]))
                return dict(event_agreement=event,time_difference=time,velocity_difference=dv,
                            passes=event and time<1e-5 and dv<.05 and b['invariant_error']<1e-5)
            initial=check(a,b);extra=None;final=initial
            if not initial['passes']:
                extra=integrate(y0,2e-13);final=check(b,extra)
            record=dict(R=R,mu=float(u),phi=float(phi),raw_weight=weight,coarse=a,refined=b,
                        extra_refined=extra,initial_check=initial,final_check=final)
            records.append(record)
            print(grid,order,len(records),R,round(u,3),round(phi,3),(extra or b)['entry_time'],final['passes'],flush=True)
            out=dict(grid=grid,order=order,nmu=nmu,nphi=nphi,records=records,
                     source_hash=hashlib.sha256(path.read_bytes()).hexdigest())
            OUT.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')

from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole

mode=sys.argv[1];assert mode in ('spherical','static','rotating')
omega=37.5 if mode=='rotating' else 0.
def rotate(x,angle):
    c=np.cos(angle);s=np.sin(angle)
    return np.array([c*x[0]-s*x[1],s*x[0]+c*x[1],x[2]])

def integrate(bar,y0,tol):
    monopole=CubicSpline(np.log(bar.r),bar.coeff[:,0]/np.sqrt(4*np.pi))
    def field(t,x):
        if mode=='spherical':
            r=np.linalg.norm(x);p=float(monopole(np.log(r)))
            return p,-float(monopole(np.log(r),1))*x/r**2
        p,a=bar.evaluate(rotate(x,-omega*t))
        return float(p[0]),rotate(a[0],omega*t)
    def rhs(t,y):return np.r_[y[3:],field(t,y[:3])[1]]
    def core(t,y):return np.linalg.norm(y[:3])-.1
    def turn(t,y):return y[:3]@y[3:]/np.linalg.norm(y[:3])
    core.terminal=True;core.direction=-1;turn.terminal=True;turn.direction=1
    sol=solve_ivp(rhs,(0,.05),y0,method='DOP853',rtol=tol,atol=tol*.01,max_step=.0002,
                  events=(core,turn),dense_output=True)
    assert sol.success
    event='core_entry' if len(sol.t_events[0]) else 'first_turn' if len(sol.t_events[1]) else 'time_limit'
    t=np.linspace(0,sol.t[-1],151);y=sol.sol(t).T
    angular=np.cross(y[:,:3],y[:,3:]);pot=np.array([field(tt,row[:3])[0] for tt,row in zip(t,y)])
    E=.5*np.sum(y[:,3:]**2,axis=1)+pot
    invariant=E-omega*angular[:,2]
    return dict(event=event,time=float(t[-1]),radius=float(np.linalg.norm(y[-1,:3])),
                endpoint=y[-1].tolist(),L_final=angular[-1].tolist(),max_L=float(np.max(np.linalg.norm(angular,axis=1))),
                delta_E=float(E[-1]-E[0]),rotation_work=float(omega*(angular[-1,2]-angular[0,2])),
                invariant_error_over_220_squared=float(max(abs(invariant-invariant[0]))/220**2),evaluations=sol.nfev)

records=[];hashes={}
for order in ((64,) if mode=='spherical' else (40,64)):
    path=ROOT/f'research_work/data-cache/bar-field/bar-L{order}.npz'
    hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    bar=FastMultipole.load(path)
    for R in (1.,3.):
        for phi in (np.pi/6,np.pi/3):
            for mu in (.2,.6):
                direction=np.array([np.sqrt(1-mu*mu)*np.cos(phi),np.sqrt(1-mu*mu)*np.sin(phi),mu])
                y0=np.r_[R*direction,-50*direction]
                a=integrate(bar,y0,2e-9);b=integrate(bar,y0,2e-11)
                check=dict(event_agreement=a['event']==b['event'],radius_change=abs(a['radius']-b['radius']),
                           velocity_change=float(np.linalg.norm(np.array(a['endpoint'][3:])-b['endpoint'][3:])))
                check['passes']=check['event_agreement'] and check['radius_change']<1e-4 and check['velocity_change']<.01 and b['invariant_error_over_220_squared']<1e-5
                records.append(dict(order=order,launch_radius=R,phi=float(phi),z_over_r=mu,coarse=a,refined=b,check=check))
                print(mode,order,R,round(phi,3),mu,b['event'],b['radius'],b['L_final'],check['passes'],flush=True)
                (HERE/f'{mode}.json').write_text(json.dumps(dict(mode=mode,omega=omega,records=records,hashes=hashes),indent=2)+'\n',encoding='utf8',newline='\n')

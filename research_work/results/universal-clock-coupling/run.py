import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import sympy as s

HERE=Path(__file__).resolve().parent


def profile(mode,t,x):
    if mode=='static_fast':sign,e,de=1,.2,0.
    elif mode=='static_slow':sign,e,de=-1,.2,0.
    elif mode=='relaxing_fast':sign,e,de=1,.2*math.exp(-.2*t),-.04*math.exp(-.2*t)
    elif mode=='growing_slow':sign,e,de=-1,.2-.1*math.exp(-.2*t),.02*math.exp(-.2*t)
    else:sign,e,de=-1,.2*math.exp(-.2*t),-.04*math.exp(-.2*t)
    I=math.sin(math.pi*x)**2 if 0<x<1 else 0.
    Ix=math.pi*math.sin(2*math.pi*x) if 0<x<1 else 0.
    return 1+sign*e*I,sign*e*Ix,sign*de*I,1+sign*e


def trace(mode,launch,tolerance):
    def ode(t,y):
        q,qx,qt,_=profile(mode,t,y[0])
        return [q,qt/q,qx]
    def arrived(t,y):return y[0]-2
    arrived.terminal=True
    sol=solve_ivp(ode,(launch,launch+8),[-1.,0.,0.],method='DOP853',rtol=tolerance,atol=tolerance*.01,max_step=(.01 if tolerance<1e-11 else .04),events=arrived)
    assert sol.success and len(sol.t_events[0])==1
    arrival=sol.t_events[0][0];y=sol.y_events[0][0]
    return arrival,math.exp(-y[1]),math.exp(y[2])


# Established Hamiltonian algebra; local c retained by the common clock factor.
p,m,c,q,qx=s.symbols('p m c q qx',positive=True)
H=q*s.sqrt(c*c*p*p+m*m*c**4)
v=s.diff(H,p)/q
proper_acceleration_at_rest=s.simplify(s.diff(v,p).subs(p,0)*(-s.diff(H,q)*qx/q).subs(p,0))
assert s.simplify(proper_acceleration_at_rest+c*c*qx/q)==0
assert s.simplify(s.diff(q*c*p,p)/q-c)==0

records=[]
for mode in ['static_fast','static_slow','relaxing_fast','growing_slow','relaxing_slow']:
    for launch in [0.,2.]:
        arrival,stretch,timing=trace(mode,launch,2e-12)
        arrival_lo,stretch_lo,timing_lo=trace(mode,launch,2e-10)
        h=1e-4
        fd=(trace(mode,launch+h,2e-12)[0]-trace(mode,launch-h,2e-12)[0])/(2*h)
        bound=profile(mode,launch,0)[3]/profile(mode,arrival,0)[3]
        assert abs(stretch-timing)<2e-9
        assert abs(fd-timing)<2e-6
        assert abs(stretch-stretch_lo)<2e-8
        if mode.startswith('static'):assert abs(stretch-1)<2e-9
        elif mode in ('relaxing_fast','growing_slow'):assert 1<stretch<=bound+1e-9
        else:assert stretch<1
        records.append(dict(mode=mode,launch_time=launch,arrival_time=arrival,
            wavelength_stretch=stretch,instantaneous_event_stretch=timing,
            finite_difference_event_stretch=fd,redshift=stretch-1,
            void_endpoint_contrast_ratio=bound,tolerance_stretch_difference=abs(stretch-stretch_lo)))
result={'status':'Prescribed-background kinematic check, no dynamical source or astronomical fit.',
    'proper_acceleration_at_rest':str(proper_acceleration_at_rest),
    'local_photon_speed':'c',
    'records':records}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))

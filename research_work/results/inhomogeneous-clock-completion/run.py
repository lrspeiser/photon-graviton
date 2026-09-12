from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp,quad
HERE=Path(__file__).resolve().parent
records=[]
def run(beta,a,b,uniform=False):
    def W(x):return 1. if uniform else np.sin(np.pi*x)**2
    def Wx(x):return 0. if uniform else np.pi*np.sin(2*np.pi*x)
    def n(x,t):return 1+beta*t*W(x)
    def ray(te,dense=False):
        sol=solve_ivp(lambda x,y:[n(x,y[0]),-beta*W(x),beta*W(x)*y[2]],[a,b],[te,0.,1.],rtol=2e-12,atol=1e-13,max_step=.02,dense_output=dense)
        assert sol.success;return sol
    te=1.;sol=ray(te,True);to,q,J=sol.y[:,-1];ne=n(a,te);no=n(b,to);S=J*ne/no
    exactJ=np.exp(beta*((b-a) if uniform else (b-a)/2-(np.sin(2*np.pi*b)-np.sin(2*np.pi*a))/(4*np.pi)))
    observed_frequency_ratio=np.exp(q)*no/ne
    dt=1e-5;early=ray(te-dt).y[0,-1];late=ray(te+dt).y[0,-1]
    proper_in=quad(lambda t:1/n(a,t),te-dt,te+dt,epsabs=1e-14)[0]
    proper_out=quad(lambda t:1/n(b,t),early,late,epsabs=1e-14)[0]
    measured=proper_out/proper_in
    grid=np.linspace(a,b,501);time=sol.sol(grid)[0]
    speed_error=max(abs((1/n(x,t))/(1/n(x,t))-1) for x,t in zip(grid,time))
    force=max(abs(beta*t*Wx(x)/n(x,t)) for x,t in zip(grid,time))
    assert abs(J/exactJ-1)<1e-9 and abs(S*observed_frequency_ratio-1)<1e-9 and abs(measured/S-1)<1e-7
    if uniform:assert abs(S-1)<1e-9
    return dict(beta=beta,start=a,end=b,profile='uniform_control' if uniform else 'periodic',emission_time=te,arrival_time=float(to),coordinate_event_stretch=float(J),endpoint_clock_factor=float(ne/no),measured_stretch=float(S),measured_redshift=float(S-1),measured_frequency_ratio=float(observed_frequency_ratio),proper_interval_stretch=float(measured),extra_coordinate_delay=float(to-te-(b-a)),relative_EM_GW_delay=0.,max_local_speed_identity_error=float(speed_error),max_holding_acceleration_in_c0_squared_per_length=float(force),jacobian_check_error=float(abs(J/exactJ-1)))
for beta in [0.,.02,.2]:
    for a,b in [(0.,1.),(0.,2.),(0.,5.),(.25,1.25),(.25,1.)]:records.append(run(beta,a,b))
    records.append(run(beta,0.,2.,True))
(HERE/'results.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf8',newline='\n')
for r in records:print(r['profile'],r['beta'],r['start'],r['end'],'z=',round(r['measured_redshift'],8),'force=',round(r['max_holding_acceleration_in_c0_squared_per_length'],5))

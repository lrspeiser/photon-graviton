from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp,quad

HERE=Path(__file__).resolve().parent
rows=[]
for beta in [1e-6,1e-5]:
    for f in [0.,.1,1.]:
        for gamma in [0.,1.]:
            def phi(x,t):return -beta*(1+f*t/20)/np.sqrt(1+x*x)
            def ray(te):
                def rhs(x,y):
                    n=np.exp(-(1+gamma)*phi(x,y[0]));pt=-beta*f/(20*np.sqrt(1+x*x))
                    return [n,-(1+gamma)*n*pt]
                sol=solve_ivp(rhs,[-10,10],[te,0.],rtol=1e-12,atol=1e-14,max_step=.025)
                assert sol.success;return sol.y[:,-1]
            to,logJ=ray(0.);endpoint=phi(10,to)-phi(-10,0);logS=logJ+endpoint;S=np.exp(logS)
            h=.001;early=ray(-h)[0];late=ray(h)[0]
            din=quad(lambda t:np.exp(phi(-10,t)),-h,h,epsabs=1e-14)[0]
            dout=quad(lambda t:np.exp(phi(10,t)),early,late,epsabs=1e-14)[0]
            err=abs(dout/din/S-1);assert err<1e-7
            approx=beta*f*((1+gamma)*np.arcsinh(10)/10-1/np.sqrt(101))
            rows.append(dict(beta=beta,f=f,gamma=gamma,arrival=float(to),log_coordinate_stretch=float(logJ),log_endpoint_factor=float(endpoint),measured_z=float(np.expm1(logS)),clock_error=float(err),first_order_logS=float(approx),first_order_error=float(abs(logS-approx)),shared_messenger_delay=0.))
gravity=[]
for beta in [1e-6,1e-5]:
    m=beta
    density=lambda r:3*m/(4*np.pi*(1+r*r)**2.5)
    enclosed=quad(lambda r:4*np.pi*r*r*density(r),0,1,epsabs=1e-18)[0]
    vc2=beta/(2*np.sqrt(2));assert abs(enclosed/vc2-1)<1e-9
    U=quad(lambda r:.5*4*np.pi*r*r*density(r)*(-m/np.sqrt(1+r*r)),0,np.inf,epsabs=1e-24)[0]
    Uexact=-3*np.pi*m*m/32;assert abs(U/Uexact-1)<1e-9
    for gamma in [0.,1.]:
        angle=quad(lambda s:(1+gamma)*beta/(s*s+2)**1.5,-np.inf,np.inf,epsabs=1e-18)[0]
        expected=(1+gamma)*beta;assert abs(angle/expected-1)<1e-9
        gravity.append(dict(beta=beta,gamma=gamma,vc_over_c=float(np.sqrt(vc2)),vc_km_s=float(299792.458*np.sqrt(vc2)),deflection_radians=float(angle),deflection_arcsec=float(angle*180/np.pi*3600),binding_energy=float(U),binding_energy_over_rest_energy=float(U/m)))
(HERE/'results.json').write_text(json.dumps(dict(rays=rows,gravity_snapshots=gravity,scope='Prescribed growing well and weak-field gravity/clock diagnostics; no formation or observational fit.'),indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(dict(rays=rows,gravity_snapshots=gravity),indent=2))

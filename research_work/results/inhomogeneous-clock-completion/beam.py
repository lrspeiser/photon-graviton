"""Paraxial beam distances in the retained plane-symmetric lapse profile."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
rows=[]
for z,te in [(0.,1.),(.1,0.),(.1,1.),(1.,0.),(1.,1.)]:
    beta=2*np.log1p(z)
    def lapse(x,t):return 1+beta*t*np.sin(np.pi*x)**2
    def nx(x,t):return beta*t*np.pi*np.sin(2*np.pi*x)
    def rhs(x,y):
        n=lapse(x,y[0]);return [n,-beta*np.sin(np.pi*x)**2,1/(n*np.exp(y[1]))]
    sol=solve_ivp(rhs,[0,1],[te,0.,0.],rtol=2e-11,atol=1e-12,max_step=.02)
    assert sol.success
    to,logomega,I=sol.y[:,-1];po=np.exp(logomega);S=1/po
    DG=I;DA=po*I;DL=S*DG
    def tilted(start,end,time,p,theta):
        py=p*np.sin(theta)
        def ode(x,y):
            t,px,transverse=y;n=lapse(x,t);p2=px*px+py*py
            return [n*np.sqrt(p2)/px,p2*nx(x,t)/(n*px),py/px]
        s=solve_ivp(ode,[start,end],[time,p*np.cos(theta),0.],rtol=2e-11,atol=1e-12,max_step=.01)
        assert s.success
        return abs(s.y[2,-1]/theta)
    forward=tilted(0,1,te,1.,1e-5);backward=tilted(1,0,to,po,1e-5)
    assert abs(forward/DG-1)<1e-7 and abs(backward/DA-1)<1e-7
    rows.append(dict(target_redshift=z,emission_time=te,beta=beta,arrival_time=float(to),
        measured_stretch=float(S),coordinate_distance=1.,
        emission_area_distance_DG=float(DG),angular_distance_DA=float(DA),luminosity_distance_DL=float(DL),
        flux_relative_to_euclidean_S_squared=float(1/DG**2),
        magnitude_change_relative_to_euclidean_S_squared=float(5*np.log10(DG)),
        forward_tilt_distance=float(forward),backward_tilt_distance=float(backward),
        reciprocity_residual=float(DL/(S*S*DA)-1)))
(HERE/'beam-results.json').write_text(json.dumps(dict(scope='Prescribed slab optics, no galaxy fit or field source',rows=rows),indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(rows,indent=2))

"""Independent radial energy quadrature for the spherical control."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
f=np.load(ROOT/'research_work/data-cache/bar-field/bar-L64.npz')
spline=CubicSpline(np.log(f['r']),f['coeff'][:,0]/np.sqrt(4*np.pi))
def phi(r):return spline(np.log(np.maximum(r,f['r'][0])))

results=[]
for R in (1.,3.):
    E=float(phi(R)+.5*50**2)
    apo=brentq(lambda r:float(phi(r))-E,R,100.)
    pair=[]
    for n in (128,256):
        z,w=roots_legendre(n)
        def travel(r):
            limit=np.arcsin(min(r/apo,1.))
            theta=(1+z)*limit/2;radius=apo*np.sin(theta)
            speed=np.sqrt(2*(E-phi(radius)))
            return float(np.sum(w*apo*np.cos(theta)/speed)*limit/2)
        center=travel(R);width=travel(.1);quarter=travel(apo)
        records=[]
        for T in (.05,.1,.25):
            dwell=0.
            for j in range(int((T+width)/(2*quarter))+2):
                t=center+2*j*quarter
                dwell+=max(0.,min(t+width,T)-max(t-width,0.))
            records.append(dict(T=T,central_residence_fraction=dwell/T,
                                previously_entered_fraction=max(T-(center-width),0.)/T))
        pair.append(dict(nodes=n,apocenter=apo,first_entry=center-width,center_to_apocenter_time=quarter,statistics=records))
    error=max(abs(a['central_residence_fraction']-b['central_residence_fraction']) for a,b in zip(pair[0]['statistics'],pair[1]['statistics']))
    assert error<1e-6
    results.append(dict(R=R,coarse=pair[0],refined=pair[1],fraction_change=error))
(HERE/'spherical-residence.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
print(results)

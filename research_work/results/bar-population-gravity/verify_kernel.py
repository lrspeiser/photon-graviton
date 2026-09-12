"""Independent analytic ring and time-weighted Kepler checks."""
from pathlib import Path
import json
import numpy as np
from kernel import integrate_kernel

T=2*np.pi;steps=np.linspace(0,T,101)
ring=lambda t:np.c_[np.cos(t),np.sin(t),np.zeros_like(t)]
points=np.array([[0,0,0],[0,0,.5],[0,0,2.]])
r=integrate_kernel(ring,steps,T,points)
z=points[:,2];expected_p=-1/np.sqrt(1+z*z)
expected_a=np.c_[np.zeros_like(z),np.zeros_like(z),-z/(1+z*z)**1.5]
dp=float(max(abs(np.array(r['potential'])-expected_p)))
da=float(np.max(abs(np.array(r['acceleration'])-expected_a)))
assert max(dp,da)<1e-12

def eccentric(t):
    E=t.copy()
    for _ in range(20):E-=(E-.7*np.sin(E)-t)/(1-.7*np.cos(E))
    return np.c_[np.cos(E)-.7,np.sqrt(1-.7**2)*np.sin(E),np.zeros_like(E)]
e=integrate_kernel(eccentric,steps,T,np.array([[0.,0.,0.]]))
error=abs(e['potential'][0]+1)
assert error<1e-10
result=dict(ring_potential_error=dp,ring_acceleration_error=da,kepler_inverse_radius_error=error,
            ring=r,kepler=e,passes=True)
Path(__file__).with_name('analytic-checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n')
print(result)

"""Compare the unchanged finite-stream force kernel in two source coordinates."""
from pathlib import Path
import json
import numpy as np
from scipy.special import roots_legendre
HERE=Path(__file__).resolve().parent
TARGET=np.array([.137,.219]);L=1.;AREA=4.

def polar(target,n):
    # Radial boundary of the square as seen from the interior test point.
    corners=np.array([[-1,-1],[-1,1],[1,-1],[1,1]])-target
    angles=np.sort(np.mod(np.arctan2(corners[:,1],corners[:,0]),2*np.pi))
    bounds=np.r_[angles,angles[0]+2*np.pi]
    z,w=roots_legendre(n);potential=0.;force=np.zeros(2)
    for lo,hi in zip(bounds[:-1],bounds[1:]):
        t=(lo+hi)/2+(hi-lo)/2*z;weight=w*(hi-lo)/2
        co=np.cos(t);si=np.sin(t)
        rx=np.where(co>0,(1-target[0])/co,(-1-target[0])/co)
        ry=np.where(si>0,(1-target[1])/si,(-1-target[1])/si)
        radius=np.minimum(rx,ry)
        potential+=np.sum(weight*(-.5*radius**2*np.arcsinh(L/radius)-.5*L*np.sqrt(radius**2+L**2)+.5*L**2)/L)/AREA
        force+=np.array([sum(weight*co*np.arcsinh(radius/L)),sum(weight*si*np.arcsinh(radius/L))])/AREA
    return dict(potential=float(potential),force=force.tolist(),kernel_nodes=4*n)

def tensor(n):
    z,w=roots_legendre(n)
    x,y=np.meshgrid(z-TARGET[0],z-TARGET[1]);weight=np.outer(w,w)/AREA
    r=np.hypot(x,y)
    p=-np.sum(weight*np.arcsinh(L/r)/L)
    f=np.array([np.sum(weight*x/(r*r*np.sqrt(r*r+L*L))),np.sum(weight*y/(r*r*np.sqrt(r*r+L*L)))])
    return dict(potential=float(p),force=f.tolist(),kernel_nodes=n*n)

reference=polar(TARGET,128);coarse=polar(TARGET,64)
dp=abs(reference['potential']-coarse['potential'])/abs(reference['potential'])
df=float(np.linalg.norm(np.array(reference['force'])-coarse['force'])/np.linalg.norm(reference['force']))
assert max(dp,df)<1e-10
gradient_checks=[]
for h in (1e-4,5e-5):
    fd=[]
    for i in range(2):
        e=np.eye(2)[i]*h
        fd.append(-(polar(TARGET+e,128)['potential']-polar(TARGET-e,128)['potential'])/(2*h))
    error=float(np.linalg.norm(np.array(fd)-reference['force'])/max(np.linalg.norm(reference['force']),.01))
    gradient_checks.append(dict(h=h,force=fd,scaled_error=error))
assert gradient_checks[-1]['scaled_error']<1e-7
rows=[]
for n in (6,8,12,24,48,96):
    r=tensor(n)
    r.update(n=n,potential_relative_error=abs(r['potential']-reference['potential'])/abs(reference['potential']),
             force_relative_error=float(np.linalg.norm(np.array(r['force'])-reference['force'])/np.linalg.norm(reference['force'])))
    rows.append(r)
out=dict(target=TARGET.tolist(),reference=reference,reference_coarse=coarse,
         potential_reference_change=dp,force_reference_change=df,gradient_checks=gradient_checks,tensor=rows,passes=True)
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
print(json.dumps(out,indent=2))

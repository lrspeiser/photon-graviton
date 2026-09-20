"""Energy flux compatible with LR3's exact discrete Hamiltonian allocation."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from local_rotor import curl
from local_transfer import rhs, energy

ROOT=Path(__file__).resolve().parent


def site_energy(state,h,eps,direct):
    a,pi,q,p=state;k=p-direct*a-eps*np.cross(curl(a,h),q)
    value=np.sum(pi*pi+.04*a*a+k*k+q*q,axis=-1)/2
    for j in range(3):
        da=(np.roll(a,-1,j)-a)/h;dq=(np.roll(q,-1,j)-q)/h
        value+=np.sum(da*da+.25*dq*dq,axis=-1)/2
    return value


def face_flux(state,h,eps,direct):
    a,pi,q,p=state;k=p-direct*a-eps*np.cross(curl(a,h),q);s=np.cross(q,k)
    result=[]
    for j in range(3):
        pip=np.roll(pi,-1,j);kp=np.roll(k,-1,j);sp=np.roll(s,-1,j)
        da=(np.roll(a,-1,j)-a)/h;dq=(np.roll(q,-1,j)-q)/h
        result.append(-np.sum(da*pip+.25*dq*kp,axis=-1)+eps/2*(np.cross(pi,sp)+np.cross(pip,s))[...,j])
    return np.stack(result,axis=-1)


def divergence(flux,h):
    return sum((flux[...,j]-np.roll(flux[...,j],1,j))/h for j in range(3))


def outward_power(flux,mask,h):
    weight=mask.astype(float)
    return float(sum(np.sum((weight-np.roll(weight,-1,j))*flux[...,j]) for j in range(3))*h*h)


def main():
    out=ROOT/'discrete-flux-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20261005);rows=[]
    for n,eps,direct in itertools.product((8,12),(0.,2.),(0.,1.)):
        h=8/n;x=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)),indexing='ij'),axis=-1)
        masks=[np.sum(x*x,axis=-1)<1.5**2,np.max(abs(x),axis=-1)<2]
        for fixture in range(3):
            state=.03*rng.normal(size=(4,n,n,n,3));flow=rhs(state,h,eps,direct)
            flux=face_flux(state,h,eps,direct);div=divergence(flux,h)
            scale=max(1.,float(np.max(abs(div))))
            powers=[outward_power(flux,m,h) for m in masks]
            face_errors=[abs(power-float(np.sum(div[m])*h**3))/max(1.,abs(power)) for m,power in zip(masks,powers)]
            residuals=[];balances=[]
            for delta in (1e-6,5e-7):
                hdot=(site_energy(state+delta*flow,h,eps,direct)-site_energy(state-delta*flow,h,eps,direct))/(2*delta)
                residuals.append(float(np.max(abs(hdot+div)))/scale)
                balances.extend(abs(float(np.sum(hdot[m])*h**3)+power)/max(1.,abs(power)) for m,power in zip(masks,powers))
            reference=energy(state,h,eps,direct)
            total_error=abs(float(np.sum(site_energy(state,h,eps,direct))*h**3)-reference)/max(1.,abs(reference))
            periodic=abs(float(np.sum(div)*h**3))
            rows.append(dict(n=n,epsilon=eps,direct=direct,fixture=fixture,local_residuals=residuals,region_residuals=balances,
                             face_errors=face_errors,periodic_residual=periodic,energy_error=total_error,
                             passed=bool(max(residuals+balances)<1e-7 and max(face_errors+[periodic,total_error])<1e-12)))
    filenames=('discrete_flux.py','local_transfer.py','local_rotor.py')
    result=dict(seed=20261005,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in filenames},cases=rows,passed=all(r['passed'] for r in rows))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'],cases=len(rows),maximum_local=max(max(r['local_residuals']) for r in rows),maximum_region=max(max(r['region_residuals']) for r in rows))))
    assert result['passed']


if __name__=='__main__':main()

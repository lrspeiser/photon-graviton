"""Resolved continuum identity, not a finite-difference surface-flux ledger."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent


def derivative(a,axis):
    n=a.shape[0]; shape=[1]*a.ndim;shape[axis]=n
    wave=(2*np.pi*np.fft.fftfreq(n,d=8/n)).reshape(shape)
    return np.fft.ifft(np.fft.fft(a,axis=axis)*(1j*wave),axis=axis).real


def gradients(a):return [derivative(a,j) for j in range(3)]


def curl(a):
    d=gradients(a)
    return np.stack((d[1][...,2]-d[2][...,1],d[2][...,0]-d[0][...,2],d[0][...,1]-d[1][...,0]),axis=-1)


def density(state,eps,direct):
    a,pi,q,p=state;k=p-direct*a-eps*np.cross(curl(a),q)
    value=np.sum(pi*pi+.04*a*a+k*k+q*q,axis=-1)
    for da,dq in zip(gradients(a),gradients(q)):value+=np.sum(da*da+.25*dq*dq,axis=-1)
    return .5*value


def flow_and_flux(state,eps,direct):
    a,pi,q,p=state;b=curl(a);k=p-direct*a-eps*np.cross(b,q)
    da=gradients(a);dq=gradients(q)
    la=sum(derivative(da[j],j) for j in range(3));lq=sum(derivative(dq[j],j) for j in range(3))
    flow=np.stack((pi,la-.04*a+direct*k+eps*curl(np.cross(q,k)),k,.25*lq-q-eps*np.cross(b,k)))
    flux=-np.stack([np.sum(pi*da[j]+.25*k*dq[j],axis=-1) for j in range(3)],axis=-1)+eps*np.cross(pi,np.cross(q,k))
    return flow,flux


def main():
    out=ROOT/'local-energy-flux-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20261004);rows=[]
    modes=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[2,1,0],[1,0,2]])
    for n,eps,direct in itertools.product((24,32),(0.,2.),(0.,1.)):
        x=np.stack(np.meshgrid(*(np.arange(n)*8/n for _ in range(3)),indexing='ij'),axis=-1)
        phase=2*np.pi/8*np.einsum('...j,kj->...k',x,modes)
        for fixture in range(3):
            sine=.003*rng.normal(size=(4,8,3));cosine=.003*rng.normal(size=(4,8,3))
            state=np.einsum('...k,skj->s...j',np.sin(phase),sine)+np.einsum('...k,skj->s...j',np.cos(phase),cosine)
            flow,flux=flow_and_flux(state,eps,direct)
            divergence=sum(derivative(flux[...,j],j) for j in range(3))
            residuals=[];integrals=[]
            for delta in (1e-6,5e-7):
                hdot=(density(state+delta*flow,eps,direct)-density(state-delta*flow,eps,direct))/(2*delta)
                residuals.append(float(np.max(abs(hdot+divergence))/max(1.,np.max(abs(divergence)))))
                integrals.append(float(np.sum(hdot)*(8/n)**3))
            rows.append(dict(n=n,epsilon=eps,direct=direct,fixture=fixture,residuals=residuals,integrated_time_derivatives=integrals,
                             passed=bool(max(residuals)<1e-7 and max(abs(np.array(integrals)))<1e-7)))
    result=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seed=20261004,modes=modes.tolist(),cases=rows,passed=all(r['passed'] for r in rows))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'],cases=len(rows),maximum_residual=max(max(r['residuals']) for r in rows),maximum_integral=max(max(abs(np.array(r['integrated_time_derivatives']))) for r in rows))))
    assert result['passed']


if __name__=='__main__':main()

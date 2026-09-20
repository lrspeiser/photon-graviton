"""SE-H finite-time homogeneous tangent test, not full spatial stability."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
J=np.block([[np.zeros((3,3)),np.eye(3)],[-np.eye(3),np.zeros((3,3))]])


def energy(y,chi,k0):
    phi,x,z=y[:3];p=y[3:];u=.08*phi
    return np.exp(4*u)*np.dot(p,p)/2+.2**2*phi**2/2+.75**2*np.exp(2*chi*u)*(z-k0*np.tanh(u/.001)*x)**2/2


def rhs(y,chi,k0):
    phi,x,z=y[:3];p=y[3:];u=.08*phi;a=np.exp(4*u)
    k=k0*np.tanh(u/.001);ku=k0/.001*(1-np.tanh(u/.001)**2)
    m=.75**2*np.exp(2*chi*u);r=z-k*x
    return np.concatenate((a*p,[-2*.08*a*np.dot(p,p)-.2**2*phi-.08*m*(chi*r*r-ku*x*r),m*k*r,-m*r]))


def jacobian(y,chi,k0):
    columns=[]
    for j in range(6):
        perturbed=y.astype(complex);perturbed[j]+=1e-25j
        columns.append(rhs(perturbed,chi,k0).imag/1e-25)
    return np.column_stack(columns)


def integrate(y,chi,k0,tangent=False):
    def ode(t,z):
        if not tangent:return rhs(z,chi,k0)
        return np.concatenate((rhs(z[:6],chi,k0),(jacobian(z[:6],chi,k0)@z[6:].reshape(6,6)).ravel()))
    initial=np.concatenate((y,np.eye(6).ravel())) if tangent else y
    sol=solve_ivp(ode,(0,20),initial,method='DOP853',rtol=1e-10,atol=1e-12,t_eval=np.linspace(0,20,201))
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def main():
    out=ROOT/'homogeneous-stability-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20260925);controls=[]
    for i in range(48):
        y=rng.normal(0,.01,6);v=rng.normal(size=6);v/=np.linalg.norm(v)
        chi=[0,50,200][i%3];k0=[0,.5][i%2]
        numerical=(energy(y+1e-7*v,chi,k0)-energy(y-1e-7*v,chi,k0))/2e-7
        error=abs(numerical-np.dot(-J@rhs(y,chi,k0),v))
        generator=jacobian(y,chi,k0);sym=float(np.max(abs(generator.T@J+J@generator)))
        controls.append(dict(error=float(error),symplectic_generator_error=sym,passed=bool(error<1e-7 and sym<1e-10)))
    rows=[]
    for chi,k0,channel in itertools.product([0,50,200],[0,.5],['X','Y']):
        y=np.zeros(6);y[1 if channel=='X' else 2]=.01
        sol=integrate(y,chi,k0,True);maps=sol.y[6:].T.reshape(-1,6,6)
        singular=np.linalg.svd(maps,compute_uv=False)
        energies=np.array([energy(z,chi,k0) for z in sol.y[:6].T])
        drift=float(np.max(abs(energies-energies[0])))
        sym=float(np.max(abs(maps[-1].T@J@maps[-1]-J)))
        direction=rng.normal(size=6);direction/=np.linalg.norm(direction);prediction=maps[-1]@direction
        shadows=[]
        for eps in (1e-7,5e-8):
            plus=integrate(y+eps*direction,chi,k0).y[:,-1]
            minus=integrate(y-eps*direction,chi,k0).y[:,-1]
            err=float(np.linalg.norm((plus-minus)/(2*eps)-prediction)/max(1e-12,np.linalg.norm(prediction)))
            shadows.append(dict(epsilon=eps,relative_error=err))
        passed=drift<1e-11+1e-8*abs(energies[0]) and sym<1e-6 and all(s['relative_error']<1e-4 for s in shadows)
        row=dict(chi=chi,k0=k0,channel=channel,initial_energy=float(energies[0]),energy_drift=drift,
                 symplectic_error=sym,maximum_amplification=float(singular[:,0].max()),
                 endpoint_singular_values=singular[-1].tolist(),finite_time_rate=float(np.log(singular[-1,0])/20),
                 shadows=shadows,passed=bool(passed))
        rows.append(row)
        np.savez_compressed(out/f'chi{chi}-mix{k0}-{channel}.npz',time=sol.t,state=sol.y[:6],tangent=maps,energy=energies)
        print(chi,k0,channel,'amplification',row['maximum_amplification'],'passed',row['passed'],flush=True)
    summary=dict(passed=all(r['passed'] for r in rows+controls),count=len(rows),
                 commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                 sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['homogeneous_stability.py','homogeneous-stability-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=summary,controls=controls,rows=rows),indent=2)+'\n')
    print(json.dumps(summary))
    if not summary['passed']:raise RuntimeError('SE-H gate failed; retain first archive')


if __name__=='__main__':main()

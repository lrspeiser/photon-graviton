"""SE-Rot local and grid interaction derivatives; no full evolution claim."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent


def interaction(Q,P,B,coupling):
    K=P-coupling*np.cross(B,Q)
    mass=1+(np.dot(K,K)+np.dot(Q,Q))/2
    return mass,Q+coupling*np.cross(B,K),K,-coupling*np.cross(Q,K)


def curl(A,h):
    def d(a,j):return (np.roll(a,-1,j)-np.roll(a,1,j))/(2*h)
    return np.array([d(A[2],1)-d(A[1],2),d(A[0],2)-d(A[2],0),d(A[1],0)-d(A[0],1)])


def weights(x,q):
    delta=x-q;t=np.maximum(1-np.sum(delta**2,axis=-1)/1.2**2,0)
    raw=t**3;dr=6*t[...,None]**2*delta/1.2**2;norm=raw.sum()
    return raw/norm,dr/norm-raw[...,None]*dr.sum(axis=(0,1,2))/norm**2


def main():
    out=ROOT/'rotor-checks-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20260928);local=[];grid=[]
    for i in range(96):
        Q=rng.normal(0,.03,3);P=rng.normal(0,.03,3);B=rng.normal(0,.1,3);e=(0.,.4,2.,10.)[i%4]
        v=rng.normal(size=9);v/=np.linalg.norm(v);step=1e-6
        M,dQ,dP,dB=interaction(Q,P,B,e)
        numeric=(interaction(Q+step*v[:3],P+step*v[3:6],B+step*v[6:],e)[0]-interaction(Q-step*v[:3],P-step*v[3:6],B-step*v[6:],e)[0])/(2*step)
        err=abs(numeric-np.dot(np.concatenate((dQ,dP,dB)),v))
        torque=float(np.linalg.norm(np.cross(Q,dQ)+np.cross(P,dP)+np.cross(B,dB)))
        local.append(dict(index=i,coupling=e,mass=M,derivative_error=err,torque_identity_error=torque,passed=bool(err<1e-8 and torque<1e-12 and M>=1)))
    n=12;h=8/n;line=np.arange(n)*h-4;x=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
    for i in range(48):
        A=rng.normal(0,.02,(3,n,n,n));q=rng.uniform(-.5,.5,3);Q=rng.normal(0,.03,3);P=rng.normal(0,.03,3);e=(0.,.4,2.,10.)[i%4]
        W,dW=weights(x,q);Bgrid=curl(A,h);B=np.sum(W*Bgrid,axis=(1,2,3))
        M,dQ,dP,dB=interaction(Q,P,B,e)
        gradA=curl(W*dB[:,None,None,None]/h**3,h)
        gradq=np.sum(dW*np.einsum('jxyz,j->xyz',Bgrid,dB)[...,None],axis=(0,1,2))
        v=rng.normal(size=A.size+9);v/=np.linalg.norm(v);vA=v[:A.size].reshape(A.shape);vq=v[A.size:A.size+3];vQ=v[A.size+3:A.size+6];vP=v[-3:]
        def mass(sign):
            step=sign*1e-6;w,_=weights(x,q+step*vq)
            field=curl(A+step*vA,h);b=np.sum(w*field,axis=(1,2,3))
            return interaction(Q+step*vQ,P+step*vP,b,e)[0]
        numeric=(mass(1)-mass(-1))/2e-6
        analytic=np.sum(gradA*vA)*h**3+np.dot(gradq,vq)+np.dot(dQ,vQ)+np.dot(dP,vP)
        err=abs(numeric-analytic)
        C=rng.normal(0,.02,A.shape);adjoint=abs(np.sum(curl(A,h)*C-A*curl(C,h))*h**3)
        grid.append(dict(index=i,coupling=e,derivative_error=err,curl_adjoint_error=adjoint,passed=bool(err<1e-8 and adjoint<1e-12)))
    result=dict(passed=all(r['passed'] for r in local+grid),local_count=len(local),grid_count=len(grid),
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['rotor_checks.py','rotor-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=result,local=local,grid=grid),indent=2)+'\n')
    print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()

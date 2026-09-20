"""SE-DT1: full energy derivatives, frozen propagation, initial source kick."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from density_transfer import DensityTransfer
from joint_transfer import JointTransfer
from check_joint_transfer import state as random_state

ROOT=Path(__file__).resolve().parent


def principal(C,cq,eps,zeta,Q,n,beta):
    a=C*C;b=float(beta@n);T=float(n@Q)*np.eye(3)-np.outer(n,Q);D=np.column_stack((zeta*n,eps*T))
    I=np.eye(4);J=np.eye(3);Z43=np.zeros((4,3));Z34=Z43.T
    M=np.block([[-b*I,-a*D.T,I,b*D.T],[-a*D,-b*J,b*D,cq*cq*J],
                [a*I,Z43,-b*I,Z43],[Z34,a*J,Z34,-b*J]])
    S=np.block([[a*I,Z43,-b*I,Z43],[Z34,a*J,Z34,-b*J],[-b*I,Z43,I,Z43],[Z34,-b*J,Z34,cq*cq*J]])
    return M,S


def main():
    out=ROOT/'density-transfer-v1';out.mkdir(exist_ok=False);rng=np.random.default_rng(20261008);derivatives=[];waves=[];kicks=[]
    for n,zeta in itertools.product((8,12),(0.,.5,2.)):
        model=DensityTransfer(n=n,zeta=zeta);base=JointTransfer(n=n)
        for fixture in range(8):
            y=random_state(rng,model);v=rng.normal(size=y.shape);v/=np.linalg.norm(v)
            fd,pd,xd,ppd=model.unpack(model.evaluate(y));df,dp,dx,dpp=model.unpack(v)
            expected=float(model.dv*np.sum(fd*dp-pd*df)+np.sum(xd*dpp-ppd*dx))
            numerical=[(model.evaluate(y+d*v,False)['energy']-model.evaluate(y-d*v,False)['energy'])/(2*d) for d in (1e-6,5e-7)]
            errors=[abs(value-expected)/max(1.,abs(expected)) for value in numerical];difference=abs(numerical[0]-numerical[1])/max(1.,abs(expected))
            reduction=None
            if zeta==0:reduction=float(max(np.max(abs(model.evaluate(y)-base.evaluate(y))),abs(model.evaluate(y,False)['energy']-base.evaluate(y,False)['energy'])))
            derivatives.append(dict(n=n,zeta=zeta,fixture=fixture,analytic=expected,numerical=numerical,errors=errors,step_difference=difference,reduction_error=reduction,
                                    passed=bool(max(errors+[difference])<2e-7 and (reduction is None or reduction<1e-12))))
    for C,eps,zeta,fraction in itertools.product((.5,2.),(0.,2.),(0.,.5,2.),(0.,.9)):
        for fixture in range(8):
            n,Q,beta=rng.normal(size=(3,3));n/=np.linalg.norm(n);Q*=.03/np.linalg.norm(Q);beta*=fraction*.5*C/np.linalg.norm(beta)
            M,S=principal(C,.5,eps,zeta,Q,n,beta);spectrum=np.linalg.eigvals(M);speeds=-spectrum.real;b=float(beta@n)
            minimum=float(np.linalg.eigvalsh(S)[0]);symmetry=float(np.max(abs(S@M-M.T@S))/max(1.,np.max(abs(S@M))));imaginary=float(np.max(abs(spectrum.imag)))
            longitudinal=[];direction=np.concatenate(([0.],n))
            for sign in (-1.,1.):
                vector=np.concatenate((direction,np.zeros(3),sign*C*direction,np.zeros(3)))
                longitudinal.append(float(np.max(abs(M@vector-(-b+sign*C)*vector))))
            coverage=float(max(0.,min(speeds)-(b-C),(b+C)-max(speeds)))
            waves.append(dict(C=C,epsilon=eps,zeta=zeta,beta_fraction=fraction,fixture=fixture,minimum_S=minimum,symmetry=symmetry,imaginary=imaginary,longitudinal=longitudinal,coverage=coverage,
                              passed=minimum>0 and symmetry<1e-12 and imaginary<1e-10 and max(longitudinal)<1e-10 and coverage<1e-10))
    for n,zeta,speed in itertools.product((8,12),(0.,.5,2.),(0.,.0007,.2)):
        model=DensityTransfer(n=n,zeta=zeta,masses=[1.]*6);h=model.h;theta=np.arange(6)*np.pi/3
        q=.7*np.column_stack((np.cos(theta),np.sin(theta),np.zeros(6)));p=speed*np.column_stack((-np.sin(theta),np.cos(theta),np.zeros(6)))
        initial=np.concatenate((np.zeros(2*model.nf),q.ravel(),p.ravel()))
        xyz=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)),indexing='ij'),axis=-1);rho=np.zeros((n,n,n));current=np.zeros((n,n,n,3))
        for position,momentum in zip(q,p):
            d=abs((xyz-position)/h);basis=np.where(d<1,(4-6*d*d+3*d**3)/6,np.where(d<2,(2-d)**3/6,0.));w=np.prod(basis,axis=-1)/model.dv
            pp=momentum@momentum;E=np.sqrt(1+pp);rho+=w*(E+pp/E);current+=w[...,None]*momentum
        gradient=np.stack([(np.take(rho,(np.arange(n)+1)%n,axis=j)-np.take(rho,(np.arange(n)-1)%n,axis=j))/(2*h) for j in range(3)],axis=-1)
        predicted=model.direct*model.kappa*model.cq*model.eta*current+zeta*model.g*gradient
        errors=[];observed=[]
        for dt in (1e-4,5e-5):
            k1=model.evaluate(initial);k2=model.evaluate(initial+.5*dt*k1);k3=model.evaluate(initial+.5*dt*k2);k4=model.evaluate(initial+dt*k3)
            final=initial+dt*(k1+2*k2+2*k3+k4)/6;f,_,_,_=model.unpack(final);estimate=6*np.moveaxis(f[4:7],0,-1)/dt**3
            errors.append(float(np.max(abs(estimate-predicted))/max(1.,np.max(abs(predicted)))));observed.append(float(np.max(abs(estimate))))
        kicks.append(dict(n=n,zeta=zeta,speed=speed,predicted_max=float(np.max(abs(predicted))),observed_max=observed,errors=errors,passed=max(errors)<1e-4))
    files=('density_transfer.py','check_density_transfer.py','joint_transfer.py','point_probes.py','local_rotor.py','check_joint_transfer.py')
    result=dict(seed=20261008,commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
                derivatives=derivatives,principal=waves,kicks=kicks,passed=all(r['passed'] for r in derivatives+waves+kicks))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'],derivatives=len(derivatives),principal=len(waves),kicks=len(kicks),max_gradient=max(max(r['errors']) for r in derivatives),max_kick=max(max(r['errors']) for r in kicks),rest_kicks=[r for r in kicks if r['speed']==0.])))
    assert result['passed']


if __name__=='__main__':main()

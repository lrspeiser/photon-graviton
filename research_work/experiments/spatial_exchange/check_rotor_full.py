"""SE-Rot2 full Hamiltonian derivatives and zero-coupling scalar reduction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from rotor_model import ExchangeEvolution as Rotor
from model import ExchangeEvolution as Scalar

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'rotor-full-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20260929);rows=[]
    for i in range(48):
        cfg=dict(n=12,length=10.,radius=1.2,chi=(0.,200.)[(i//4)%2],rotor_coupling=(0.,.4,2.,10.)[i%4],spin=(1,-1)[(i//8)%2])
        model=Rotor(**cfg);y=model.initial();f,pi,q,p,Q,P=model.unpack(y)
        env=np.exp(-np.sum(model.x**2,axis=-1)/2)
        f[:]=rng.normal(0,.03,f.shape)*env;pi[:]=rng.normal(0,.02,pi.shape)*env
        q[:]+=rng.normal(0,.03,q.shape);Q[:]+=rng.normal(0,.03,Q.shape);P[:]+=rng.normal(0,.03,P.shape)
        fd,pd,qd,ppd,Qd,Pd=model.unpack(model.rhs(y))
        gradient=np.concatenate(((-model.dv*(pd+model.gamma*fd)).ravel(),(model.dv*fd).ravel(),(-ppd).ravel(),qd.ravel(),(-Pd).ravel(),Qd.ravel(),[0.]))
        v=rng.normal(size=y.size);v[-1]=0;v/=np.linalg.norm(v);step=2e-5
        numerical=(model.energy(y+step*v)-model.energy(y-step*v))/(2*step);exact=np.dot(gradient,v)
        error=abs(numerical-exact)/max(1.,abs(exact));metrics=model.metrics(y)
        reduction=None
        if cfg['rotor_coupling']==0:
            Q[:,1:]=0;P[:,1:]=0
            old=Scalar(n=12,length=10.,radius=1.2,chi=cfg['chi'],emission=0.)
            oldy=np.concatenate((y[:2*model.nf+36],Q[:,0],P[:,0],[y[-1]]))
            rd=model.rhs(y);od=old.rhs(oldy);rf,rpi,rq,rp,rQ,rP=model.unpack(rd);of,opi,oq,op,oQ,oP=old.unpack(od)
            reduction=max(abs(model.energy(y)-old.energy(oldy)),float(np.max(abs(rd[:2*model.nf+36]-od[:2*model.nf+36]))),float(np.max(abs(rQ[:,0]-oQ))),float(np.max(abs(rP[:,0]-oP))),float(np.max(abs(rQ[:,1:]))),float(np.max(abs(rP[:,1:]))))
        rows.append(dict(index=i,config=cfg,gradient_error=error,minimum_mass=metrics['minimum_mass'],cone_error=metrics['cone_error'],reduction_error=reduction,
                         passed=bool(error<2e-6 and metrics['minimum_mass']>0 and metrics['cone_error']<1e-10 and (reduction is None or reduction<1e-12))))
    result=dict(passed=all(r['passed'] for r in rows),count=len(rows),
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['rotor_model.py','check_rotor_full.py','rotor-full-protocol.md','model.py']})
    (out/'results.json').write_text(json.dumps(dict(summary=result,rows=rows),indent=2)+'\n')
    print(json.dumps(result));print('max gradient',max(r['gradient_error'] for r in rows),'max reduction',max(r['reduction_error'] or 0 for r in rows))
    assert result['passed']


if __name__=='__main__':main()

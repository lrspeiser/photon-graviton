"""SE-JT1 full canonical pairing and zero-feedback reductions."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from joint_transfer import JointTransfer
from local_transfer import rhs as local_rhs
from local_rotor import laplacian

ROOT=Path(__file__).resolve().parent


def state(rng,model):
    fields=.03*rng.normal(size=2*model.nf)
    positions=rng.uniform(-.6,.6,size=(model.count,3));momenta=rng.normal(size=(model.count,3))*np.array([.2,.2,.01])[:,None]
    return np.concatenate((fields,positions.ravel(),momenta.ravel()))


def main():
    out=ROOT/'joint-transfer-v1';out.mkdir(exist_ok=False);rng=np.random.default_rng(20261007);rows=[];reductions=[]
    for n,eps,direct in itertools.product((8,12),(0.,2.),(0.,1.)):
        model=JointTransfer(n=n,epsilon=eps,direct=direct)
        for fixture in range(6):
            y=state(rng,model);v=rng.normal(size=y.shape);v/=np.linalg.norm(v)
            fd,pd,xd,ppd=model.unpack(model.evaluate(y));df,dp,dx,dpp=model.unpack(v)
            exact=float(model.dv*np.sum(-pd*df+fd*dp)+np.sum(-ppd*dx+xd*dpp))
            numerics=[]
            for delta in (1e-6,5e-7):numerics.append((model.evaluate(y+delta*v,False)['energy']-model.evaluate(y-delta*v,False)['energy'])/(2*delta))
            errors=[abs(value-exact)/max(1.,abs(exact)) for value in numerics];difference=abs(numerics[0]-numerics[1])/max(1.,abs(exact));d=model.evaluate(y,False)
            rows.append(dict(n=n,epsilon=eps,direct=direct,fixture=fixture,analytic=exact,numerical=numerics,errors=errors,step_difference=difference,
                             energy=d['energy'],minimum_density=d['minimum_density'],cone_error=d['cone_error'],
                             passed=bool(max(errors+[difference])<2e-7 and d['energy']>0 and d['cone_error']<1e-12)))
    for eps,direct in itertools.product((0.,2.),(0.,1.)):
        model=JointTransfer(n=8,epsilon=eps,direct=direct,g=0.,eta=0.)
        for fixture in range(3):
            y=state(rng,model);f,m,q,p=model.unpack(y);fd,pd,xd,ppd=model.unpack(model.evaluate(y))
            local=np.stack([np.moveaxis(s,0,-1) for s in (f[1:4],m[1:4],f[4:7],m[4:7])]);expected=local_rhs(local,model.h,eps,direct)
            actual=np.stack([np.moveaxis(s,0,-1) for s in (fd[1:4],pd[1:4],fd[4:7],pd[4:7])])
            error=float(max(np.max(abs(actual-expected)),np.max(abs(fd[0]-m[0])),np.max(abs(pd[0]-laplacian(f[0],model.h)+model.omega**2*f[0]))))
            speed=p/np.sqrt(model.masses**2+np.sum(p*p,axis=1))[:,None]
            particle_error=float(max(np.max(abs(xd-speed)),np.max(abs(ppd))))
            reductions.append(dict(epsilon=eps,direct=direct,fixture=fixture,field_error=error,particle_error=particle_error,passed=error<1e-12 and particle_error<1e-12))
    result=dict(seed=20261007,commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                hashes={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('joint_transfer.py','check_joint_transfer.py','point_probes.py','local_rotor.py','local_transfer.py')},
                derivatives=rows,reductions=reductions,passed=all(r['passed'] for r in rows+reductions))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'],derivatives=len(rows),reductions=len(reductions),maximum_derivative_error=max(max(r['errors']) for r in rows),maximum_reduction_error=max(max(r['field_error'],r['particle_error']) for r in reductions))))
    assert result['passed']


if __name__=='__main__':main()

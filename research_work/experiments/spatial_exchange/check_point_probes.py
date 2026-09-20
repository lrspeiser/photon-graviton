"""SE-PI point-interpolation derivative and continuum-consistency controls."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from point_probes import evaluate,kernel

ROOT=Path(__file__).resolve().parent


def analytic(q):
    env=np.exp(-np.sum(q*q,axis=0)/8);u=-.0032*env
    A=.1*env*np.array([-q[1],q[0],0*q[0]])
    return np.exp(u),np.exp(2*u),.04*np.exp(2*u)*A/np.sqrt(1+.08**2*np.sum(A*A,axis=0))


def grid(n):
    line=np.arange(n)*16/n-8
    return analytic(np.array(np.meshgrid(line,line,line,indexing='ij')))


def exact(q,p,m):
    al,z,beta=analytic(q);E=np.sqrt(m*m+z*np.dot(p,p));velocity=al*z*p/E+beta
    force=[]
    for j in range(3):
        pert=q.astype(complex);pert[j]+=1e-25j;a,b,c=analytic(pert)
        force.append(-(a*np.sqrt(m*m+b*np.dot(p,p))+np.dot(c,p)).imag/1e-25)
    return np.concatenate((velocity,force))


def main():
    out=ROOT/'point-probe-controls-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20261001);rows=[];grids={n:grid(n) for n in (24,32)}
    for i in range(72):
        n=(24,32)[(i//3)%2];m=(0.,.1,1.)[i%3];q=rng.uniform(-2,2,3);p=rng.normal(size=3)
        arrays=grids[n];value=evaluate(q,p,m,16,*arrays);v=rng.normal(size=6);v/=np.linalg.norm(v);eps=1e-5
        numerical=(evaluate(q+eps*v[:3],p+eps*v[3:],m,16,*arrays)['energy']-evaluate(q-eps*v[:3],p-eps*v[3:],m,16,*arrays)['energy'])/(2*eps)
        derivative_error=abs(numerical-np.dot(np.concatenate((-value['force'],value['velocity'])),v))
        speed=np.linalg.norm(value['velocity']-value['bbar'])/value['cbar'];cone=max(0.,speed-1);light=abs(speed-1) if m==0 else 0.
        _,w,dw=kernel(q,n,16);partition=abs(w.sum()-1);dsum=float(np.max(abs(dw.sum(axis=(1,2,3)))))
        uniform=evaluate(q,p,m,16,np.ones((n,n,n)),np.ones((n,n,n)),np.zeros((3,n,n,n)))
        flat_force=float(np.max(abs(uniform['force'])))
        passed=derivative_error<1e-7 and cone<1e-12 and light<1e-12 and partition<1e-12 and dsum<1e-12 and flat_force<1e-12
        rows.append(dict(index=i,n=n,mass=m,derivative_error=derivative_error,cone_excess=cone,light_error=light,partition_error=partition,derivative_sum_error=dsum,uniform_force=flat_force,passed=bool(passed)))
    tests=[];points=[np.array(q) for q in [(.23,.37,.19),(1.13,-.67,.41),(-1.71,.83,-.29)]];p=np.array([.7,.2,-.1])
    for q in points:
        for m in (0.,.1,1.):tests.append(dict(position=q.tolist(),mass=m,errors=[]))
    for n in (24,48,96):
        arrays=grid(n)
        for test in tests:
            q=np.array(test['position']);m=test['mass'];value=evaluate(q,p,m,16,*arrays)
            error=float(np.linalg.norm(np.concatenate((value['velocity'],value['force']))-exact(q,p,m)))
            test['errors'].append(dict(n=n,error=error))
    for test in tests:test['passed']=test['errors'][-1]['error']<test['errors'][0]['error']/8
    result=dict(passed=all(r['passed'] for r in rows+tests),controls=len(rows),continuum_fixtures=len(tests),
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['point_probes.py','check_point_probes.py','point-probe-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=result,rows=rows,continuum=tests),indent=2)+'\n')
    print(json.dumps(result));print('minimum coarse/fine error ratio',min(t['errors'][0]['error']/t['errors'][-1]['error'] for t in tests));assert result['passed']


if __name__=='__main__':main()

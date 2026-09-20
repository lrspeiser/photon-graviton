"""SE-O seeded observer convention controls; preserve first execution."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from probes import coefficients,evaluate
from observers import observer,photon_energy

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'observer-controls-v2';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20260927);rows=[];uniform=[]
    for i in range(48):
        n=(24,32)[i%2];radius=(.9,1.2)[(i//2)%2];v=np.array((0.,0.,0.) if (i//4)%2==0 else (.1,-.03,.02))
        line=np.arange(n)*16/n-8;xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=0)
        env=np.exp(-np.sum(xyz*xyz,axis=0)/8);f=np.zeros((6,n,n,n));f[0]=-.04*env;f[1]=-.1*xyz[1]*env;f[2]=.1*xyz[0]*env
        arrays=coefficients(f);q=rng.uniform(-2,2,3);k=rng.normal(size=3);obs=observer(q,v,radius,16,arrays);p=obs['momentum']
        eps=1e-5
        derivative=(evaluate(q,p,1+eps,radius,16,*arrays)['energy']-evaluate(q,p,1-eps,radius,16,*arrays)['energy'])/(2*eps)
        mass_error=abs(derivative-obs['rate']);euler=abs(obs['energy']-np.dot(p,obs['velocity'])-obs['rate'])
        measured=photon_energy(q,k,obs,radius,16,arrays);scaling=abs(photon_energy(q,2*k,obs,radius,16,arrays)-2*measured)
        rows.append(dict(index=i,n=n,radius=radius,clock_rate=obs['rate'],velocity_residual=obs['residual'],solver_success=obs['solver_success'],mass_derivative_error=mass_error,euler_error=euler,photon_energy=measured,scaling_error=scaling,
                         passed=bool(obs['residual']<1e-10 and obs['rate']>0 and mass_error<1e-8 and euler<1e-12 and scaling<1e-12)))
    for shifted in (False,True):
        f=np.zeros((6,24,24,24))
        if shifted:f[0]=-.02;f[1]=.1;f[2]=-.05
        arrays=coefficients(f);obs=observer(np.zeros(3),np.zeros(3),.9,16,arrays)
        al,z,beta=arrays;expected=float(al[0,0,0]*np.sqrt(1-np.sum(beta[:,0,0,0]**2)/z[0,0,0]**2))
        rate_error=abs(obs['rate']-expected);energy_error=0.
        if not shifted:energy_error=abs(photon_energy(np.zeros(3),np.array([1.,2.,3.]),obs,.9,16,arrays)-np.sqrt(14))
        uniform.append(dict(shifted=shifted,rate_error=rate_error,energy_error=energy_error,passed=bool(rate_error<1e-10 and energy_error<1e-10)))
    summary=dict(passed=all(r['passed'] for r in rows+uniform),fixtures=len(rows),uniform=len(uniform),
                 commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                 sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['observers.py','check_observers.py','probes.py','observer-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=summary,rows=rows,uniform=uniform),indent=2)+'\n')
    print(json.dumps(summary));assert summary['passed']


if __name__=='__main__':main()

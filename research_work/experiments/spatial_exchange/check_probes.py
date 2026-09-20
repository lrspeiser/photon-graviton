"""Preregistered probe derivative, cone, scaling and uniform-flight controls."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from probes import coefficients,evaluate,derivative

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'probe-controls-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20260926);rows=[];flights=[]
    for i in range(72):
        n=(24,32)[(i//6)%2];radius=(.9,1.2)[(i//3)%2];mass=(0.,.1,1.)[i%3]
        line=np.arange(n)*16/n-8;xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=0)
        env=np.exp(-np.sum(xyz*xyz,axis=0)/8);fields=np.zeros((6,n,n,n))
        fields[0]=-.04*env;fields[1]=-.1*xyz[1]*env;fields[2]=.1*xyz[0]*env
        coeff=coefficients(fields);q=rng.uniform(-2,2,3);p=rng.normal(size=3)
        value=evaluate(q,p,mass,radius,16,*coeff);direction=rng.normal(size=6);direction/=np.linalg.norm(direction)
        eps=1e-5
        plus=evaluate(q+eps*direction[:3],p+eps*direction[3:],mass,radius,16,*coeff)['energy']
        minus=evaluate(q-eps*direction[:3],p-eps*direction[3:],mass,radius,16,*coeff)['energy']
        exact=np.dot(-value['force'],direction[:3])+np.dot(value['velocity'],direction[3:])
        error=abs((plus-minus)/(2*eps)-exact)
        relative_speed=np.linalg.norm(value['velocity']-value['bbar'])/value['cbar']
        cone=max(0.,relative_speed-1);scaling=0.;equality=0.
        if mass==0:
            equality=abs(relative_speed-1)
            for scale in (.5,2):
                other=evaluate(q,scale*p,mass,radius,16,*coeff)
                scaling=max(scaling,float(np.max(abs(other['velocity']-value['velocity']))),float(np.max(abs(other['force']-scale*value['force']))),abs(other['energy']-scale*value['energy']))
        rows.append(dict(index=i,mass=mass,n=n,radius=radius,derivative_error=error,cone_excess=cone,
                         light_speed_error=equality,scaling_error=scaling,
                         passed=bool(error<1e-7 and cone<1e-12 and equality<1e-12 and scaling<1e-12)))
    for mass in (0.,1.):
        for shifted in (False,True):
            fields=np.zeros((6,24,24,24));fields[0]=-.02
            if shifted:fields[1]=.1;fields[2]=-.05
            coeff=coefficients(fields);state=np.array([[-2.,0.,0.,.4,.1,0.]])
            initial=state.copy();value=evaluate(state[0,:3],state[0,3:],mass,.9,16,*coeff)
            for step in range(200):
                k1=derivative(state,[mass],.9,16,coeff);k2=derivative(state+.01*k1,[mass],.9,16,coeff)
                k3=derivative(state+.01*k2,[mass],.9,16,coeff);k4=derivative(state+.02*k3,[mass],.9,16,coeff)
                state+=.02*(k1+2*k2+2*k3+k4)/6
            position=float(np.linalg.norm(state[0,:3]-initial[0,:3]-4*value['velocity']))
            momentum=float(np.linalg.norm(state[0,3:]-initial[0,3:]))
            flights.append(dict(mass=mass,shifted=shifted,position_error=position,momentum_error=momentum,passed=position<1e-10 and momentum<1e-12))
    summary=dict(passed=all(r['passed'] for r in rows+flights),fixtures=len(rows),flights=len(flights),
                 commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                 sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['probes.py','check_probes.py','probe-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=summary,rows=rows,flights=flights),indent=2)+'\n')
    print(json.dumps(summary));assert summary['passed']


if __name__=='__main__':main()

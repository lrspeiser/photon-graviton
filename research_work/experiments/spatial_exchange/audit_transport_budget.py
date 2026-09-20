"""Reconstruct region energies and integrate archived stage powers separately."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def independent_site(state,h):
    a,pi,q,p=state;n=a.shape[0];i=np.arange(n)
    forward=lambda f,j:(np.take(f,(i+1)%n,axis=j)-f)/h
    center=lambda f,j:(np.take(f,(i+1)%n,axis=j)-np.take(f,(i-1)%n,axis=j))/(2*h)
    d=[center(a,j) for j in range(3)]
    b=np.stack((d[1][...,2]-d[2][...,1],d[2][...,0]-d[0][...,2],d[0][...,1]-d[1][...,0]),axis=-1)
    k=p-a-2*np.cross(b,q)
    values=np.sum(pi*pi+.04*a*a+k*k+q*q,axis=-1)/2
    for j in range(3):values+=np.sum(forward(a,j)**2+.25*forward(q,j)**2,axis=-1)/2
    return values


def main():
    out=ROOT/'transport-budget-v1';summary=json.loads((out/'summary.json').read_text());manifest=json.loads((out/'manifest.json').read_text())
    checks=[hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==v for f,v in manifest['hashes'].items()]
    maximum=0.
    for row in summary['runs']:
        archive=np.load(out/(row['name']+'.npz'));trace=archive['trace'];stages=archive['stage_powers'];n=row['n'];h=8/n
        x=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)),indexing='ij'),axis=-1);masks=[np.sum(x*x,axis=-1)<r*r for r in (1.2,1.5,2.)]
        checks.append(row==json.loads((out/(row['name']+'.json')).read_text()))
        for key,t in (('initial',0),('final',-1)):
            density=independent_site(archive[key],h)
            energies=np.array([np.sum(density[m])*h**3 for m in masks]);error=float(np.max(abs(energies-trace[t,1:4])));maximum=max(maximum,error);checks.append(error<1e-12)
            if key=='initial':checks.append(abs(np.sum(density)*h**3-row['initial_total_energy'])<1e-12)
        increments=np.einsum('skj,k->sj',stages,np.array([1.,2.,2.,1.]))*row['dt']/6
        cumulative=np.vstack((np.zeros(3),np.cumsum(increments,axis=0)))
        checks.extend((np.isfinite(stages).all(),np.max(abs(cumulative-trace[:,4:7]))<1e-12,
                       np.max(abs(cumulative[-1]-row['net_outward_energy']))<1e-12,
                       np.max(abs(cumulative[-1]/row['initial_total_energy']-row['net_outward_fractions']))<1e-12))
        residual=float(np.max(abs(trace[:,1:4]-trace[0,1:4]+cumulative))/row['initial_total_energy'])
        replay=float(np.max(abs(archive['final']-np.load(ROOT/'local-transfer-v1'/(row['name']+'.npz'))['final'])))
        checks.extend((abs(residual-row['balance_error'])<1e-12,replay==row['replay_error'],row['passed']==(residual<1e-6 and replay<1e-12)))
    coarse,fine=summary['runs'][:2];error=abs(coarse['net_outward_energy'][1]-fine['net_outward_energy'][1])/max(abs(fine['net_outward_energy'][1]),1e-12)
    checks.extend((abs(error-summary['time_relative_difference'])<1e-12,summary['passed']==(all(r['passed'] for r in summary['runs']) and error<.001)))
    result=dict(passed=bool(all(checks)),checks=len(checks),maximum_endpoint_error=maximum,campaign_passed=summary['passed'],
                scope='Endpoint regional energies independently reconstructed; integrated transport recomputed from archived powers, without independent full-stage fields.')
    (ROOT/'transport-budget-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()

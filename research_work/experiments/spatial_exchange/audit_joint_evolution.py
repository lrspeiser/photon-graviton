"""Independent joint endpoint Hamiltonian and full-grid particle reconstruction."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
NAMES=('free','baseline','combined','slow-baseline','slow-combined','time','space','rotation')


def reconstruct(y,cfg):
    n=cfg['n'];h=10/n;dv=h**3;nf=7*n**3
    f=y[:nf].reshape(7,n,n,n);m=y[nf:2*nf].reshape(f.shape);q,p=y[2*nf:].reshape(2,7,3)
    index=np.arange(n)
    grad=[(np.take(f,(index+1)%n,axis=j+1)-np.take(f,(index-1)%n,axis=j+1))/(2*h) for j in range(3)]
    A,Pi,Q,P=[np.moveaxis(s,0,-1) for s in (f[1:4],m[1:4],f[4:7],m[4:7])]
    B=np.stack((grad[1][3]-grad[2][2],grad[2][1]-grad[0][3],grad[0][2]-grad[1][1]),axis=-1)
    K=P-cfg['direct']*A-cfg['epsilon']*np.cross(B,Q)
    g=0. if cfg['free'] else .08;eta=0. if cfg['free'] else .08;C=np.exp(2*g*f[0]);a=C*C
    beta=.25*eta*C[...,None]*A/np.sqrt(1+eta**2*np.sum(A*A,axis=-1))[...,None]
    effective=m.copy();effective[4:7]=np.moveaxis(K,-1,0)
    density=.5*a*np.sum(effective*effective,axis=0)+.02*np.sum(f[:4]**2,axis=0)+.5*np.sum(f[4:7]**2,axis=0)
    eq=.5*a*np.sum(K*K,axis=-1)+.5*np.sum(Q*Q,axis=-1)
    momentum=np.zeros((n,n,n,3))
    for j in range(3):
        plus=(np.take(f,(index+1)%n,axis=j+1)-f)/h
        density+=.5*np.sum(plus[:4]**2,axis=0)+.125*np.sum(plus[4:7]**2,axis=0)-beta[...,j]*np.sum(effective*grad[j],axis=0)
        eq+=.125*np.sum(plus[4:7]**2,axis=0)-beta[...,j]*np.sum(effective[4:7]*grad[j][4:7],axis=0)
        momentum[...,j]=-np.sum(m*grad[j],axis=0)
    xyz=np.stack(np.meshgrid(*(index*h-5 for _ in range(3)),indexing='ij'),axis=-1)
    angular=np.sum(np.cross(xyz,momentum)+np.cross(A,Pi)+np.cross(Q,P),axis=(0,1,2))*dv+np.sum(np.cross(q,p),axis=0)
    momtotal=np.sum(momentum,axis=(0,1,2))*dv+np.sum(p,axis=0)
    matter=0.;velocity=None;cone=0.
    for i,mass in enumerate([1.]*6+[0.]):
        distance=abs((xyz-q[i])/h);basis=np.where(distance<1,(4-6*distance**2+3*distance**3)/6,np.where(distance<2,(2-distance)**3/6,0.))
        w=np.prod(basis,axis=-1);E=np.sqrt(mass*mass+C*(p[i]@p[i]))
        matter+=float(np.sum(w*(np.sqrt(C)*E+np.sum(beta*p[i],axis=-1))))
        bbar=np.sum(w[...,None]*beta,axis=(0,1,2));v=p[i]*np.sum(w*np.sqrt(C)*C/E)+bbar;cbar=np.sum(w*C)
        cone=max(cone,abs(np.linalg.norm(v-bbar)-cbar) if mass==0 else max(0.,np.linalg.norm(v-bbar)-cbar))
        if i==6:velocity=v
    field=float(np.sum(density)*dv)
    return dict(energy=field+matter,field=field,matter=matter,KQ_energy=float(np.sum(eq)*dv),angular=angular,momentum=momtotal,photon_velocity=velocity,cone_error=float(cone))


def main():
    out=ROOT/'joint-evolution-v1';manifest=json.loads((out/'manifest.json').read_text())
    checks=[hashlib.sha256((ROOT/f).read_bytes()).hexdigest()==v for f,v in manifest['hashes'].items()];rows=[];maximum=0.
    for name in NAMES:
        if not (out/(name+'.json')).exists():continue
        row=json.loads((out/(name+'.json')).read_text());rows.append(row);raw=np.load(out/(name+'.npz'));trace=row['trace']
        checks.append(len(trace)==round(2/row['dt'])+1 and trace[-1]['time']==2.)
        for key,t in (('initial',0),('final',-1)):
            d=reconstruct(raw[key],row)
            error=max(float(np.max(abs(np.array(d[k])-trace[t][k]))) for k in d);maximum=max(maximum,error);checks.append(error<1e-10)
            n=row['n'];nf=7*n**3;checks.append(np.max(abs(raw[key][2*nf:].reshape(2,7,3)-raw['particles'][t]))<1e-12)
        checks.append(np.max(abs(raw['initial'][:14*row['n']**3]))==0.)
        if name in ('baseline','slow-baseline'):
            n=row['n'];nf=7*n**3;fields=raw['final'][:2*nf].reshape(2,7,n,n,n);checks.append(np.max(abs(fields[:,4:7]))<1e-12)
        e=max(abs(t['energy']-trace[0]['energy']) for t in trace)/trace[0]['energy'];j0=np.array(trace[0]['angular'])
        j=max(np.linalg.norm(np.array(t['angular'])-j0) for t in trace)/max(np.linalg.norm(j0),1e-8)
        expected=dict(finite=bool(np.isfinite(raw['final']).all()),energy=e<1e-5,angular=bool(j<.01),cone=max(t['cone_error'] for t in trace)<1e-10,edge=max(t['edge_amplitude'] for t in trace)<1e-5)
        if row['free']:
            q0,p0=raw['particles'][0];q,p=raw['particles'][-1];v0=p0/np.sqrt(np.array([1.]*6+[0.])**2+np.sum(p0*p0,axis=1))[:,None]
            err=float(max(np.max(abs(q-q0-2*v0)),np.max(abs(p-p0))));checks.append(abs(err-row['free_error'])<1e-12)
            expected['free']=np.max(abs(raw['final'][:14*row['n']**3]))<1e-12 and err<1e-10
        checks.extend((expected==row['checks'],all(expected.values())==row['passed'],abs(e-row['energy_drift'])<1e-12,abs(j-row['angular_drift'])<1e-12))
        rot=np.array(row['rotation']);v=np.array(trace[-1]['photon_velocity']);v0=np.array(trace[0]['photon_velocity'])
        bend=float(np.arctan2(v@rot[:,1],v@rot[:,0])-np.arctan2(v0@rot[:,1],v0@rot[:,0]));checks.append(abs(bend-row['bend'])<1e-12)
    complete=(out/'summary.json').exists();campaign=None
    if complete:
        summary=json.loads((out/'summary.json').read_text());checks.append(summary['runs']==[{k:v for k,v in r.items() if k!='trace'} for r in rows]);by={r['name']:r for r in rows};base=by['combined'];comparisons=[]
        for name,limit in (('time',.001),('space',.05)):
            error=abs(base['bend']-by[name]['bend'])/max(abs(by[name]['bend']),1e-8);comparisons.append(dict(name=name,relative_error=error,passed=error<limit))
        v=np.array(by['rotation']['trace'][-1]['photon_velocity'])@np.array(by['rotation']['rotation']);u=np.array(base['trace'][-1]['photon_velocity'])
        error=float(np.arctan2(np.linalg.norm(np.cross(v,u)),v@u)/max(abs(base['bend']),1e-8));comparisons.append(dict(name='rotation',relative_error=error,passed=error<.02))
        checks.append(comparisons==summary['comparisons']);campaign=summary['passed'];checks.append(campaign==all(r['passed'] for r in rows+comparisons))
    result=dict(passed=bool(all(checks)),checks=len(checks),completed_runs=len(rows),declared_runs=8,complete=complete,campaign_passed=campaign,maximum_endpoint_error=maximum)
    (ROOT/'joint-evolution-audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()

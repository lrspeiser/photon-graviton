"""CC-2 independent endpoint energy reconstruction and evidence audit."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import hashlib
import subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=HERE/'evolution-v1'


def independent_energy(state,spec):
    n=spec.get('n',32);length=spec.get('length',12);dx=length/n;size=4*n**3
    count=6 if spec['probe']=='none' else 7
    f=state[:size].reshape((4,)+(n,)*3);pi=state[size:2*size].reshape(f.shape)
    q=state[2*size:2*size+3*count].reshape(-1,3);p=state[2*size+3*count:-1].reshape(-1,3)
    alpha=np.exp(.08*f[0]);eta=spec['eta'];a=f[1:]
    beta=alpha*.5*eta*a/np.sqrt(1+eta*eta*np.sum(a*a,axis=0))
    field=np.sum(alpha*np.sum(pi*pi,axis=0)/2+.2**2*np.sum(f*f,axis=0)/2)
    for axis in range(3):
        forward=(np.roll(f,-1,axis=axis+1)-f)/dx
        backward=(f-np.roll(f,1,axis=axis+1))/dx
        field+=np.sum(alpha*np.sum(forward**2+backward**2,axis=0)/4)
        field-=np.sum(beta[axis]*np.sum(pi*(forward+backward)/2,axis=0))
    field*=dx**3
    c=np.arange(n)*dx-length/2
    x,y,z=np.meshgrid(c,c,c,indexing='ij')
    matter=0
    for i in range(count):
        radius=spec.get('radius',.9)
        weights=np.maximum(1-((x-q[i,0])**2+(y-q[i,1])**2+(z-q[i,2])**2)/radius**2,0)**3
        weights/=weights.sum()
        mass=1 if i<6 else (0 if spec['probe']=='photon' else 1e-4)
        e=np.sqrt(mass*mass+p[i]@p[i])
        matter+=np.sum(weights*(alpha*e+sum(beta[j]*p[i,j] for j in range(3))))
    return np.array([field+matter,field,matter,field+matter+state[-1]])


def main():
    read=lambda name:json.loads((OUT/name).read_text(encoding='utf8'))
    sha=lambda data:hashlib.sha256(data).hexdigest()
    checks=[]
    def check(name,error,tol=1e-12):checks.append(dict(name=name,error=float(error),tolerance=tol,passed=bool(error<=tol)))
    for name,expected in read('hashes.json').items():check('archive '+name,int(sha((OUT/name).read_bytes())!=expected),0)
    manifest=read('manifest.json')
    for name,expected in manifest['hashes'].items():
        committed=subprocess.check_output(['git','show',manifest['git_head']+':'+name],cwd=ROOT)
        check('pinned '+name,int(sha(committed)!=expected),0)
        check('working '+name,int(sha((ROOT/name).read_bytes())!=expected),0)
    rows=read('runs.json');check('evolution count',abs(len(rows)-31),0)
    errors=[]
    for row in rows:
        name=row['spec']['id'];archive=np.load(OUT/(name+'.npz'))
        energy=independent_energy(archive['final_state'],row['spec'])
        last=row['metrics'][-1];expected=np.array([last[k] for k in ('total','field','matter','ledger')])
        err=np.max(abs(energy-expected));errors.append(float(err));check('energy '+name,err)
        times=np.array(row['times']);check('times '+name,max(abs(times[0]),abs(times[-1]-2.5)),1e-12)
        check('path count '+name,abs(len(archive['q'])-len(times)),0)
        dr=max(abs(m['ledger']-row['metrics'][0]['ledger']) for m in row['metrics'])/max(1,abs(row['metrics'][0]['ledger']))
        check('energy gate reconstruction '+name,abs(dr-row['energy_drift']))
        check('ledger monotonic '+name,max(0,-np.diff([m['absorbed'] for m in row['metrics']]).min()),1e-12)
    comparisons=read('comparisons.json')
    # Reconstruct the declared comparisons from trajectory arrays, not pass labels.
    byid={r['spec']['id']:r for r in rows}
    for c in comparisons:
        if c['kind']=='rotation':
            name='rotating-eta0.08-photon' if c['n']==32 else 'rotating-space'
            q=np.load(OUT/(name+'.npz'))['q']
            qr=np.load(OUT/(f'rotated-n{c["n"]}.npz'))['q']
            axis=np.array([1.,2.,3.])/np.sqrt(14);angle=np.pi/3
            x,y,z=axis;cross=np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
            rotation=np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis,axis)+np.sin(angle)*cross
            position=np.linalg.norm(q-qr@rotation,axis=2).max()
            passed=position<=.01
        else:
            name=c['source']+'-eta0.08-photon';other=c['source']+'-'+c['kind']
            q=np.load(OUT/(name+'.npz'))['q'];qr=np.load(OUT/(other+'.npz'))['q']
            position=np.linalg.norm(q-qr[::1 if c['kind']=='domain' else 2],axis=2).max()
            fine=byid[other]['metrics'][-1]['field']
            energy_error=abs(byid[name]['metrics'][-1]['field']-fine)/max(fine,1e-10)
            check('comparison energy '+other,abs(energy_error-c['field_energy_relative']))
            passed=position<=(1e-4 if c['kind']=='domain' else .01) and energy_error<=(1e-3 if c['kind']=='domain' else .05)
        check('comparison position '+str(c),abs(position-c['position']))
        check('comparison gate '+str(c),int(bool(passed)!=c['passed']),0)
    # Explicit deposition/interpolation adjoint test at an evolved nonzero state.
    from evolution import Evolution
    model=Evolution()
    state=np.load(OUT/'rotating-eta0.08-photon.npz')['final_state'].copy()
    f,pi,q,p,alpha,a,s,b,beta,gp,gm,gc,k,j,w,e,abar,bbar,v,force,density=model.ingredients(state)
    esource=np.einsum('k,kxyz->xyz',e,w)/model.dv
    psource=np.einsum('ki,kxyz->ixyz',p,w)/model.dv
    source=np.empty_like(f)
    source[0]=model.g*alpha*(esource+np.sum(b*psource,axis=0))
    source[1:]=alpha*model.kappa*model.eta*(psource/s-model.eta**2*a*np.sum(a*psource,axis=0)/s**3)
    rng=np.random.default_rng(20260923);direction=rng.normal(size=f.shape)
    expected=np.sum(source*direction)*model.dv
    shifted=[]
    for sign in (-1,1):
        trial=state.copy();trial[:model.nf]+=sign*1e-5*direction.ravel()
        shifted.append(independent_energy(trial,byid['rotating-eta0.08-photon']['spec'])[2])
    check('source deposition interpolation adjoint',abs((shifted[1]-shifted[0])/2e-5-expected),2e-7)
    summary=dict(audit_passed=all(c['passed'] for c in checks),checks=checks,evolutions=len(rows),
                 evolution_gate_passes=sum(r['passed'] for r in rows),comparison_passes=sum(c['passed'] for c in comparisons),
                 comparison_count=len(comparisons),maximum_independent_energy_error=max(errors),
                 maximum_energy_drift=max(r['energy_drift'] for r in rows),maximum_cone_residual=max(r['max_cone_residual'] for r in rows),
                 maximum_momentum_residual=max(r['momentum_residual'] for r in rows),maximum_angular_residual=max(r['angular_residual'] for r in rows),
                 minimum_boundary_clearance=min(r['conservative_boundary_clearance'] for r in rows),
                 maximum_edge_amplitude=max(m['edge_amplitude'] for r in rows for m in r['metrics']),
                 selected=[dict(id=r['spec']['id'],field=r['metrics'][-1]['field'],matter_loss=r['metrics'][0]['matter']-r['metrics'][-1]['matter'],
                                absorbed=r['metrics'][-1]['absorbed'],angle=r['metrics'][-1]['probe_angle'],circulation=r['metrics'][-1]['circulation_z'])
                           for r in rows if r['spec']['probe']=='photon'])
    (HERE/'evolution-audit.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf8',newline='\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('checks','selected')},indent=2))
    if not summary['audit_passed']:raise AssertionError('Independent endpoint/provenance audit failed')


if __name__=='__main__':main()

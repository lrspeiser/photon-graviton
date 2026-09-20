"""SE-T full scalar coupling source; independent energy and implemented RHS checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_emitter import reconstruct
from model import ExchangeEvolution as OldModel
from emitter_model import ExchangeEvolution as NewModel

ROOT=Path(__file__).resolve().parent


def analyze(directory,cfg):
    path=directory/(cfg['name']+'.npz');state=np.load(path)['final']
    cls=NewModel if directory.name=='emitter-v1' else OldModel
    model=cls(**{k:v for k,v in cfg.items() if k not in ['name','dt']})
    f,pi,q,p,Q,P=model.unpack(state);u=.08*f[0];a=np.exp(4*u);z=np.exp(2*u);alpha=np.exp(u)
    beta=z*.04*f[1:4]/np.sqrt(1+.08**2*np.sum(f[1:4]**2,axis=0))
    current=np.zeros((3,)+u.shape);wave_current=np.zeros_like(current)
    for j in range(3):
        central=(np.roll(f,-1,j+1)-np.roll(f,1,j+1))/(2*model.dx)
        current[j]=np.sum(pi*central,axis=0);wave_current[j]=np.sum(pi[4:]*central[4:],axis=0)
    k=model.mix*np.tanh(u/.001);ku=model.mix/.001*(1-np.tanh(u/.001)**2)
    m2=.75**2*np.exp(2*model.chi*u);r=f[5]-k*f[4]
    interaction=m2*(model.chi*r*r-ku*f[4]*r)
    field=float(np.sum(2*a*np.sum(pi*pi,axis=0)-2*np.sum(beta*current,axis=0)+interaction)*model.dv)
    wave=float(np.sum(2*a*np.sum(pi[4:]**2,axis=0)-2*np.sum(beta*wave_current,axis=0)+interaction)*model.dv)
    theta=cfg.get('emitter_angle',0.);matter=0.
    for i in range(6):
        w=np.maximum(1-np.sum((model.x-q[i])**2,axis=-1)/model.radius**2,0)**3;w/=w.sum()
        bbar=np.sum(w*(np.cos(theta)*f[4]+np.sin(theta)*f[5]))
        mass=1+(P[i]**2+(Q[i]-cfg.get('emission',.4)*bbar)**2)/2
        pp=np.dot(p[i],p[i]);e=np.sqrt(mass*mass+z*pp)
        matter+=float(np.sum(w*(alpha*(e+z*pp/e)+2*np.einsum('jxyz,j->xyz',beta,p[i]))))
    total=matter+field;restoring=float(.2**2*np.sum(f[0])*model.dv/.08)
    fd,pd,*_=model.unpack(model.rhs(state))
    from_rhs=float(-np.sum(pd[0]+model.gamma*fd[0])*model.dv/.08-restoring)
    errors=[]
    for eps in [1e-7,5e-8]:
        plus=state.copy();minus=state.copy()
        model.unpack(plus)[0][0]+=eps/.08;model.unpack(minus)[0][0]-=eps/.08
        derivative=(reconstruct(plus,cfg)['total']-reconstruct(minus,cfg)['total'])/(2*eps)-restoring
        errors.append(dict(epsilon=eps,error=abs(derivative-total)))
    passed=abs(from_rhs-total)<1e-10+1e-10*abs(total) and all(e['error']<1e-7+1e-5*abs(total) for e in errors)
    return dict(campaign=directory.name,name=cfg['name'],matter_source=matter,field_source=field,
                wave_source=wave,other_field_source=field-wave,total_coupling_source=total,
                wave_fraction=wave/total if abs(total)>1e-15 else None,
                restoring_source=restoring,integrated_rhs_source=from_rhs,rhs_error=abs(from_rhs-total),
                derivative_checks=errors,passed=passed,state_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    rows=[];expected=0
    for name in ['evidence-v1','emitter-v1']:
        directory=ROOT/name;manifest=json.loads((directory/'manifest.json').read_text())
        expected+=len(manifest['configurations'])
        for cfg in manifest['configurations']:
            if (directory/(cfg['name']+'.json')).exists():rows.append(analyze(directory,cfg))
    summary=dict(passed=all(r['passed'] for r in rows),completed=len(rows),declared=expected,complete=len(rows)==expected,
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['model.py','emitter_model.py','audit_emitter.py','total-source-protocol.md']})
    (ROOT/'total-source-results.json').write_text(json.dumps(dict(summary=summary,rows=rows),indent=2)+'\n')
    lines=['# SE-T: total scalar-source accounting','',f"Completed states: {len(rows)}/{expected}. Checks pass: {summary['passed']}.",'',
           '| Campaign/case | Matter source | Other fields | X/Y waves | Total coupling source | Wave fraction |',
           '|---|---:|---:|---:|---:|---:|']
    for r in rows:
        fraction='undefined' if r['wave_fraction'] is None else f"{100*r['wave_fraction']:.5g}%"
        lines.append(f"| {r['campaign']}/{r['name']} | {r['matter_source']:.8g} | {r['other_field_source']:.8g} | {r['wave_source']:.8g} | {r['total_coupling_source']:.8g} | {fraction} |")
    lines += ['', 'These integrated scalar coupling derivatives are not mass, force or lensing.',
              'The phi restoring term is recorded separately in JSON and is not included in',
              'the table total. Field geometry and time evolution must still be solved.',
              'A large wave-sector ratio can coexist with a small share of the total source.',
              'The full energy derivative is independently reconstructed with full-grid',
              'source weights at two perturbation sizes, and checked against the integrated',
              'implemented scalar equation. All raw-state hashes and errors are preserved.',
              'Partial diagnostic success is not campaign completion or observational success.',
              'This evaluates the candidate Hamiltonian, using established differentiation;',
              'older gravity formulas are not acceptance targets.']
    (ROOT/'total-source-report.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary));assert summary['passed']


if __name__=='__main__':main()

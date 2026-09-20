"""Run the preregistered CC-1 local/stencil audit; preserve every result."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import hashlib
import subprocess
import platform
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import model as M

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def save(path, obj):
    def plain(x):
        if isinstance(x, dict): return {k:plain(v) for k,v in x.items()}
        if isinstance(x, (list,tuple,np.ndarray)): return [plain(v) for v in x]
        if isinstance(x,np.generic): return x.item()
        return x
    path.write_text(json.dumps(plain(obj),indent=2,allow_nan=False)+'\n',encoding='utf8',newline='\n')


def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def local_audit():
    rng=np.random.default_rng(20260920)
    rows=[]
    for index in range(600):
        g,eta=rng.choice([-.3,-.08,0,.08,.3],2)
        kappa=rng.choice([0,.25,.5,.9])
        f=rng.uniform(-4,4,4);p=rng.normal(size=3);wave=rng.normal(size=3)
        mass=[0,.01,1][index%3]
        args=(g,eta,kappa)
        h,v,source,hessian=M.particle(f,p,mass,*args)
        alpha,beta,db=M.geometry(f,*args)
        photon=M.particle(f,p,0,*args)
        fdv=np.zeros(3);fds=np.zeros(4);fdh=np.zeros((3,3))
        step=1e-5
        for j in range(3):
            dp=np.eye(3)[j]*step
            plus=M.particle(f,p+dp,mass,*args);minus=M.particle(f,p-dp,mass,*args)
            fdv[j]=(plus[0]-minus[0])/(2*step)
            fdh[:,j]=(plus[1]-minus[1])/(2*step)
        for j in range(4):
            df=np.eye(4)[j]*step
            fds[j]=(M.particle(f+df,p,mass,*args)[0]-M.particle(f-df,p,mass,*args)[0])/(2*step)
        block=np.eye(4)*alpha;block[0,1:]=-beta;block[1:,0]=-beta
        eigmin=np.linalg.eigvalsh(block).min()
        shift=np.dot(beta,wave);omega=.2
        matrix=np.array([[-1j*shift,alpha],[-alpha*np.dot(wave,wave)-omega**2,-1j*shift]])
        actual=np.linalg.eigvals(matrix)
        frequency=np.sqrt(alpha**2*np.dot(wave,wave)+alpha*omega**2)
        expected=np.array([-1j*(shift+frequency),-1j*(shift-frequency)])
        dispersion=min(np.max(abs(actual-expected)),np.max(abs(actual-expected[::-1])))
        rot,_=np.linalg.qr(rng.normal(size=(3,3)))
        if np.linalg.det(rot)<0: rot[:,0]*=-1
        rotated=M.particle(np.concatenate(([f[0]],rot@f[1:])),rot@p,mass,*args)
        ra,rb,_=M.geometry(np.concatenate(([f[0]],rot@f[1:])),*args)
        residuals={
            'velocity_gradient':np.max(abs(fdv-v))/max(1,np.max(abs(v))),
            'field_source_gradient':np.max(abs(fds-source))/max(1,np.max(abs(source))),
            'momentum_hessian':np.max(abs(fdh-hessian))/max(1,np.max(abs(hessian))),
            'photon_cone':abs(np.linalg.norm(photon[1]-beta)/alpha-1),
            'massive_cone_excess':max(0,np.linalg.norm(v-beta)/alpha-1),
            'momentum_convexity':max(0,-np.linalg.eigvalsh(hessian).min()),
            'energy_positivity':0 if h>0 and eigmin>0 else 1,
            'dispersion':dispersion,
            'rotation':max(abs(rotated[0]-h),np.max(abs(rotated[1]-rot@v)),abs(ra-alpha),np.max(abs(rb-rot@beta))),
            'photon_energy_scaling':np.max(abs(M.particle(f,2*p,0,*args)[1]-photon[1])),
        }
        limits={k:(2e-5 if k in ('velocity_gradient','field_source_gradient','momentum_hessian') else
                   1e-10 if k=='dispersion' else 1e-11 if k=='rotation' else 1e-12) for k in residuals}
        # With the wrong field drift sign, a ray usually misses the field's cone.
        wrong_cone=abs(np.linalg.norm(photon[1]+beta)/alpha-1)
        rows.append(dict(index=index,g=g,eta=eta,kappa=kappa,F=f,p=p,mass=mass,wavevector=wave,
                         energy=h,alpha=alpha,beta=beta,minimum_field_kinetic_eigenvalue=eigmin,
                         maximum_fourier_real_part=np.max(abs(actual.real)),wrong_sign_cone_residual=wrong_cone,
                         residuals=residuals,limits=limits,passed=all(residuals[k]<=limits[k] for k in residuals)))
    return rows


def source_audit():
    rows=[];summaries=[]
    # Gauss-Legendre integration of polynomial radial moments is independent of grid normalization.
    nodes,weights=np.polynomial.legendre.leggauss(64)
    r=(nodes+1)/2;w=weights/2
    norm=4*np.pi*np.sum(w*r*r*(1-r*r)**3)
    r2=4*np.pi*np.sum(w*r**4*(1-r*r)**3)/norm
    analytic=dict(normalization_error=abs(norm-64*np.pi/315),r2_error=abs(r2-3/11))
    for n in (32,48,64):
        dx=8/n;c=np.arange(n)*dx-4
        x=np.stack(np.meshgrid(c,c,c,indexing='ij'),axis=-1)
        for radius in (.6,.9,1.2):
            samples=[]
            for index in range(12):
                rot=M.rotation(index*np.pi/6)
                q=rot@np.array([.3,.2,.1]);center=rot@np.array([.7,-.2,.4])
                field=np.exp(-np.sum((x-center)**2,axis=-1)/2)
                weights,dw=M.spherical_weights(x,q,radius)
                sample=np.sum(weights*field);samples.append(sample)
                gradient=np.sum(dw*field[...,None],axis=(0,1,2))
                fd=[]
                for j in range(3):
                    dq=np.eye(3)[j]*1e-5
                    plus=M.spherical_weights(x,q+dq,radius)[0]
                    minus=M.spherical_weights(x,q-dq,radius)[0]
                    fd.append(np.sum((plus-minus)*field)/(2e-5))
                points=np.array([[.1,.2,.3],[-.2,.1,.4],[.8,-.1,0]])
                unrot=np.maximum(1-np.sum(points**2,axis=1)/radius**2,0)**3
                rotated=np.maximum(1-np.sum((points@rot.T)**2,axis=1)/radius**2,0)**3
                residuals=dict(weight_sum=abs(weights.sum()-1),derivative_sum=np.max(abs(dw.sum(axis=(0,1,2)))),
                               force=np.max(abs(gradient-np.array(fd))),
                               deposition_adjoint=abs(np.sum(field*(1.3*weights/dx**3))*dx**3-1.3*sample),
                               continuous_rotation=np.max(abs(unrot-rotated)))
                rows.append(dict(n=n,radius=radius,rotation_index=index,sample=sample,
                                 residuals=residuals,passed=all(v<=(2e-5 if k=='force' else 1e-12) for k,v in residuals.items())))
            spread=max(abs(np.array(samples)-samples[0]))/abs(samples[0])
            summaries.append(dict(n=n,radius=radius,reference_sample=samples[0],rotation_spread=spread,
                                  declared_gate_applies=n==64,passed=bool(spread<=.02) if n==64 else None))
            print(f'source N={n}, a={radius}: rotation spread={spread:.5g}',flush=True)
    return analytic,rows,summaries


def main():
    out=HERE/'evidence-v1';out.mkdir(exist_ok=False)
    files=sorted(list(HERE.glob('*.py'))+list(HERE.glob('*.md')))
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
         utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),numpy=np.__version__,
         hashes={p.relative_to(ROOT).as_posix():digest(p) for p in files}))
    local=local_audit();save(out/'local.json',local)
    print(f'Local samples: {sum(r["passed"] for r in local)}/600 pass',flush=True)
    analytic,source,rotation=source_audit()
    save(out/'sources.json',dict(analytic=analytic,cases=source,rotation=rotation))
    summary=dict(local_samples=len(local),local_passes=sum(r['passed'] for r in local),
                 local_max_residuals={k:max(r['residuals'][k] for r in local) for k in local[0]['residuals']},
                 source_cases=len(source),source_passes=sum(r['passed'] for r in source),
                 fine_rotation_passes=sum(r['passed'] is True for r in rotation if r['declared_gate_applies']),
                 maximum_fine_rotation_spread=max(r['rotation_spread'] for r in rotation if r['n']==64),
                 opposite_sign_control_detections=sum(r['wrong_sign_cone_residual']>1e-8 for r in local),
                 nonzero_drift_samples=sum(np.linalg.norm(r['beta'])>1e-12 for r in local),
                 max_fourier_real_part=max(r['maximum_fourier_real_part'] for r in local),
                 radial_integrals=analytic)
    summary['passed']=bool(summary['local_passes']==600 and summary['source_passes']==108 and
                           summary['fine_rotation_passes']==3 and max(analytic.values())<=1e-12)
    save(out/'summary.json',summary)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__': main()

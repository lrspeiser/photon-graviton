"""Frozen shared-drift compatibility; not a complete coupled evolution."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent


def matrices(C,cq,eps,Q,n,beta):
    a=C*C;b=float(beta@n);I=np.eye(3);Z=np.zeros((3,3));T=float(n@Q)*I-np.outer(n,Q)
    M=np.block([[-b*I,-a*eps*T.T,I,eps*b*T.T],[-a*eps*T,-b*I,eps*b*T,cq*cq*I],
                [a*I,Z,-b*I,Z],[Z,a*I,Z,-b*I]])
    S=np.block([[a*I,Z,-b*I,Z],[Z,a*I,Z,-b*I],[-b*I,Z,I,Z],[Z,-b*I,Z,cq*cq*I]])
    return M,S


def particle(p,m,C,beta):
    return np.sqrt(C)*np.sqrt(m*m+C*(p@p))+beta@p


def main():
    out=ROOT/'shifted-transfer-v1';out.mkdir(exist_ok=False);rng=np.random.default_rng(20261006);rows=[]
    for C,cq,eps,qnorm,fraction in itertools.product((.5,1.,2.),(.25,.5,1.),(0.,2.,10.),(.03,.3),(0.,.5,.9)):
        for fixture in range(4):
            n,Q,beta=rng.normal(size=(3,3));n/=np.linalg.norm(n);Q*=qnorm/np.linalg.norm(Q);beta*=fraction*cq*C/np.linalg.norm(beta)
            M,S=matrices(C,cq,eps,Q,n,beta);eigen=np.linalg.eigvals(M);speeds=-eigen.real;b=float(beta@n)
            symmetry=float(np.max(abs(S@M-M.T@S))/max(1.,np.max(abs(S@M))))
            minimum=float(np.linalg.eigvalsh(S)[0]);imaginary=float(np.max(abs(eigen.imag)))
            longitudinal=[]
            for sign in (-1.,1.):
                vector=np.concatenate((n,np.zeros(3),sign*C*n,np.zeros(3)))
                longitudinal.append(float(np.max(abs(M@vector-(-b+sign*C)*vector))))
            coverage=float(max(0.,min(speeds)-(b-C),(b+C)-max(speeds)))
            particles=[]
            for mass in (0.,1.):
                p=rng.normal(size=3);velocity=np.sqrt(C)*C*p/np.sqrt(mass*mass+C*(p@p))+beta
                errors=[]
                for delta in (1e-6,5e-7):
                    numeric=np.array([(particle(p+delta*np.eye(3)[j],mass,C,beta)-particle(p-delta*np.eye(3)[j],mass,C,beta))/(2*delta) for j in range(3)])
                    errors.append(float(np.max(abs(numeric-velocity))/max(1.,np.max(abs(velocity)))))
                speed=float(np.linalg.norm(velocity-beta))
                particles.append(dict(mass=mass,p=p.tolist(),velocity=velocity.tolist(),errors=errors,relative_drift_speed=speed,
                                      passed=bool(max(errors)<1e-7 and speed<=C+1e-12)))
            rows.append(dict(C=C,c_Q=cq,epsilon=eps,Q_norm=qnorm,beta_fraction=fraction,fixture=fixture,n=n.tolist(),Q=Q.tolist(),beta=beta.tolist(),
                             minimum_symmetrizer_eigenvalue=minimum,symmetry_error=symmetry,imaginary_error=imaginary,longitudinal_residuals=longitudinal,
                             cone_coverage_error=coverage,min_speed=float(min(speeds)),max_speed=float(max(speeds)),particles=particles,
                             passed=bool(minimum>0 and symmetry<1e-12 and imaginary<1e-10 and max(longitudinal)<1e-10 and coverage<1e-10 and all(r['passed'] for r in particles))))
    result=dict(seed=20261006,commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),cases=rows,passed=all(r['passed'] for r in rows))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(passed=result['passed'],cases=len(rows),min_S=min(r['minimum_symmetrizer_eigenvalue'] for r in rows),
                         max_symmetry=max(r['symmetry_error'] for r in rows),max_longitudinal=max(max(r['longitudinal_residuals']) for r in rows),
                         max_coverage=max(r['cone_coverage_error'] for r in rows),max_particle_error=max(max(p['errors']) for r in rows for p in r['particles']))))
    assert result['passed']


if __name__=='__main__':main()

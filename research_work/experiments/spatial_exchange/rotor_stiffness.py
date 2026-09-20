"""SE-Rot3 frozen-source vector Hessian; not full-system propagation."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.sparse.linalg import LinearOperator,eigsh
from rotor_checks import curl,weights

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'rotor-stiffness-v1';out.mkdir(exist_ok=False)
    rng=np.random.default_rng(20260930);rows=[]
    for n in (12,24):
        h=8/n;line=np.arange(n)*h-4;x=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
        theta=np.arange(6)*np.pi/3;positions=.7*np.column_stack((np.cos(theta),np.sin(theta),np.zeros(6)))
        kernels=np.array([weights(x,q)[0] for q in positions]);size=3*n**3
        Q=np.array([np.sqrt(.001),0.,0.]);T=np.dot(Q,Q)*np.eye(3)-np.outer(Q,Q)
        for coupling in (0.,2.,10.,100.):
            def apply(v):
                A=v.reshape(3,n,n,n);result=.04*A.copy()
                for axis in (1,2,3):result+=(2*A-np.roll(A,1,axis)-np.roll(A,-1,axis))/h**2
                if coupling:
                    B=curl(A,h);sampled=np.einsum('ixyz,jxyz->ij',kernels,B)
                    deposit=np.einsum('ixyz,ij->jxyz',kernels,sampled@T.T)/h**3
                    result+=coupling**2*curl(deposit,h)
                return result.ravel()
            operator=LinearOperator((size,size),matvec=apply,dtype=float)
            a=rng.normal(size=size);b=rng.normal(size=size)
            symmetry=abs(np.dot(a,apply(b))-np.dot(apply(a),b))/max(1.,np.linalg.norm(a)*np.linalg.norm(b))
            values,vectors=eigsh(operator,k=1,which='LA',tol=1e-9,v0=rng.normal(size=size),maxiter=3000)
            eigenvalue=float(values[0]);vector=vectors[:,0]
            residual=float(np.linalg.norm(apply(vector)-eigenvalue*vector)/max(1.,abs(eigenvalue)))
            bound=float(.04+12/h**2+3*coupling**2/h**2*.001*np.sum(kernels**2)/h**3)
            exact_error=abs(eigenvalue-(.04+12/h**2)) if coupling==0 else None
            passed=residual<1e-7 and symmetry<1e-10 and .04-1e-7<=eigenvalue<=bound+1e-7 and (exact_error is None or exact_error<1e-7)
            row=dict(n=n,coupling=coupling,largest_eigenvalue=eigenvalue,upper_bound=bound,residual=residual,symmetry_error=float(symmetry),checkerboard_error=exact_error,
                     conservative_half_limit=float(np.sqrt(2)/np.sqrt(bound)),passed=bool(passed))
            rows.append(row);np.savez_compressed(out/f'n{n}-coupling{coupling}.npz',eigenvector=vector)
            print(n,coupling,'lambda',eigenvalue,'half timestep',row['conservative_half_limit'],'passed',passed,flush=True)
    result=dict(passed=all(r['passed'] for r in rows),count=len(rows),
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['rotor_stiffness.py','rotor_checks.py','rotor-stiffness-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=result,rows=rows),indent=2)+'\n')
    print(json.dumps(result));assert result['passed']


if __name__=='__main__':main()

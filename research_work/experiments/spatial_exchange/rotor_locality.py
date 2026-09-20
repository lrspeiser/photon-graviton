"""SE-Rot4: separated-cell response of the finite rotor sampling kernel."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from rotor_checks import curl,weights

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'rotor-locality-v1';out.mkdir(exist_ok=False);rows=[]
    for n in (24,48):
        h=8/n;line=np.arange(n)*h-4;xyz=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
        W,_=weights(xyz,np.zeros(3));derivative=(np.roll(W,-1,0)-np.roll(W,1,0))/(2*h)
        center=n//2;offset=round((2/3)/h);source=(center+offset,center,center);target=(center-offset,center,center)
        for coupling in (0.,2.,10.):
            def response(start,end):
                A=np.zeros((3,n,n,n));A[(2,)+start]=1
                B=np.sum(W*curl(A,h),axis=(1,2,3));T=np.diag([0.,.001,.001])
                added=coupling**2*curl(W*(T@B)[:,None,None,None]/h**3,h)
                local=.04*A.copy()
                for axis in (1,2,3):local+=(2*A-np.roll(A,1,axis)-np.roll(A,-1,axis))/h**2
                return float(added[(2,)+end]),float(local[(2,)+end])
            actual,local=response(source,target);reverse,_=response(target,source)
            analytic=float(coupling**2*.001*derivative[source]*derivative[target]/h**3)
            error=abs(actual-analytic);reciprocity=abs(actual-reverse)
            passed=error<=1e-12+1e-10*abs(actual) and reciprocity<1e-12 and local==0 and (coupling!=0 or actual==0)
            rows.append(dict(n=n,coupling=coupling,separation=2*offset*h,stencil_response=local,rotor_response=actual,
                             response_per_cell_volume=actual/h**3,factorized_response=analytic,error=error,
                             reciprocity_error=reciprocity,nonlocal_entry=actual!=0,passed=passed))
    scaling=[]
    for n in (24,48):
        selected=[r for r in rows if r['n']==n and r['coupling']>0]
        error=abs(selected[1]['rotor_response']/selected[0]['rotor_response']-25)
        scaling.append(dict(n=n,error=error,passed=error<1e-10))
    summary=dict(passed=all(r['passed'] for r in rows+scaling),count=len(rows),
                 finite_kernel_nonlocality_detected=all(r['nonlocal_entry'] for r in rows if r['coupling']>0),
                 commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                 sources={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['rotor_locality.py','rotor_checks.py','rotor-locality-protocol.md']})
    (out/'results.json').write_text(json.dumps(dict(summary=summary,rows=rows,scaling=scaling),indent=2)+'\n')
    print(json.dumps(summary));print(json.dumps(rows));assert summary['passed']


if __name__=='__main__':main()

"""Resolve failing path comparisons with independent per-orbit integrations."""
import argparse
from pathlib import Path
import json
import hashlib
import time
import numpy as np
from scipy.integrate import solve_ivp
import run as parent

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('model',choices=['ordinary','axisymmetric_extra','full'])
args=parser.parse_args()
field=parent.Field(args.model)
original=np.load(parent.CACHE/f'{args.model}.npz')['trajectory']
print('Rechecking previous tolerance to identify failing paths',flush=True)
low,_=parent.integrate(field,2e-11)
dx=np.max(np.linalg.norm(original[:,:,:3]-low[:,:,:3],axis=2),axis=0)
dv=np.max(np.linalg.norm(original[:,:,3:]-low[:,:,3:],axis=2),axis=0)
failed=np.flatnonzero((dx>=1e-4)|(dv>=.01))
print('Failing paths',failed.tolist(),dx.tolist(),dv.tolist(),flush=True)
retained=original.copy();details=[]
for j in failed:
    outputs=[];attempts=[]
    for tol in [2e-13,2.3e-14]:
        def rhs(t,y):
            s=y.reshape(1,6);x,p=s[:,:3],s[:,3:]
            a=field.evaluate(x)[1]
            return np.c_[p-parent.OMEGA*parent.cross(x),a-parent.OMEGA*parent.cross(p)].ravel()
        started=time.monotonic()
        sol=solve_ivp(rhs,[0,.25],parent.Y0[j],t_eval=parent.TIMES,method='DOP853',rtol=tol,atol=tol*.01)
        assert sol.success
        outputs.append(sol.y.T)
        attempts.append(dict(rtol=tol,nfev=sol.nfev,seconds=time.monotonic()-started))
        print('Individual',int(j),attempts[-1],flush=True)
    errx=float(np.max(np.linalg.norm(outputs[1][:,:3]-outputs[0][:,:3],axis=1)))
    errv=float(np.max(np.linalg.norm(outputs[1][:,3:]-outputs[0][:,3:],axis=1)))
    passed=bool(errx<1e-4 and errv<.01)
    details.append(dict(probe=int(j),attempts=attempts,position_difference_kpc=errx,velocity_difference_kms=errv,passed=passed))
    retained[:,j]=outputs[-1]
np.savez_compressed(parent.CACHE/f'{args.model}-individual-refinement.npz',times=parent.TIMES,trajectory=retained)
summary=parent.summary(field,retained)
out=dict(model=args.model,scope='Numerical path refinement only; identical physical initial states and field.',
         original_batch_refinement_dx=dx.tolist(),original_batch_refinement_dv=dv.tolist(),
         failing_paths=failed.tolist(),individual_checks=details,
         all_selected_path_checks_pass=all(d['passed'] for d in details),
         maximum_scaled_Jacobi_drift=max(s['Jacobi_drift_over_220_squared'] for s in summary),
         probes=summary,
         original_trajectory_sha256=hashlib.sha256((parent.CACHE/f'{args.model}.npz').read_bytes()).hexdigest(),
         code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         no_data_or_parameters_fitted=True)
(HERE/f'{args.model}-individual-refinement.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(details,indent=2),flush=True)

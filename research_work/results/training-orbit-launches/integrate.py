"""Integrate provisional training-derived orbit basis; retain every failure."""
from pathlib import Path
import hashlib
import importlib.util
import json
import time
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DRIVER=HERE.parent/'full-bar-orbits/run.py'
spec=importlib.util.spec_from_file_location('verified_orbit',DRIVER)
orbit=importlib.util.module_from_spec(spec);spec.loader.exec_module(orbit)
CACHE=ROOT/'research_work/data-cache/training-orbit-launches'
SOURCE=CACHE/'launches.parquet'
seeds=pd.read_parquet(SOURCE)
assert len(seeds)==72 and (seeds.holdout_role=='training').all()
initial=seeds[['x_kpc','y_kpc','z_kpc','vx_kms','vy_kms','vz_kms']].to_numpy(float)
TIMES=np.linspace(0,.25,501)
field=orbit.Field('full')


def save(rows,complete=False):
    data=dict(scope='Provisional training-derived full-field orbit library; no population fit.',
              model='full_bar_empirical_response',seeds=72,completed=complete,
              duration_kpc_per_kms=.25,sample_times=501,
              tolerances=[2e-9,2e-11,2e-13],position_gate_kpc=1e-4,velocity_gate_kms=.01,
              Jacobi_gate_over_220_squared=1e-5,rows=rows,
              source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              driver_sha256=hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
              input_field_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in field.paths},
              holdouts_opened=False,likelihood_fitted=False,association_uncertainty_resolved=False)
    name='integration-results.json' if complete else 'integration-progress.json'
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')


def main():
    rows=[]
    trajectories=np.full((len(TIMES),len(seeds),6),np.nan)
    for j,y0 in enumerate(initial):
        started=time.monotonic();attempts=[];previous=None;retained=None;error=None
        for tol in [2e-9,2e-11,2e-13]:
            def rhs(t,y):
                state=y.reshape(1,6);x,p=state[:,:3],state[:,3:]
                a=field.evaluate(x)[1]
                return np.c_[p-orbit.OMEGA*orbit.cross(x),a-orbit.OMEGA*orbit.cross(p)].ravel()
            try:
                sol=solve_ivp(rhs,[0,.25],y0,t_eval=TIMES,method='DOP853',rtol=tol,atol=tol*.01)
                if not sol.success:raise RuntimeError(sol.message)
            except (ValueError,RuntimeError) as exc:
                error=str(exc);break
            y=sol.y.T
            item=dict(rtol=tol,nfev=sol.nfev)
            if previous is not None:
                dx=float(np.max(np.linalg.norm(y[:,:3]-previous[:,:3],axis=1)))
                dv=float(np.max(np.linalg.norm(y[:,3:]-previous[:,3:],axis=1)))
                item.update(position_difference_kpc=dx,velocity_difference_kms=dv,passed=bool(dx<1e-4 and dv<.01))
            attempts.append(item);retained=y
            if item.get('passed',False):break
            previous=y
        row=dict(seed_index=j,source_id=str(seeds.iloc[j].source_id),R_stratum=int(seeds.iloc[j].R_stratum),
                 height_stratum=int(seeds.iloc[j].height_stratum),distance_disagreement_flag=bool(seeds.iloc[j].parallax_interval_disjoint_5),
                 attempts=attempts,error=error,seconds=time.monotonic()-started)
        if error is None and retained is not None:
            trajectories[:,j]=retained
            # Evaluate in small chunks to avoid large harmonic work arrays.
            pot=np.concatenate([field.evaluate(q[:,:3])[0] for q in np.array_split(retained,32)])
            x,p=retained[:,:3],retained[:,3:]
            J=.5*np.sum(p*p,axis=1)+pot-orbit.OMEGA*(x[:,0]*p[:,1]-x[:,1]*p[:,0])
            drift=float(np.max(abs(J-J[0]))/220**2)
            row.update(Jacobi_drift_over_220_squared=drift,
                       integration_gate_pass=attempts[-1].get('passed',False),
                       numerical_checks_pass=bool(attempts[-1].get('passed',False) and drift<1e-5),
                       radius_min_kpc=float(np.linalg.norm(x,axis=1).min()),radius_max_kpc=float(np.linalg.norm(x,axis=1).max()),
                       max_absolute_height_kpc=float(abs(x[:,2]).max()))
        else:row['numerical_checks_pass']=False
        rows.append(row);save(rows)
        print(j+1,'/72',row['numerical_checks_pass'],round(row['seconds'],1),'s',error or '',flush=True)
    path=CACHE/'full-field-orbits.npz'
    np.savez_compressed(path,times=TIMES,trajectories=trajectories,source_id=seeds.source_id.to_numpy())
    save(rows,True)
    print('Passed',sum(r['numerical_checks_pass'] for r in rows),'of72',flush=True)


if __name__=='__main__':main()

"""Matched ordinary/extra-potential duration diagnostics on fixed training starts."""
from pathlib import Path
import importlib.util
import json
import sys
import time
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

exp=load('paired_expansion',HERE.parent/'orbit-library-expansion/run.py')
old=load('paired_old_duration',HERE.parent/'orbit-duration-audit/run.py')
orbit=old.orbit
CACHE=ROOT/'research_work/data-cache'
OUT=CACHE/'paired-orbit-duration';OUT.mkdir(parents=True,exist_ok=True)
INDICES=list(range(0,36,4))
TIMES=np.linspace(0,1,2001)
SOURCE=CACHE/'orbit-library-expansion/launches.parquet'
digest=exp.digest

def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def main(kind):
    assert kind in ['full','ordinary']
    seeds=pd.read_parquet(SOURCE)
    selection=json.loads((HERE.parent/'orbit-library-expansion/selection.json').read_text())
    assert digest(SOURCE)==selection['launches_sha256']
    assert len(seeds)==36 and (seeds.holdout_role=='training').all()
    selected=seeds.iloc[INDICES]
    assert len(selected[['R_stratum','height_stratum']].drop_duplicates())==9
    assert not selected.parallax_interval_disjoint_5.any()
    field=orbit.Field(kind)
    paths=[Path(__file__),HERE/'protocol.md',SOURCE,HERE.parent/'full-bar-orbits/run.py',
           HERE.parent/'orbit-duration-audit/run.py',HERE.parent/'orbit-library-expansion/run.py',*field.paths]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in paths}
    rows=[]
    for index in INDICES:
        started=time.monotonic();attempts=[];previous=None;retained=None;error=None
        y0=seeds.iloc[index][exp.NAMES].to_numpy(float)
        def rhs(t,y):
            x,p=y[None,:3],y[None,3:]
            a=field.evaluate(x)[1]
            return np.c_[p-orbit.OMEGA*orbit.cross(x),a-orbit.OMEGA*orbit.cross(p)].ravel()
        for tol in [2e-9,2e-11,2e-13]:
            try:
                sol=solve_ivp(rhs,[0,1],y0,method='DOP853',t_eval=TIMES,rtol=tol,atol=tol*.01)
                if not sol.success:raise RuntimeError(sol.message)
            except (ValueError,RuntimeError) as exc:
                error=str(exc);break
            retained=sol.y.T;attempt=dict(rtol=tol,nfev=sol.nfev)
            if previous is not None:
                dx=float(np.linalg.norm(retained[:,:3]-previous[:,:3],axis=1).max())
                dv=float(np.linalg.norm(retained[:,3:]-previous[:,3:],axis=1).max())
                attempt.update(position_difference_kpc=dx,velocity_difference_kms=dv,passed=dx<1e-4 and dv<.01)
            attempts.append(attempt)
            print(kind,index,'tolerance',tol,'pass',attempt.get('passed'),flush=True)
            if attempt.get('passed',False):break
            previous=retained
        row=dict(seed_index=index,source_id=str(seeds.iloc[index].source_id),
            R_stratum=int(seeds.iloc[index].R_stratum),height_stratum=int(seeds.iloc[index].height_stratum),
            attempts=attempts,error=error,numerical_pass=False,complete_trajectory=retained is not None)
        if retained is not None:
            path=OUT/f'{kind}-{index}.npz'
            np.savez_compressed(path,times=TIMES,trajectory=retained,source_id=seeds.iloc[index].source_id)
            row['trajectory_sha256']=digest(path)
            pot=np.concatenate([field.evaluate(q[:,:3])[0] for q in np.array_split(retained,128)])
            x,p=retained[:,:3],retained[:,3:]
            J=.5*np.sum(p*p,axis=1)+pot-orbit.OMEGA*(x[:,0]*p[:,1]-x[:,1]*p[:,0])
            drift=float(np.max(abs(J-J[0]))/220**2)
            row.update(Jacobi_drift_over_220_squared=drift,
                numerical_pass=bool(error is None and attempts[-1].get('passed',False) and drift<1e-5),
                radius_min_kpc=float(np.linalg.norm(x,axis=1).min()),radius_max_kpc=float(np.linalg.norm(x,axis=1).max()),
                max_abs_height_kpc=float(abs(x[:,2]).max()))
            angle=np.unwrap(np.arctan2(x[:,1],x[:,0]))
            row['bar_angle_range_turns']=float(np.ptp(angle)/(2*np.pi))
            row['occupation_windows']=[]
            for stop in [500,1000,2000]:
                a,b=retained[1:stop//2+1],retained[stop//2+1:stop+1]
                row['occupation_windows'].append(dict(duration_kpc_per_kms=float(TIMES[stop]),
                    spatial_TV=old.tv(a,b),R_z_TV=old.tv(a,b,False),half_sampling_spatial_TV=old.tv(a[1::2],b[1::2])))
        row['seconds']=time.monotonic()-started;rows.append(row)
        save(f'{kind}-integration.json',dict(kind=kind,rows=rows,indices=INDICES,completed=len(rows)==9,
            duration_kpc_per_kms=1.,source_hashes=hashes,holdouts_opened=False,stellar_likelihood=False))
        print(kind,index,'COMPLETE',row['numerical_pass'],round(row['seconds'],1),'s',error or '',flush=True)
    for p in paths:assert digest(p)==hashes[str(p.relative_to(ROOT))]

if __name__=='__main__':main(sys.argv[1])

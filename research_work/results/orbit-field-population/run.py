"""Numerical field versus solver sensitivity of one retained orbit population."""
from pathlib import Path
import importlib.util
import json
import time
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('field_phase',HERE.parent/'orbit-phase-duration/run.py')
phase=importlib.util.module_from_spec(spec);spec.loader.exec_module(phase)
spec=importlib.util.spec_from_file_location('field_population',HERE.parent/'orbit-phase-duration/population_resolution.py')
pop=importlib.util.module_from_spec(spec);spec.loader.exec_module(pop)
CACHE=phase.ROOT/'research_work/data-cache/orbit-field-population';CACHE.mkdir(parents=True,exist_ok=True)
scales=pop.scales

def moments(y):
    values=phase.smooth.prior.features(y);R,z=values[:,:2].T
    labels=np.searchsorted([0,3.5,5,9,np.inf],R,side='right')-1
    out=[]
    for r in range(4):
        for high in [False,True]:
            sel=(labels==r)&((abs(z)>=.2)==high);v=values[sel,3:]
            out.append(dict(radial_bin=r,offplane=high,count=int(sel.sum()),
                mean_velocity_kms=v.mean(axis=0).tolist() if len(v) else None,
                dispersion_kms=v.std(axis=0).tolist() if len(v) else None))
    return out

def main():
    olddir=HERE.parent/'orbit-phase-duration';base=json.loads((olddir/'results.json').read_text())
    previous=json.loads((olddir/'population-resolution.json').read_text())
    assert base['completed'] and previous['completed']
    for data in [base,previous]:
        for p,sha in data['source_hashes'].items():assert phase.digest(phase.ROOT/p)==sha
    hashes={str(p.relative_to(phase.ROOT)):phase.digest(p) for p in [Path(__file__),HERE/'protocol.md',olddir/'results.json',olddir/'population-resolution.json']}
    highpath=phase.CACHE/'full-12.npz';assert phase.digest(highpath)==base['models']['full']['trajectory_sha256']
    hashes[str(highpath.relative_to(phase.ROOT))]=phase.digest(highpath)
    with np.load(highpath) as f:high=f['trajectory']
    lower=phase.orbit.Field('lower_field');upper=phase.orbit.Field('full');ordinary=phase.orbit.Field('ordinary')
    for p in lower.paths+upper.paths+ordinary.paths:hashes[str(p.relative_to(phase.ROOT))]=phase.digest(p)
    def rhs(t,y):
        x,p=y[None,:3],y[None,3:];a=lower.evaluate(x)[1]
        return np.c_[p-phase.orbit.OMEGA*phase.orbit.cross(x),a-phase.orbit.OMEGA*phase.orbit.cross(p)].ravel()
    paths=[];attempts=[]
    for tol in [2e-11,2e-13]:
        started=time.monotonic()
        sol=solve_ivp(rhs,[0,4],high[0],method='DOP853',t_eval=phase.TIMES,rtol=tol,atol=tol*.01)
        assert sol.success,sol.message;y=sol.y.T;paths.append(y)
        path=CACHE/f'lower-12-{tol}.npz';np.savez_compressed(path,times=phase.TIMES,trajectory=y)
        phi=np.concatenate([lower.evaluate(q[:,:3])[0] for q in np.array_split(y,128)])
        J=.5*np.sum(y[:,3:]**2,axis=1)+phi-phase.orbit.OMEGA*(y[:,0]*y[:,4]-y[:,1]*y[:,3])
        attempts.append(dict(rtol=tol,nfev=sol.nfev,seconds=time.monotonic()-started,trajectory_sha256=phase.digest(path),
            Jacobi_drift_over_220_squared=float(np.max(abs(J-J[0]))/220**2)))
        print('lower field',attempts[-1],flush=True)
    low=paths[-1];coarse=paths[0]
    dx=float(np.linalg.norm(low[:,:3]-coarse[:,:3],axis=1).max());dv=float(np.linalg.norm(low[:,3:]-coarse[:,3:],axis=1).max())
    force=[]
    for trajectory in [low,high]:
        errors=[]
        for q in np.array_split(trajectory[:,:3],128):
            bary=ordinary.evaluate(q)[1];a=upper.evaluate(q)[1]-bary;b=lower.evaluate(q)[1]-bary
            errors.extend(np.linalg.norm(a-b,axis=1)/np.linalg.norm(a,axis=1))
        force.append(dict(maximum=float(max(errors)),median=float(np.median(errors))))
    oldrow=next(r for r in previous['rows'] if r['kind']=='full')
    rows=[]
    for name,indices in [(f'first_{d}',np.arange(d,2000*d+1,d)) for d in [1,2,4]]+[('last_unit',np.arange(6001,8001))]:
        a,b,c=low[indices],high[indices],coarse[indices]
        fieldgrams=phase.smooth.kernel_grams(np.concatenate([a,b])[None,:,:],scales)
        solvergrams=phase.smooth.kernel_grams(np.concatenate([a,c])[None,:,:],scales)
        earlier=next(w for w in oldrow['windows'] if w['window']==name)
        cells=[]
        for l,h in zip(moments(a),moments(b)):
            eligible=min(l['count'],h['count'])>=50
            cells.append(dict(lower=l,upper=h,comparison_eligible=eligible,
                lower_minus_upper_mean_velocity_kms=(np.array(l['mean_velocity_kms'])-h['mean_velocity_kms']).tolist() if eligible else None,
                lower_minus_upper_dispersion_kms=(np.array(l['dispersion_kms'])-h['dispersion_kms']).tolist() if eligible else None))
        row=dict(window=name,field_occupation_TV=pop.tv(a,b),lower_solver_occupation_TV=pop.tv(a,c),
                 cells=cells,smooth=[])
        for (sx,sp),G,H,old in zip(scales,fieldgrams,solvergrams,earlier['smooth']):
            row['smooth'].append(dict(position_scale_kpc=sx,velocity_scale_kms=sp,
                field_occupation_MMD=float(np.sqrt(max(0,G[0,0]))),
                lower_solver_occupation_MMD=float(np.sqrt(max(0,H[0,0]))),
                upper_solver_occupation_MMD=old['solver_occupation_MMD']))
        rows.append(row);print(name,'field TV',row['field_occupation_TV'],'MMD',row['smooth'],flush=True)
    for p,sha in hashes.items():assert phase.digest(phase.ROOT/p)==sha
    result=dict(scope='Field-grid and solver sensitivity of one targeted candidate orbit; not a stellar likelihood.',
        attempts=attempts,lower_field_position_difference_kpc=dx,lower_field_velocity_difference_kms=dv,
        lower_field_strict_path_pass=bool(dx<1e-4 and dv<.01 and attempts[-1]['Jacobi_drift_over_220_squared']<1e-5),
        extra_force_refinement_on_lower_and_upper_paths=force,rows=rows,source_hashes=hashes,
        velocity_order=['radial','rotational','vertical'],holdouts_opened=False,physical_parameters_changed=False,completed=True)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

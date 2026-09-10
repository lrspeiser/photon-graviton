"""Paired long integration of a retained, diagnostically selected training orbit."""
from pathlib import Path
import importlib.util
import hashlib
import json
import time
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
paired=load('phase_paired',HERE.parent/'paired-orbit-duration/run.py')
smooth=load('phase_smooth',HERE.parent/'orbit-stability-resolution/run.py')
orbit=paired.orbit
CACHE=ROOT/'research_work/data-cache/orbit-phase-duration';CACHE.mkdir(parents=True,exist_ok=True)
TIMES=np.linspace(0,4,8001)
def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(obj):
    (HERE/'results.json').write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def describe(states):
    theta=np.unwrap(np.arctan2(states[:,1],states[:,0]));dt=np.diff(theta)
    return dict(radius_range_kpc=[float(np.hypot(states[:,0],states[:,1]).min()),float(np.hypot(states[:,0],states[:,1]).max())],
        max_abs_height_kpc=float(abs(states[:,2]).max()),bar_angle_span_turns=float(np.ptp(theta)/(2*np.pi)),
        net_bar_angle_turns=float((theta[-1]-theta[0])/(2*np.pi)),sampled_angular_reversals=int(np.sum(dt[1:]*dt[:-1]<0)))

def main():
    proof=json.loads((HERE.parent/'orbit-stability-resolution/results.json').read_text())
    for p,sha in proof['source_hashes'].items():assert digest(ROOT/p)==sha
    files=[Path(__file__),HERE/'protocol.md',HERE.parent/'orbit-stability-resolution/run.py',
           HERE.parent/'orbit-stability-resolution/results.json',HERE.parent/'paired-orbit-duration/run.py',
           HERE.parent/'full-bar-orbits/run.py']
    hashes={str(p.relative_to(ROOT)):digest(p) for p in files}
    result=dict(scope='Targeted long-duration numerical phase coverage; no observational likelihood.',seed_index=12,
        duration_kpc_per_kms=4.,models={},source_hashes=hashes,holdouts_opened=False,physical_parameters_changed=False,
        actual_galactic_history_simulated=False,completed=False)
    for kind in ['ordinary','full']:
        oldpath=paired.OUT/f'{kind}-12.npz';hashes[str(oldpath.relative_to(ROOT))]=digest(oldpath)
        with np.load(oldpath) as f:old=f['trajectory'];y0=old[0];source_id=str(f['source_id'])
        field=orbit.Field(kind)
        for p in field.paths:hashes[str(p.relative_to(ROOT))]=digest(p)
        def rhs(t,y):
            x,p=y[None,:3],y[None,3:];a=field.evaluate(x)[1]
            return np.c_[p-orbit.OMEGA*orbit.cross(x),a-orbit.OMEGA*orbit.cross(p)].ravel()
        previous=None;states=None;row=dict(source_id=source_id,attempts=[],numerical_pass=False,error=None)
        result['models'][kind]=row
        for tol in [2e-9,2e-11,2e-13]:
            start=time.monotonic()
            try:
                sol=solve_ivp(rhs,[0,4],y0,t_eval=TIMES,method='DOP853',rtol=tol,atol=tol*.01)
                if not sol.success:raise RuntimeError(sol.message)
            except (ValueError,RuntimeError) as e:
                row['error']=str(e);break
            states=sol.y.T;attempt=dict(rtol=tol,nfev=sol.nfev,seconds=time.monotonic()-start)
            if previous is not None:
                dx=float(np.linalg.norm(states[:,:3]-previous[:,:3],axis=1).max())
                dv=float(np.linalg.norm(states[:,3:]-previous[:,3:],axis=1).max())
                attempt.update(position_difference_kpc=dx,velocity_difference_kms=dv,passed=dx<1e-4 and dv<.01)
            row['attempts'].append(attempt);save(result);print(kind,attempt,flush=True)
            if attempt.get('passed',False):break
            previous=states
        if states is None:continue
        path=CACHE/f'{kind}-12.npz';np.savez_compressed(path,times=TIMES,trajectory=states)
        row['trajectory_sha256']=digest(path)
        phi=np.concatenate([field.evaluate(q[:,:3])[0] for q in np.array_split(states,128)])
        J=.5*np.sum(states[:,3:]**2,axis=1)+phi-orbit.OMEGA*(states[:,0]*states[:,4]-states[:,1]*states[:,3])
        drift=float(np.max(abs(J-J[0]))/220**2)
        dx=float(np.linalg.norm(states[:2001,:3]-old[:,:3],axis=1).max())
        dv=float(np.linalg.norm(states[:2001,3:]-old[:,3:],axis=1).max())
        row.update(Jacobi_drift_over_220_squared=drift,previous_first_unit_position_difference_kpc=dx,
            previous_first_unit_velocity_difference_kms=dv,full_path_description=describe(states))
        gate=bool(row['error'] is None and row['attempts'][-1].get('passed',False) and drift<1e-5 and dx<1e-4 and dv<.01)
        if kind=='full':
            lower=orbit.Field('lower_field');ordinary=orbit.Field('ordinary')
            for p in lower.paths+ordinary.paths:hashes[str(p.relative_to(ROOT))]=digest(p)
            errors=[]
            for q in np.array_split(states[:,:3],128):
                base=ordinary.evaluate(q)[1];a=field.evaluate(q)[1]-base;b=lower.evaluate(q)[1]-base
                errors.extend(np.linalg.norm(a-b,axis=1)/np.linalg.norm(a,axis=1))
            row['max_extra_force_refinement_fraction']=float(max(errors));gate=gate and max(errors)<.01
        row['numerical_pass']=bool(gate);row['windows']=[]
        if gate:
            windows=[(f'first_{duration}',states[duration:2000*duration+1:duration]) for duration in [1,2,4]]
            windows.append(('last_unit',states[6001:8001]))
            for name,s in windows:
                assert s.shape==(2000,6)
                grams=smooth.kernel_grams(s[None,:,:],[(.5,50.),(2.,100.)])
                hist=smooth.tvs(s[None,:,:],{'single':np.ones(1)})
                entry=dict(window=name,description=describe(s),histogram_TV={k:v['single'] for k,v in hist.items()},
                    smooth=[dict(position_scale_kpc=sx,velocity_scale_kms=sp,chronological_MMD=float(np.sqrt(max(0,g[0,0]))),
                        interleaved_MMD=float(np.sqrt(max(0,g[1,1])))) for (sx,sp),g in zip([(.5,50.),(2.,100.)],grams)])
                row['windows'].append(entry)
        save(result);print(kind,'GATE',gate,'windows',row['windows'],flush=True)
    for p,sha in hashes.items():assert digest(ROOT/p)==sha
    result['completed']=True;save(result)

if __name__=='__main__':main()

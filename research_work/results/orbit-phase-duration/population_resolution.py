"""Post-failure comparison of solver occupations, with no replacement pass gate."""
from pathlib import Path
import importlib.util
import json
import time
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('phase_extension',HERE/'run.py')
ext=importlib.util.module_from_spec(spec);spec.loader.exec_module(ext)
scales=[(.5,50.),(2.,100.)]

def tv(a,b):
    values=np.vstack([ext.smooth.prior.features(a),ext.smooth.prior.features(b)])
    bins=np.stack([np.searchsorted(ext.smooth.prior.EDGES[j],values[:,j],side='right')-1 for j in range(6)],axis=1)
    labels=np.ravel_multi_index(bins.T,tuple(len(e)-1 for e in ext.smooth.prior.EDGES))
    _,inv=np.unique(labels,return_inverse=True);n=int(inv.max()+1)
    return float(.5*np.abs(np.bincount(inv[:len(a)],minlength=n)/len(a)-np.bincount(inv[len(a):],minlength=n)/len(b)).sum())

def main():
    base=json.loads((HERE/'results.json').read_text());assert base['completed']
    for p,sha in base['source_hashes'].items():assert ext.digest(ext.ROOT/p)==sha
    hashes={str(p.relative_to(ext.ROOT)):ext.digest(p) for p in [Path(__file__),HERE/'results.json',HERE/'population-resolution-protocol.md']}
    rows=[]
    for kind,row in base['models'].items():
        if row['numerical_pass']:continue
        assert row['attempts'][-1]['rtol']==2e-13 and row['attempts'][-2]['rtol']==2e-11
        path=ext.CACHE/f'{kind}-12.npz';assert ext.digest(path)==row['trajectory_sha256']
        hashes[str(path.relative_to(ext.ROOT))]=ext.digest(path)
        with np.load(path) as f:fine=f['trajectory']
        field=ext.orbit.Field(kind)
        def rhs(t,y):
            x,p=y[None,:3],y[None,3:];a=field.evaluate(x)[1]
            return np.c_[p-ext.orbit.OMEGA*ext.orbit.cross(x),a-ext.orbit.OMEGA*ext.orbit.cross(p)].ravel()
        started=time.monotonic()
        sol=solve_ivp(rhs,[0,4],fine[0],method='DOP853',t_eval=ext.TIMES,rtol=2e-11,atol=2e-13)
        assert sol.success;coarse=sol.y.T
        dx=float(np.linalg.norm(coarse[:,:3]-fine[:,:3],axis=1).max())
        dv=float(np.linalg.norm(coarse[:,3:]-fine[:,3:],axis=1).max())
        np.testing.assert_allclose([dx,dv],[row['attempts'][-1]['position_difference_kpc'],row['attempts'][-1]['velocity_difference_kms']],rtol=1e-7,atol=1e-9)
        path=ext.CACHE/f'{kind}-12-rtol-2e-11.npz';np.savez_compressed(path,times=ext.TIMES,trajectory=coarse)
        entry=dict(kind=kind,rtol=2e-11,seconds=time.monotonic()-started,nfev=sol.nfev,trajectory_sha256=ext.digest(path),windows=[])
        for name,indices in [(f'first_{d}',np.arange(d,2000*d+1,d)) for d in [1,2,4]]+[('last_unit',np.arange(6001,8001))]:
            a,b=coarse[indices],fine[indices]
            # First Gram block compares entire solver occupations; individual
            # chronological Grams compare early and late within each solver.
            between=ext.smooth.kernel_grams(np.concatenate([a,b])[None,:,:],scales)
            within=ext.smooth.kernel_grams(np.stack([a,b]),scales)
            window=dict(window=name,max_position_difference_kpc=float(np.linalg.norm(a[:,:3]-b[:,:3],axis=1).max()),
                max_velocity_difference_kms=float(np.linalg.norm(a[:,3:]-b[:,3:],axis=1).max()),solver_occupation_TV=tv(a,b),smooth=[])
            for (sx,sp),G,H in zip(scales,between,within):
                window['smooth'].append(dict(position_scale_kpc=sx,velocity_scale_kms=sp,
                    solver_occupation_MMD=float(np.sqrt(max(0,G[0,0]))),
                    coarse_early_late_MMD=float(np.sqrt(max(0,H[0,0]))),fine_early_late_MMD=float(np.sqrt(max(0,H[1,1])))))
            entry['windows'].append(window)
        rows.append(entry);print(kind,entry['windows'],flush=True)
    for p,sha in hashes.items():assert ext.digest(ext.ROOT/p)==sha
    result=dict(scope='Post-failure population sensitivity to solver accuracy, not a replacement gate.',rows=rows,source_hashes=hashes,
        primary_path_gate_changed=False,observational_significance_claimed=False,holdouts_opened=False,completed=True)
    (HERE/'population-resolution.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

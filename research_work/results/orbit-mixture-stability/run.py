"""Certified convex common-weight stability diagnostic, not an observed fit."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.optimize import linprog
from scipy import sparse

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache/paired-orbit-duration'
SOURCE=HERE.parent/'paired-orbit-duration'

def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def features(states):
    x,y,z,px,py,pz=states.T;R=np.hypot(x,y)
    assert np.all(R>0)
    vr=(x*px+y*py)/R;vp=(x*py-y*px)/R
    return np.c_[R,z,np.arctan2(y,x),vr,vp,pz]

EDGES=[np.array(x,float) for x in [
    [0,.5,3.5,5,9,20,40,80,np.inf],[-np.inf,-1.5,-.5,-.2,0,.2,.5,1.5,np.inf],
    np.linspace(-np.pi,np.pi,9),[-np.inf,-150,-100,-50,0,50,100,150,np.inf],
    [-np.inf,0,100,150,200,250,300,np.inf],[-np.inf,-150,-100,-50,0,50,100,150,np.inf]]]

def matrices(paths,dimensions,step):
    values=np.stack([features(y) for y in paths],axis=1)
    indices=np.stack([np.searchsorted(EDGES[j],values[:,:,j],side='right')-1 for j in range(dimensions)],axis=2)
    shape=tuple(len(e)-1 for e in EDGES[:dimensions])
    labels=np.ravel_multi_index(indices.reshape(-1,dimensions).T,shape).reshape(values.shape[:2])
    unique,inverse=np.unique(labels,return_inverse=True);labels=inverse.reshape(labels.shape)
    out=[]
    for stop in [500,1000,2000]:
        half=stop//2;a=labels[step:half+1:step];b=labels[half+step:stop+1:step]
        assert len(a)==len(b)
        def hist(block):
            return np.stack([np.bincount(block[:,j],minlength=len(unique))/len(block) for j in range(block.shape[1])],axis=1)
        A,B=hist(a),hist(b)
        np.testing.assert_allclose(A.sum(axis=0),1,atol=1e-14)
        np.testing.assert_allclose(B.sum(axis=0),1,atol=1e-14)
        D=A-B;D=D[np.any(D!=0,axis=1)]
        out.append(D)
    return out

def solve(ds,cap):
    n=ds[0].shape[1];sizes=[len(d) for d in ds];m=sum(sizes)
    D=sparse.csr_matrix(np.vstack(ds));I=sparse.eye(m,format='csr')
    # Variables are common orbit weights, per-bin absolute discrepancies, worst TV.
    zero=sparse.csr_matrix((m,1))
    blocks=[sparse.hstack([D,-I,zero]),sparse.hstack([-D,-I,zero])]
    limits=[];start=0
    for size in sizes:
        z=sparse.csr_matrix((np.full(size,.5),(np.zeros(size,int),np.arange(start,start+size))),shape=(1,m))
        limits.append(sparse.hstack([sparse.csr_matrix((1,n)),z,sparse.csr_matrix([[-1.]])]))
        start+=size
    A=sparse.vstack(blocks+limits,format='csr');b=np.zeros(A.shape[0])
    E=sparse.csr_matrix(([1.]*n,([0]*n,range(n))),shape=(1,n+m+1))
    c=np.r_[np.zeros(n+m),1.]
    upper=1. if cap is None else cap
    bounds=[(0,upper)]*n+[(0,None)]*(m+1)
    fit=linprog(c,A_ub=A,b_ub=b,A_eq=E,b_eq=[1.],bounds=bounds,method='highs')
    assert fit.success,fit.message
    w=fit.x[:n];tv=[float(.5*np.abs(d@w).sum()) for d in ds]
    # HiGHS supplies multipliers for <=, equality, lower and upper bounds.
    stationarity=c-A.T@fit.ineqlin.marginals-E.T@fit.eqlin.marginals-fit.lower.marginals-fit.upper.marginals
    dual=float(fit.eqlin.marginals[0]+upper*fit.upper.marginals[:n].sum())
    primal=float(fit.fun);gap=abs(primal-dual)
    residual=max(float(np.max(A@fit.x-b)),abs(float((E@fit.x)[0])-1),float(max(0,-fit.x.min())),float(max(0,w.max()-upper)))
    dual_sign=max(float(max(0,fit.ineqlin.marginals.max())),float(max(0,-fit.lower.marginals.min())),float(max(0,fit.upper.marginals.max())))
    assert residual<1e-7 and np.max(abs(stationarity))<1e-7 and dual_sign<1e-7 and gap<1e-7
    assert abs(max(tv)-primal)<1e-7 and abs(w.sum()-1)<1e-7
    return dict(weights=w.tolist(),worst_TV=primal,window_TV=tv,dual_lower_bound=dual,duality_gap=gap,
        primal_residual=residual,dual_stationarity_residual=float(np.max(abs(stationarity))),
        dual_sign_residual=dual_sign,effective_orbits=float(1/np.sum(w*w)),
        nonzero_weights=int((w>1e-8).sum()),passes_005_target=primal<=.05+1e-7,
        occupied_difference_rows=m,iterations=int(fit.nit))

def main():
    source=json.loads((SOURCE/'results.json').read_text());pairs=source['common_passing_pairs']
    assert pairs==[0,8,12,16,20,24,28,32]
    for path,sha in source['source_hashes'].items():assert digest(ROOT/path)==sha
    hashes={str((SOURCE/'results.json').relative_to(ROOT)):digest(SOURCE/'results.json')}
    rows=[]
    for kind in ['ordinary','full']:
        paths=[]
        for index in pairs:
            path=CACHE/f'{kind}-{index}.npz'
            hashes[str(path.relative_to(ROOT))]=digest(path)
            with np.load(path) as f:paths.append(f['trajectory'])
        for dimensions,name in [(2,'R_z'),(3,'spatial'),(6,'position_velocity')]:
            ds=matrices(paths,dimensions,1);coarse=matrices(paths,dimensions,2)
            for cap in [None,.25]:
                result=solve(ds,cap);thin=solve(coarse,cap)
                long=solve(ds[-1:],cap);long_thin=solve(coarse[-1:],cap)
                w=np.array(result['weights']);thin_w=np.array(thin['weights'])
                result.update(kind=kind,partition=name,weight_cap=cap,
                    equal_weights_TV=[float(.5*np.abs(d@np.full(8,.125)).sum()) for d in ds],
                    half_sampling_refit=thin,
                    half_sampling_frozen_weights_TV=[float(.5*np.abs(d@w).sum()) for d in coarse],
                    full_sampling_thin_weights_TV=[float(.5*np.abs(d@thin_w).sum()) for d in ds])
                result.update(longest_only=long,longest_only_half_sampling=long_thin,
                    longest_only_frozen_weights_half_sampling_TV=float(.5*np.abs(coarse[-1]@np.array(long['weights'])).sum()))
                rows.append(result)
                print(kind,name,cap,'best worst TV',round(result['worst_TV'],6),'longest only',round(long['worst_TV'],6),flush=True)
    # Analytic controls: identical distributions and cancelling two-orbit drifts.
    for ds,expected in [([np.zeros((2,8))]*3,0.),([np.tile([[.2],[-.2]],(1,8))]*3,.2),
                        ([np.tile([[.2,-.2],[-.2,.2]],(1,4))]*3,0.)]:
        for cap in [None,.25]:assert abs(solve(ds,cap)['worst_TV']-expected)<1e-8
    for p,sha in hashes.items():assert digest(ROOT/p)==sha
    data=dict(scope='Best possible finite-bin common-weight occupation stability; not a stellar fit or gravity ranking.',
        pairs=pairs,rows=rows,numerical_TV_target=.05,windows_kpc_per_kms=[.25,.5,1.],
        velocity_edges_kms=[[str(x) for x in e] for e in EDGES[3:]],
        cap_is_diagnostic_not_prior=True,holdouts_opened=False,stationarity_proved=False,
        gravitational_parameters_fitted=False,source_hashes=hashes,
        code_sha256=digest(Path(__file__)),protocol_sha256=digest(HERE/'protocol.md'),analytic_controls=6)
    (HERE/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

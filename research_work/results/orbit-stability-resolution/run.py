"""Exact empirical kernel discrepancies for existing paired orbit populations."""
from pathlib import Path
import importlib.util
import json
import hashlib
import numpy as np
from scipy.spatial.distance import cdist

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PRIOR=HERE.parent/'orbit-mixture-stability'
spec=importlib.util.spec_from_file_location('resolution_prior',PRIOR/'run.py')
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)

def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def coefficients(norbit,ntime):
    assert ntime%2==0
    C=np.zeros((norbit*ntime,2*norbit))
    half=ntime//2
    for j in range(norbit):
        C[j*ntime:(j+1)*ntime,j]=np.r_[np.ones(half),-np.ones(half)]/half
        C[j*ntime:(j+1)*ntime,norbit+j]=np.tile([1.,-1.],half)/half
    return C

def kernel_grams(paths,scales):
    norbit,ntime,_=paths.shape
    points=paths.reshape(-1,6);C=coefficients(norbit,ntime)
    grams=[np.zeros((2*norbit,2*norbit)) for _ in scales]
    for start in range(0,len(points),128):
        stop=min(start+128,len(points));block=points[start:stop]
        dx=cdist(block[:,:3],points[:,:3],metric='sqeuclidean')
        dp=cdist(block[:,3:],points[:,3:],metric='sqeuclidean')
        for index,(sx,sp) in enumerate(scales):
            K=np.exp(-.5*(dx/sx**2+dp/sp**2)).reshape(len(block),norbit,ntime)
            chronological=(K[:,:,:ntime//2].sum(axis=2)-K[:,:,ntime//2:].sum(axis=2))/(ntime//2)
            interleaved=(K[:,:,::2].sum(axis=2)-K[:,:,1::2].sum(axis=2))/(ntime//2)
            grams[index]+=C[start:stop].T@np.c_[chronological,interleaved]
    for gram in grams:
        assert np.max(abs(gram-gram.T))<1e-10
        gram[:]=(gram+gram.T)/2
        assert np.linalg.eigvalsh(gram).min()>-1e-10
    return grams

def tvs(paths,weights):
    vals=np.stack([prior.features(p) for p in paths])
    bins=np.stack([np.searchsorted(prior.EDGES[j],vals[:,:,j],side='right')-1 for j in range(6)],axis=2)
    labels=np.ravel_multi_index(bins.reshape(-1,6).T,tuple(len(e)-1 for e in prior.EDGES)).reshape(vals.shape[:2])
    _,labels=np.unique(labels,return_inverse=True);labels=labels.reshape(vals.shape[:2])
    n=int(labels.max()+1);out={}
    for name,(a,b) in {'chronological':(labels[:,:1000],labels[:,1000:]),'interleaved':(labels[:,::2],labels[:,1::2])}.items():
        D=np.stack([(np.bincount(x,minlength=n)-np.bincount(y,minlength=n))/len(x) for x,y in zip(a,b)],axis=1)
        out[name]={key:float(.5*np.abs(D@w).sum()) for key,w in weights.items()}
    return out

def main():
    base=json.loads((PRIOR/'results.json').read_text())
    assert digest(PRIOR/'run.py')==base['code_sha256']
    for name,sha in base['source_hashes'].items():assert digest(ROOT/name)==sha
    # Exact small finite Gram check, rather than an implementation-mirroring test.
    toy=np.zeros((1,2,6));toy[0,1,0]=2.
    gram=kernel_grams(toy,[(1.,1.)])[0]
    expected=2*(1-np.exp(-2))
    np.testing.assert_allclose(gram,np.full((2,2),expected),atol=1e-13)
    np.testing.assert_allclose(kernel_grams(np.zeros((1,2,6)),[(1.,1.)])[0],0,atol=1e-13)
    toy=np.random.default_rng(174).normal(size=(2,4,6));flat=toy.reshape(-1,6)
    direct=np.exp(-.5*cdist(flat,flat,metric='sqeuclidean'))
    C=coefficients(2,4)
    np.testing.assert_allclose(kernel_grams(toy,[(1.,1.)])[0],C.T@direct@C,atol=1e-13)
    rows=[];hashes={str((PRIOR/'results.json').relative_to(ROOT)):digest(PRIOR/'results.json')}
    for kind in ['ordinary','full']:
        paths=[]
        for i in base['pairs']:
            p=prior.CACHE/f'{kind}-{i}.npz';hashes[str(p.relative_to(ROOT))]=digest(p)
            with np.load(p) as f:paths.append(f['trajectory'][1:])
        paths=np.stack(paths);assert paths.shape==(8,2000,6)
        oldrow=next(r for r in base['rows'] if r['kind']==kind and r['partition']=='position_velocity' and r['weight_cap']==.25)
        weights={'equal':np.full(8,.125),'frozen_capped_longest':np.array(oldrow['longest_only']['weights'])}
        hist=tvs(paths,weights)
        assert abs(hist['chronological']['frozen_capped_longest']-oldrow['longest_only']['worst_TV'])<1e-8
        scales=[(.5,50.),(2.,100.)]
        grams=kernel_grams(paths,scales)
        for (sx,sp),gram in zip(scales,grams):
            entry=dict(kind=kind,position_scale_kpc=sx,velocity_scale_kms=sp,
                       gram_min_eigenvalue=float(np.linalg.eigvalsh(gram).min()),
                       weights={k:v.tolist() for k,v in weights.items()},histogram_TV=hist,comparisons={})
            for name,sl in [('chronological',slice(0,8)),('interleaved',slice(8,16))]:
                G=gram[sl,sl]
                entry['comparisons'][name]=dict(individual_MMD=np.sqrt(np.maximum(0,np.diag(G))).tolist(),
                    mixtures={k:float(np.sqrt(max(0,w@G@w))) for k,w in weights.items()},gram=G.tolist())
            rows.append(entry)
            print(kind,sx,sp,entry['comparisons']['chronological']['mixtures'],entry['comparisons']['interleaved']['mixtures'],flush=True)
    for name,sha in hashes.items():assert digest(ROOT/name)==sha
    result=dict(scope='Resolution audit of empirical orbit occupations, not equilibrium proof or observed likelihood.',
        rows=rows,pairs=base['pairs'],analytic_controls=2,direct_multi_orbit_control=True,source_hashes=hashes,
        code_sha256=digest(Path(__file__)),protocol_sha256=digest(HERE/'protocol.md'),
        previous_histogram_reproduced=True,holdouts_opened=False,new_orbits_integrated=False,
        statistical_significance_claimed=False,new_weights_fitted=False)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.optimize import least_squares
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=ROOT/'redshift_paper/all_164_groups.csv'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
rows=list(csv.DictReader(source.open()))
d=np.array([float(r['catalog_distance_mpc']) for r in rows]);z=np.array([float(r['observed_cmb_z']) for r in rows])
regions=np.array([int((float(r['ra_deg'])%360)//90)+4*int(float(r['dec_deg'])>=0) for r in rows])
tiles=np.array([r['sky_tile'] for r in rows])
assert all(len(set(regions[tiles==tile]))==1 for tile in set(tiles))
C=299792.458;names=['linear','constant','smooth'];out={n:np.zeros(len(d)) for n in names};folds=[]
def pred(name,p,x):
    depth=p[0]*x/C
    if name=='linear':return depth
    if name=='smooth':depth+=(p[1]-p[0])*x*x/(200*C)
    return np.expm1(depth)
def metrics(v):return dict(n=len(v),rms=float(np.sqrt(np.mean(v*v))),mae=float(np.mean(abs(v))),median_absolute=float(np.median(abs(v))),bias=float(v.mean()))
for region in np.unique(regions):
    held=regions==region;train=~held
    for name in names:
        fit=least_squares(lambda p:C*(pred(name,p,d[train])-z[train]),[75.]*(2 if name=='smooth' else 1),bounds=(0,150),xtol=1e-12,ftol=1e-12,gtol=1e-10)
        assert fit.success
        out[name][held]=pred(name,fit.x,d[held])
        folds.append(dict(region=int(region),model=name,parameters=fit.x.tolist(),boundary_hit=bool(np.any(fit.x<1e-5)|np.any(fit.x>150-1e-5)),training_distance_range=[float(d[train].min()),float(d[train].max())],held_distance_range=[float(d[held].min()),float(d[held].max())],held_outside_training_range=int(np.sum((d[held]<d[train].min())|(d[held]>d[train].max()))),**metrics(C*(out[name][held]-z[held]))))
res={name:C*(out[name]-z) for name in names};rng=np.random.default_rng(2026090904)
diff={n:[] for n in ['linear','smooth']}
for _ in range(2000):
    idx=np.concatenate([np.flatnonzero(regions==r) for r in rng.choice(np.unique(regions),len(np.unique(regions)),replace=True)])
    for n in diff:diff[n].append(metrics(res[n][idx])['rms']-metrics(res['constant'][idx])['rms'])
results=dict(protocol_sha256=hashlib.sha256((HERE/'coarse-sky-protocol.md').read_bytes()).hexdigest(),region_counts={str(r):int(sum(regions==r)) for r in np.unique(regions)},out_of_fold={n:metrics(res[n]) for n in names},difference_from_constant_95={n:np.percentile(v,[2.5,97.5]).tolist() for n,v in diff.items()},parameter_ranges={n:dict(minimum=np.min([f['parameters'] for f in folds if f['model']==n],axis=0).tolist(),maximum=np.max([f['parameters'] for f in folds if f['model']==n],axis=0).tolist()) for n in names})
(HERE/'coarse-sky-results.json').write_text(json.dumps(results,indent=2)+'\n',newline='\n')
(HERE/'coarse-sky-folds.json').write_text(json.dumps(folds,indent=2)+'\n',newline='\n')
with (HERE/'coarse-sky-predictions.csv').open('w',newline='') as f:
    w=csv.writer(f,lineterminator='\n');w.writerow(['pgc','group_pgc','region','distance_mpc','observed_z']+[n+'_oof_z' for n in names])
    for i,r in enumerate(rows):w.writerow([r['pgc'],r['group_pgc'],regions[i],d[i],z[i]]+[out[n][i] for n in names])
print(json.dumps(results,indent=2))

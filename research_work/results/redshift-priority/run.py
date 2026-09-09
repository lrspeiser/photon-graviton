from pathlib import Path
import csv, json, hashlib
import numpy as np
from scipy.optimize import least_squares

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
source=ROOT/'redshift_paper/all_164_groups.csv'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='8a2044337ecfe108e56c9592d03d053d48169a1ef0c34405437a34f69a2844a0'
rows=list(csv.DictReader(source.open()))
d=np.array([float(r['catalog_distance_mpc']) for r in rows])
z=np.array([float(r['observed_cmb_z']) for r in rows])
tiles=np.array([int(r['sky_tile']) for r in rows]); unique=np.unique(tiles)
c=299792.458
names=['linear','constant','smooth']
def predict(name,p,x):
    if name=='linear': return p[0]*x/c
    depth=p[0]*x/c
    if name=='smooth':depth+=(p[1]-p[0])*x*x/(200*c)
    return np.expm1(depth)
def fit(name,mask):
    f=least_squares(lambda p:c*(predict(name,p,d[mask])-z[mask]),[75.]*(2 if name=='smooth' else 1),bounds=(0,150),xtol=1e-12,ftol=1e-12,gtol=1e-10)
    assert f.success
    return f.x
def stats(r):
    return dict(n=len(r),rms=float(np.sqrt(np.mean(r*r))),mae=float(np.mean(abs(r))),median_absolute=float(np.median(abs(r))),bias=float(r.mean()))
preds={n:np.zeros(len(d)) for n in names}; folds=[]
for tile in unique:
    held=tiles==tile; training=~held
    assert not set(tiles[held])&set(tiles[training])
    for name in names:
        p=fit(name,training);preds[name][held]=predict(name,p,d[held])
        folds.append(dict(tile=int(tile),model=name,parameters=p.tolist(),**stats(c*(preds[name][held]-z[held]))))
res={n:c*(preds[n]-z) for n in names}
rng=np.random.default_rng(2026090903);diff={n:[] for n in ['linear','smooth']}
for _ in range(2000):
    indices=np.concatenate([np.flatnonzero(tiles==t) for t in rng.choice(unique,len(unique),replace=True)])
    for name in diff:diff[name].append(stats(res[name][indices])['rms']-stats(res['constant'][indices])['rms'])
baseline=c*(np.expm1(.0002488993286382367*d)-z)
ra=np.deg2rad([float(r['ra_deg']) for r in rows]);dec=np.deg2rad([float(r['dec_deg']) for r in rows])
features={'distance':d,'sky_x':np.cos(dec)*np.cos(ra),'sky_y':np.cos(dec)*np.sin(ra),'sky_z':np.sin(dec)}
bins=[]
for indices in np.array_split(np.argsort(d,kind='stable'),4):
    bins.append(dict(distance_min=float(d[indices].min()),distance_max=float(d[indices].max()),**stats(baseline[indices])))
result=dict(protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),groups=len(rows),sky_tiles=len(unique),out_of_fold={n:stats(res[n]) for n in names},difference_from_constant_95={n:np.percentile(v,[2.5,97.5]).tolist() for n,v in diff.items()},baseline_correlations={n:float(np.corrcoef(x,baseline)[0,1]) for n,x in features.items()},baseline_distance_quartiles=bins,environment_inputs_available=False,scope='Previously exposed exploratory sky-tile cross-validation; nearby tiles and calibrations can remain correlated.')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
(HERE/'folds.json').write_text(json.dumps(folds,indent=2)+'\n',newline='\n')
with (HERE/'predictions.csv').open('w',newline='') as f:
    w=csv.writer(f,lineterminator='\n');w.writerow(['pgc','group_pgc','sky_tile','distance_mpc','observed_z']+[n+'_oof_z' for n in names])
    for i,r in enumerate(rows):w.writerow([r['pgc'],r['group_pgc'],tiles[i],d[i],z[i]]+[preds[n][i] for n in names])
print(json.dumps(result,indent=2))

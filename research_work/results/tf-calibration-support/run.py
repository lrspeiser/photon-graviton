from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.spatial import ConvexHull
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
parent=ROOT/'research_work/results/dense-distance-source-audit/observer-features.csv'
calpath=ROOT/'research_work/results/indicator-tf-calibration/cross-validation.csv'
cal={int(r['pgc']):r for r in csv.DictReader(calpath.open()) if r['photometry_p']=='1'}
source=list(csv.DictReader(parent.open()));rows=[];excluded=[]
def number(r,k):
    try:
        value=float(r[k]);return value if value!=-100 else np.nan
    except ValueError:return np.nan
for r in source:
    i,ei,w,ew,m,ai=[number(r,k) for k in ['Inc','e_Inc','Wmx','e_Wmx','icmag','Ai']]
    reason=''
    if not np.all(np.isfinite([i,ei,w,ew,m,ai])):reason='missing_required_measurement'
    elif not (w>0 and ew>0 and ei>0 and 45<=i<=90 and ew/w<=.2):reason='declared_geometry_or_width_quality'
    pgc=int(r['PGC'])
    if reason:excluded.append(dict(pgc=pgc,name=r['Name'],reason=reason));continue
    rows.append(dict(pgc=pgc,name=r['Name'],x=np.log10(w/np.sin(np.deg2rad(i)))-2.5,m=m-ai,inclination=i,calibrator=pgc in cal))
assert len(source)==10737 and len(rows)+len(excluded)==len(source)
assert len({r['pgc'] for r in rows+excluded})==len(source)
anchors=[r for r in rows if r['calibrator']];assert len(anchors)==73
assert max(abs(r[k]-float(cal[r['pgc']][k])) for r in anchors for k in ['x','m'])<1e-12
features=['x','m','inclination'];A=np.array([[r[k] for k in features] for r in anchors]);P=np.array([[r[k] for k in features] for r in rows])
lo=A.min(axis=0);hi=A.max(axis=0)
for j,r in enumerate(rows):
    for k,col in enumerate(features):r['outside_'+col+'_range']=bool(P[j,k]<lo[k]-1e-12 or P[j,k]>hi[k]+1e-12)
for dim in [2,3]:
    scale=np.ptp(A[:,:dim],axis=0);assert np.all(scale>0)
    hull=ConvexHull((A[:,:dim]-lo[:dim])/scale)
    inside=np.all(((P[:,:dim]-lo[:dim])/scale)@hull.equations[:,:-1].T+hull.equations[:,-1]<=1e-10,axis=1)
    for j,r in enumerate(rows):r['inside_hull_'+str(dim)+'d']=bool(inside[j])
    assert all(r['inside_hull_'+str(dim)+'d'] for r in anchors)
def csvwrite(name,data):
    with (HERE/name).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
csvwrite('support.csv',rows);csvwrite('excluded.csv',excluded)
bins=[]
for lower in range(0,24,2):
    subset=[r for r in rows if lower<=r['m']<lower+2]
    bins.append(dict(lower_mag=lower,upper_mag=lower+2,parent_n=len(subset),calibrator_n=sum(r['calibrator'] for r in subset)))
assert sum(r['parent_n'] for r in bins)==len(rows)
result=dict(source_n=len(source),eligible_n=len(rows),excluded_n=len(excluded),calibrators=73,
    feature_summaries={k:dict(calibrator_min=float(lo[j]),calibrator_max=float(hi[j]),calibrator_percentiles_5_50_95=np.percentile(A[:,j],[5,50,95]).tolist(),parent_percentiles_5_50_95=np.percentile(P[:,j],[5,50,95]).tolist(),parent_outside_calibrator_range=sum(r['outside_'+k+'_range'] for r in rows)) for j,k in enumerate(features)},
    parent_inside_hull={str(dim):sum(r['inside_hull_'+str(dim)+'d'] for r in rows) for dim in [2,3]},magnitude_bins=bins,
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [parent,calpath]},
    status='Empirical feature coverage, not physical completeness, selection probability or independent validation')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8',newline='\n');print(json.dumps(result,indent=2))

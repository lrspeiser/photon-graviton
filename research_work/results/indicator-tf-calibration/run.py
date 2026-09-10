from pathlib import Path
import csv, json, hashlib
import numpy as np
from scipy.optimize import brentq
from scipy.special import lambertw

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ALPHA=0.0002488993286382367
base=ROOT/'research_work/results/dense-distance-source-audit'
paths=[base/'possible-calibrators.csv',base/'observer-features.csv',ROOT/'temporal_candidate_audit/data/cf4_table2.dat']
cal=list(csv.DictReader(paths[0].open()))
features={int(r['PGC']):r for r in csv.DictReader(paths[1].open())}
groups={int(s[:7]):int(s[8:15]) for s in paths[2].read_text().splitlines() if s[:7].strip().isdigit()}
def number(r,k):
    try:
        v=float(r[k]); return v if v!=-100 else np.nan
    except (ValueError,KeyError): return np.nan
rows=[]; exclusions=[]
for c in cal:
    pgc=int(c['pgc']); f=features[pgc]
    inc,ei,w,ew,mi,ai=[number(f,k) for k in ['Inc','e_Inc','Wmx','e_Wmx','icmag','Ai']]
    reason=''
    if not np.all(np.isfinite([inc,ei,w,ew,mi,ai])): reason='missing_required_measurement'
    elif not (w>0 and ew>0 and ei>0 and 45<=inc<=90 and ew/w<=0.2): reason='declared_geometry_or_width_quality'
    if reason:
        exclusions.append(dict(pgc=pgc,name=c['name'],reason=reason)); continue
    angle=np.deg2rad(inc)
    x=np.log10(w/np.sin(angle))-2.5
    sx=np.sqrt((ew/w)**2+(np.deg2rad(ei)/np.tan(angle))**2)/np.log(10)
    rows.append(dict(pgc=pgc,name=c['name'],method=c['method'],group=groups.get(pgc,pgc) or pgc,
        distance=float(c['distance_mpc']),sigma_mu=float(c['sigma_mag']),m=mi-ai,x=x,sigma_x=sx))
def write_csv(name,data):
    with (HERE/name).open('w',encoding='utf-8',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(data[0]),lineterminator='\n');writer.writeheader();writer.writerows(data)
write_csv('excluded.csv',exclusions)
X=np.array([[r['x'],1] for r in rows]); d=np.array([r['distance'] for r in rows]); m=np.array([r['m'] for r in rows]); g=np.array([r['group'] for r in rows])
def fit(x,y):
    coef,_,rank,_=np.linalg.lstsq(x,y,rcond=None); assert rank==2; return coef
results=[]; output=[]; root_errors=[]
for p in [0,1]:
    y=m-5*np.log10(d)-25-5*p*ALPHA*d/np.log(10)
    coef=fit(X,y); pred=np.full(len(rows),np.nan); seen=np.zeros(len(rows),int)
    for group in np.unique(g):
        test=g==group; pred[test]=X[test]@fit(X[~test],y[~test]);seen[test]+=1
    assert np.all(seen==1) and np.all(np.isfinite(pred))
    residual=y-pred; dinferred=[]
    for j,r in enumerate(rows):
        apparent=10**((m[j]-pred[j]-25)/5)
        inferred=apparent if p==0 else float(lambertw(p*ALPHA*apparent).real/(p*ALPHA))
        independent=brentq(lambda value:value*np.exp(p*ALPHA*value)-apparent,0,apparent)
        root_errors.append(abs(inferred-independent)/inferred); dinferred.append(inferred)
        output.append(dict(**r,photometry_p=p,observed_absolute_magnitude=y[j],group_excluded_prediction=pred[j],residual_mag=residual[j],group_excluded_distance=inferred,distance_ratio=inferred/d[j]))
    ratio=np.array(dinferred)/d
    results.append(dict(p=p,slope=float(coef[0]),intercept=float(coef[1]),training_rms_mag=float(np.sqrt(np.mean((y-X@coef)**2))),
        group_excluded_rms_mag=float(np.sqrt(np.mean(residual**2))),group_excluded_median_abs_mag=float(np.median(abs(residual))),
        distance_ratio_percentiles=dict(zip(['p16','p50','p84'],np.percentile(ratio,[16,50,84]).tolist())),
        median_absolute_fractional_distance_error=float(np.median(abs(ratio-1))),
        median_width_induced_mag_error=float(np.median([abs(coef[0])*r['sigma_x'] for r in rows])),
        by_method={method:dict(n=sum(r['method']==method for r in rows),rms_mag=float(np.sqrt(np.mean(residual[np.array([r['method']==method for r in rows])]**2)))) for method in sorted({r['method'] for r in rows})}))
assert max(root_errors)<1e-10
write_csv('cross-validation.csv',output)
result=dict(alpha_fixed=ALPHA,eligible=len(rows),excluded=len(exclusions),groups=len(np.unique(g)),distance_range_mpc=[float(min(d)),float(max(d))],
    source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},scenarios=results,max_relative_root_difference=max(root_errors),
    scope='Exploratory group-excluded calibration; no fresh holdout, no environment or redshift fit; not a complete likelihood.')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))

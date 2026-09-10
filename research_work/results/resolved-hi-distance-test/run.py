from pathlib import Path
import csv,json,re,hashlib
import pandas as pd
import numpy as np
from scipy.special import lambertw
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
html=ROOT/'research_work/data-cache/resolved-hi/paper.html'
tables=pd.read_html(html)
def identity(s):
    m=re.fullmatch(r'(NGC|IC)\s*0*(\d+)',s.strip(),re.I)
    return m[1].upper()+str(int(m[2])) if m else s.strip().upper()
def numbers(v): return [float(s) for s in re.findall(r'\d+(?:\.\d+)?',str(v))]
profiles={identity(str(r.iloc[0])):str(r.iloc[6]) for _,r in tables[5].iloc[1:].iterrows()}
atlas=[]
for _,r in tables[6].iloc[1:].iterrows():
    vflat=numbers(r.iloc[9]); assert len(vflat) in [0,2]
    atlas.append(dict(name=str(r.iloc[0]),identity=identity(str(r.iloc[0])),kinematic_inclination=numbers(r.iloc[5])[0],
        vmax=numbers(r.iloc[7])[0],e_vmax=numbers(r.iloc[8])[0],vflat=vflat[0] if vflat else None,
        e_vflat=vflat[1] if vflat else None,profile=profiles[identity(str(r.iloc[0]))]))
assert len(atlas)==32 and len({r['identity'] for r in atlas})==32
calpath=ROOT/'research_work/results/indicator-tf-calibration/cross-validation.csv'
cal=[r for r in csv.DictReader(calpath.open()) if r['photometry_p']=='1']
lookup={identity(r['name']):r for r in cal}; assert len(lookup)==len(cal)
matched=[]; coverage=[]
for a in atlas:
    c=lookup.get(a['identity'])
    coverage.append(dict(**a,calibrator_match=bool(c),common_comparison=bool(c and a['vflat'] is not None)))
    if c and a['vflat'] is not None:
        matched.append(dict(**a,pgc=c['pgc'],group=int(c['group']),distance=float(c['distance']),m=float(c['m']),
            x_width=float(c['x']),y=float(c['observed_absolute_magnitude']),old_residual=float(c['residual_mag'])))
def csvwrite(name,rows):
    with (HERE/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
csvwrite('atlas-coverage.csv',coverage)
alpha=0.0002488993286382367; group=np.array([r['group'] for r in matched]); y=np.array([r['y'] for r in matched])
def fit(X,y):
    a,_,rank,_=np.linalg.lstsq(X,y,rcond=None);assert rank==2;return a
results=[]; predictions=[]; errors=[]
for label in ['width','vmax','vflat']:
    x=np.array([r['x_width'] if label=='width' else np.log10(2*r[label])-2.5 for r in matched]);X=np.column_stack([x,np.ones(len(x))]);coef=fit(X,y)
    predicted=np.full(len(x),np.nan);seen=np.zeros(len(x),int)
    for g in np.unique(group):
        held=group==g;predicted[held]=X[held]@fit(X[~held],y[~held]);seen[held]+=1
    assert np.all(seen==1) and np.all(np.isfinite(predicted))
    ratios=[]
    for j,r in enumerate(matched):
        app=10**((r['m']-predicted[j]-25)/5);d=float(lambertw(alpha*app).real/alpha)
        check=brentq(lambda v:v*np.exp(alpha*v)-app,0,app);errors.append(abs(check-d)/d);ratios.append(d/r['distance'])
        predictions.append(dict(**r,predictor=label,x=x[j],group_excluded_prediction=predicted[j],residual=y[j]-predicted[j],predicted_distance=d,distance_ratio=ratios[-1]))
    results.append(dict(predictor=label,slope=float(coef[0]),intercept=float(coef[1]),rms_mag=float(np.sqrt(np.mean((y-predicted)**2))),
        median_abs_fractional_distance_error=float(np.median(abs(np.array(ratios)-1))),ratio_percentiles=np.percentile(ratios,[16,50,84]).tolist()))
assert max(errors)<1e-10
csvwrite('predictions.csv',predictions)
result=dict(atlas_rows=len(atlas),matches=sum(r['calibrator_match'] for r in coverage),common_sample=len(matched),groups=len(np.unique(group)),
    distance_range=[min(r['distance'] for r in matched),max(r['distance'] for r in matched)],
    original_73_fit_residual_rms_on_overlap=float(np.sqrt(np.mean([r['old_residual']**2 for r in matched]))),
    scenarios=results,max_relative_root_difference=max(errors),
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [html,calpath]},
    status='Exposed exploratory paired comparison, not independent holdout or complete errors/selection model')
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(result,indent=2))

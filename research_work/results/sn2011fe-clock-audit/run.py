"""Observed single-instrument tail diagnostic, not a propagation fit."""
from pathlib import Path
import urllib.request, hashlib, json, re
import numpy as np
from bs4 import BeautifulSoup
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/generated/sn2011fe-clock-audit'
CACHE.mkdir(parents=True,exist_ok=True)
URL='https://arxiv.org/html/1701.07267v2'
rawfile=CACHE/'1701.07267v2.html'
if not rawfile.exists():
    rawfile.write_bytes(urllib.request.urlopen(URL,timeout=60).read())
raw=rawfile.read_bytes();s=BeautifulSoup(raw,'html.parser')
rows=[]
for tr in s.find(id='S2.T1.2').select('tr'):
    td=tr.find_all('td',recursive=False)
    if len(td)!=7: continue
    vals=[]
    for cell in td:
        for m in cell.find_all('math'): m.replace_with(m['alttext'])
        vals.append(cell.get_text(' ',strip=True))
    if vals[2]!='LBT' or vals[3]!='R': continue
    match=re.fullmatch(r'([0-9.]+)\(([0-9.]+)\)',vals[4])
    if match is None: raise ValueError(vals)
    rows.append(dict(mjd=float(vals[0]),published_phase_days=float(vals[1]),telescope=vals[2],filter=vals[3],magnitude=float(match[1]),sigma_magnitude=float(match[2]),system=vals[5],source=vals[6]))
assert len(rows)>10 and all(r['system']=='Vega' for r in rows)
assert len({r['mjd'] for r in rows})==len(rows)
# Use actual observation dates. Reference epoch only labels windows; it
# does not rescale time and a change of origin cannot alter fitted slope.
t=np.array([r['mjd']-55814.30 for r in rows]);y=np.array([r['magnitude'] for r in rows]);err=np.array([r['sigma_magnitude'] for r in rows])
results=[]
# Exploratory windows, not a preregistered or untouched hypothesis test.
for label,mask in [('all',np.ones(len(t),dtype=bool)),('before_900_days',t<900),('after_900_days',t>=900)]:
    x=t[mask];m=y[mask];e=err[mask];center=np.average(x,weights=1/e**2)
    X=np.column_stack([np.ones(len(x)),x-center]);coef=np.linalg.lstsq(X/e[:,None],m/e,rcond=None)[0]
    covariance=np.linalg.inv((X/e[:,None]).T@(X/e[:,None]))
    pred=X@coef;resid=m-pred;chi=float(np.sum((resid/e)**2));slope=float(coef[1])
    results.append(dict(window=label,n=len(x),mjd_min=float(x.min()+55814.30),mjd_max=float(x.max()+55814.30),slope_mag_per_100_observer_days=100*slope,formal_slope_sigma=100*float(np.sqrt(covariance[1,1])),chi2_diagonal_only=chi,dof=len(x)-2,rms_magnitude=float(np.sqrt(np.mean(resid**2))),apparent_decay_stretch_if_band_traces_single_isotope=float(2.5/np.log(10)/111.3/slope),predicted_magnitudes=pred.tolist(),residual_magnitudes=resid.tolist()))
manifest=dict(url=URL,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),table='S2.T1.2',selection='Telescope LBT, filter R, retain every matching row',timing='MJD as published, no rest-frame phase rescaling',role='Exposed exploratory source diagnostic; not held-out cosmological validation')
for name,value in [('source-manifest.json',manifest),('observations.json',rows),('results.json',dict(scope='Single-band constant exponential diagnostic only; no bolometric or propagation inference',tau_days_benchmark=111.3,constant_decay_slope_mag_per_100_days=100*2.5/np.log(10)/111.3,fits=results))]:
    (HERE/name).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([{k:v for k,v in r.items() if k not in ['predicted_magnitudes','residual_magnitudes']} for r in results]))

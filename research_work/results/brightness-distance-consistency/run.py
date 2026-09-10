from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.special import lambertw
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
source=ROOT/'redshift_paper/all_164_groups.csv'
rate_source=ROOT/'research_work/results/electromagnetic-audit/protocol.json'
alpha=json.loads(rate_source.read_text())['alpha_per_mpc'];c=299792.458
rows=list(csv.DictReader(source.open()))
d=np.array([float(r['catalog_distance_mpc']) for r in rows])
mu=np.array([float(r['sbf_modulus_mag']) for r in rows])
sigma=np.array([float(r['sbf_modulus_error_mag']) for r in rows])
z=np.array([float(r['observed_cmb_z']) for r in rows])
splits=np.array([r['split'] for r in rows]);assert len(rows)==164
assert np.allclose(d,10**((mu-25)/5),rtol=1e-12)
def invert(app,p,anchor):
    if p==0:return np.array(app,copy=True)
    return lambertw(p*alpha*app*np.exp(p*alpha*anchor)).real/(p*alpha)
def stats(r):return dict(n=len(r),rms_kms=float(np.sqrt(np.mean(r*r))),bias_kms=float(np.mean(r)),mae_kms=float(np.mean(abs(r))))
baseline=np.expm1(alpha*d);output=[];scores=[];maxroot=0.;maxderiv=0.
for p,anchor in [(0,0)]+[(p,a) for p in [.5,1] for a in [0,1,10]]:
    x=invert(d,p,anchor);pred=np.expm1(alpha*x);res=c*(pred-z)
    derivative=np.log(10)/5*x/(1+p*alpha*x)
    step=1e-5
    num=(invert(10**((mu+step-25)/5),p,anchor)-invert(10**((mu-step-25)/5),p,anchor))/(2*step)
    maxderiv=max(maxderiv,float(np.max(abs(num/derivative-1))))
    assert np.max(abs(num/derivative-1))<1e-8
    assert np.max(abs(x*np.exp(p*alpha*(x-anchor))/d-1))<1e-12
    for i in range(0,len(d),20):
        root=brentq(lambda q:q*np.exp(p*alpha*(q-anchor))-d[i],0,2*max(d[i],anchor+1),xtol=1e-12)
        maxroot=max(maxroot,abs(root-x[i]));assert abs(root-x[i])<1e-9
    for i,row in enumerate(rows):
        output.append(dict(pgc=row['pgc'],historical_partition=row['split'],p=p,anchor_mpc=anchor,
            stipulated_distance_mpc=d[i],sensitivity_distance_mpc=x[i],conditional_sigma_distance_mpc=derivative[i]*sigma[i],
            observed_z=z[i],predicted_z=pred[i],residual_kms=res[i],prediction_change_kms=c*(pred[i]-baseline[i])))
    scores.append(dict(p=p,anchor_mpc=anchor,all=stats(res),
        historical_partitions={s:stats(res[splits==s]) for s in sorted(set(splits))},
        distance_change_percent_min_median_max=np.percentile(100*(x/d-1),[0,50,100]).tolist(),
        prediction_change_kms_min_median_max=np.percentile(c*(pred-baseline),[0,50,100]).tolist()))
with (HERE/'predictions.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(output[0]),lineterminator='\n');w.writeheader();w.writerows(output)
result=dict(status='Exposed-data bolometric distance sensitivity only; original distances and rate not changed.',
    alpha_per_mpc=alpha,scenarios=scores,checks=dict(max_bracket_vs_lambert_mpc=maxroot,
        max_derivative_relative_error=maxderiv,forward_distance_roundtrip='pass'),
    sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,rate_source,HERE/'protocol.md']})
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))

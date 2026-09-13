"""Necessary density-slope/anisotropy check for a separable augmented-density completion."""
from pathlib import Path
import hashlib,json
import numpy as np
from numpy.polynomial.legendre import leggauss
P=Path(__file__).resolve().parent
orbits=P/'radial-orbits-results.json';light=P.parent/'slacs-light-profile-audit/results.json'
fit=json.loads(orbits.read_text());profiles={r['Name']:r for r in json.loads(light.read_text())['rows']}
def deproject_slope(profile,x,order):
    nodes,weights=leggauss(order);total=np.zeros_like(x);slope_weight=np.zeros_like(x)
    for comp in profile['components']:
        Re=comp['R_arcsec']/profile['computed_equal_area_half_light_arcsec'];n=comp['n'];bn=comp['bn'];amp=comp['amp_at_R']
        upper=np.arccosh(np.maximum(1.,Re*(1+100/bn)**n/x))
        u=upper[:,None]*(nodes+1)/2;R=x[:,None]*np.cosh(u)
        q=(R/Re)**(1/n)
        integrand=amp*np.exp(-bn*(q-1))*bn/(n*Re)*(R/Re)**(1/n-1)
        density=np.sum(integrand*weights,axis=1)*upper/(2*np.pi)
        numerator=np.sum(integrand*(1-1/n+bn*q/n)*weights,axis=1)*upper/(2*np.pi)
        total+=density;slope_weight+=numerator
    assert np.all(total>0)
    return total,slope_weight/total

x=np.geomspace(1e-6,100,4097);cache={}
out=dict(scope='Necessary gamma>=2 beta for separable augmented density with beta0<=1/2; neither sufficient positivity nor a universal rejection test for arbitrary orbit distributions',
         radius_range_Re=[float(x[0]),float(x[-1])],rows=[],summary=[],profile_checks=[])
for row in fit['rows']:
    name=row['Name'];profile=profiles[name]
    if name not in cache:
        nu,gamma=deproject_slope(profile,x,256)
        refined,gamma2=deproject_slope(profile,x,512)
        difference=float(max(abs(gamma-gamma2)))
        assert difference<1e-6
        finite_difference=-np.gradient(np.log(refined),np.log(x),edge_order=2)
        derivative_error=float(max(abs(finite_difference[2:-2]-gamma2[2:-2])))
        assert derivative_error<.002
        cache[name]=gamma2
        out['profile_checks'].append(dict(Name=name,max_quadrature_slope_difference=difference,max_finite_difference_slope_difference=derivative_error))
    cases=[('original_bounds',row['beta0'],row['beta_infinity'])]
    if 'expanded_outer_bound_followup' in row:
        b=row['expanded_outer_bound_followup'];cases.append(('expanded_outer_bound',b['beta0'],b['beta_infinity']))
    for label,b0,bi in cases:
        assert b0<=.5
        beta=b0+(bi-b0)*x*x/(1+x*x);margin=cache[name]-2*beta;i=int(np.argmin(margin))
        out['rows'].append(dict(Name=name,population=row['population'],case=label,beta0=b0,beta_infinity=bi,
            minimum_slope_minus_twice_beta=float(margin[i]),minimum_radius_Re=float(x[i]),
            sampled_violation=bool(min(margin)<-1e-6),
            margin_at_Re=float(np.interp(0,np.log(x),margin))))
for label in ['original_bounds','expanded_outer_bound']:
    rr=[r for r in out['rows'] if r['case']==label]
    out['summary'].append(dict(case=label,cases=len(rr),sampled_violations=sum(r['sampled_violation'] for r in rr),
        minimum_margin=min(r['minimum_slope_minus_twice_beta'] for r in rr)))
out['input_sha256']={str(f.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(f.read_bytes()).hexdigest() for f in [orbits,light]}
out['theorem_reference']='https://arxiv.org/abs/1010.4301'
(P/'orbit-slope-check-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out['summary'],indent=2));print(out['profile_checks'])

"""Known fractional-response mixture applied to the exact-third reference."""
from pathlib import Path
import math,json
import numpy as np
from scipy.integrate import quad
from scipy.special import expit,gamma
P=Path(__file__).resolve().parent
q=1/3
def density(y):return math.sin(math.pi*q)/(2*math.pi*(math.cosh(q*y)+math.cos(math.pi*q)))
def response(X,limit=150):
 norm=quad(density,-limit,limit,epsabs=1e-12)[0]
 return quad(lambda y:expit(math.log(X)-y)*density(y),-limit,limit,epsabs=1e-12)[0]/norm
out=dict(scope='Inverse-designed threshold distribution, not a derivation of exponent or observed site population',exponent=q,exact_checks=[],truncations=[],release=[])
for X in np.logspace(-6,6,25):
 exact=X**q/(1+X**q);mixed=response(X)
 assert abs(mixed-exact)<1e-10
 out['exact_checks'].append(dict(X=float(X),reference=exact,mixture=mixed))
for decades in [3,6,9]:
 limit=decades*math.log(10);norm=quad(density,-limit,limit,epsabs=1e-12)[0]
 rows=[dict(X=r['X'],reference=r['reference'],mixture=response(r['X'],limit)) for r in out['exact_checks']]
 out['truncations'].append(dict(half_width_decades=decades,retained_distribution_mass=norm,max_absolute_error=max(abs(r['mixture']-r['reference']) for r in rows),rows=rows))
for u in [1.,1e3,1e6,1e9]:
 # X=1 equilibrium before source removal; each local site releases at rate lambda*t.
 val=quad(lambda y:expit(-y)*math.exp(-u*math.exp(y))*density(y),-150,math.log(750/u),epsabs=1e-12)[0]
 asym=gamma(q)*math.sin(math.pi*q)/math.pi*u**(-q)
 out['release'].append(dict(time_in_inverse_lambda=u,remaining_occupied_fraction=val,asymptotic_fraction=float(asym),fraction_of_initial_retained=2*val))
assert abs(out['release'][-1]['remaining_occupied_fraction']/out['release'][-1]['asymptotic_fraction']-1)<.001
(P/'threshold-mixture-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print([(r['half_width_decades'],r['retained_distribution_mass'],r['max_absolute_error']) for r in out['truncations']]);print(out['release'])

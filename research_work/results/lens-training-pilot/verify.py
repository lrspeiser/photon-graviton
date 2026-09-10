"""Independent Jeans kernel and analytic Hernquist virial checks."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
H=Path(__file__).resolve().parent
load=lambda name:json.loads((H/name).read_text())
a=load('predictions.json');b=load('predictions-refined.json')
key=lambda r:(r['Name'],r['model'],r['seeing_fwhm_arcsec'],r['cutoff_over_Re'])
assert [key(r) for r in a]==[key(r) for r in b]
assert len(set(r['Name'] for r in a))==33 and all(r['role']=='training' for r in a)
change=max(abs(r['theta_pred_arcsec']-s['theta_pred_arcsec']) for r,s in zip(a,b))
assert change<.001
# Dimensionless j=1/[x(1+x)^3], Hernquist GM/a²=1.
# Jeans pressure solved first and then aperture-integrated, rather than swapping
# integrals with the analytic K used in run.py.
x=np.geomspace(1e-7,1e6,100000);j=1/(x*(1+x)**3);g=1/(1+x)**2
pressure=-cumulative_trapezoid((j*g)[::-1],x[::-1],initial=0)[::-1]
errors=[]
for ap in [.1,.5,1,3,10]:
    u=np.minimum(1,(ap/x)**2);W=u/(1+np.sqrt(1-u))
    K=x**3/3*np.where(u==1,1,-np.expm1(1.5*np.log1p(-np.minimum(u,1-1e-16))))
    num1=np.trapezoid(pressure*x*x*W,x);num2=np.trapezoid(j*g*K,x)
    errors.append(abs(num1/num2-1))
assert max(errors)<1e-5
whole=quad(lambda x:x*x/(1+x)**5/3,0,np.inf)[0]/quad(lambda x:x/(1+x)**3,0,np.inf)[0]
assert abs(whole-1/18)<1e-12
out={'training_rows':len(a),'training_systems':33,'maximum_refinement_Einstein_angle_change_arcsec':change,
    'maximum_independent_Jeans_kernel_relative_difference':max(errors),'whole_aperture_Hernquist_sigma2_over_GM_by_a':whole,
    'analytic_expected':1/18,'scope':'Numerical operator tests, not a validated physical or observational likelihood.'}
(H/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))

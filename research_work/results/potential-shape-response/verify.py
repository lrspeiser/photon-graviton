"""Check deformation derivatives, numerical source signs and stored metrics."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
s=json.loads((HERE/'results.json').read_text());rows=json.loads((HERE/'predictions.json').read_text())
for name,digest in s['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest
gradient_errors=[];laplacian_errors=[];loop_errors=[]
for q in [.5,1.,1.8464986,4.]:
    def phi(R,z):return -1/np.sqrt(1+R*R+z*z/q**2)
    def force(R,z):
        u=1+R*R+z*z/q**2
        return np.array([-R,-z/q**2])/u**1.5
    for R,z in [(1.,.3),(3.,2.),(10.,5.)]:
        h=1e-4
        num=-np.array([(phi(R+h,z)-phi(R-h,z))/(2*h),(phi(R,z+h)-phi(R,z-h))/(2*h)])
        gradient_errors.append(float(np.max(abs(num-force(R,z)))))
        numerical=-((force(R+h,z)[0]-force(R-h,z)[0])/(2*h)+force(R,z)[0]/R+(force(R,z+h)[1]-force(R,z-h)[1])/(2*h))
        u=1+R*R+z*z/q**2
        exact=((2+1/q**2)*u-3*(R*R+z*z/q**4))/u**2.5
        laplacian_errors.append(abs(numerical-exact))
    work=quad(lambda R:force(R,.3)[0],1,3)[0]+quad(lambda z:force(3,z)[1],.3,2)[0]+quad(lambda R:force(R,2)[0],3,1)[0]+quad(lambda z:force(1,z)[1],2,.3)[0]
    loop_errors.append(abs(work))
assert max(gradient_errors)<1e-7 and max(laplacian_errors)<1e-7 and max(loop_errors)<1e-10
metric_errors=[]
for case in s['cases']:
    for observable,key in [('rotation_training','rotation_rms_kms'),('vertical_exposed','vertical_rms_surface_equivalent')]:
        selected=[r for r in rows if r['case']==case['case'] and r['observable']==observable]
        residual=np.array([r['predicted']-r['observed'] for r in selected])
        value=np.sqrt(np.mean(residual**2));metric_errors.append(abs(value-case['metrics'][key]))
        assert abs(value-case['metrics'][key])<1e-10
    if case['case']=='shape_only':
        assert abs(case['metrics']['rotation_rms_kms']-s['cases'][0]['metrics']['rotation_rms_kms'])<1e-12
assert all(r['negative_points']==0 for r in s['source_diagnostic'])
out=dict(max_analytic_gradient_absolute_error=max(gradient_errors),max_analytic_poisson_absolute_error=max(laplacian_errors),
    max_closed_loop_work=max(loop_errors),max_stored_metric_reconstruction_error=max(metric_errors),
    original_input_hashes_unchanged=True,finite_source_grid_negative_points=0,
    qualification='Toy analytic derivatives verify deformation calculus; the actual field source test is finite-grid only and does not prove global positivity or a companion-production mechanism.')
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(out,indent=2))

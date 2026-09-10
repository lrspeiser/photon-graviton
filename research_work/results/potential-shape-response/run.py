"""Exploratory two-parameter potential deformation; no held-out evaluation."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
from scipy.optimize import minimize_scalar, differential_evolution

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=HERE.parent/'conservative-field-completion'
spec=importlib.util.spec_from_file_location('frozen_field',SOURCE/'run.py')
field=importlib.util.module_from_spec(spec);spec.loader.exec_module(field)
stars=json.loads((HERE.parent/'cepheid-common-frame/training-bins.json').read_text())
old=json.loads((SOURCE/'predictions-refined.json').read_text())
vertical=[r for r in old if r['observable']=='vertical' and r['model']=='ordinary']
b2=np.array([r['ordinary_vc_kms']**2 for r in stars])
c2=np.array([r['completion_vc_kms']**2 for r in stars])-b2
seen=np.array([r['jeans_proxy_kms'] for r in stars])
Rv=np.array([r['R_kpc'] for r in vertical]);zv=np.full(len(Rv),1.1)
kv_seen=np.array([r['observed'] for r in vertical])
kv_b=np.array([r['predicted'] for r in vertical])
conversion=2*np.pi*field.G*1e6
P=field.P
def radial(lam):return np.sqrt(lam*b2+lam**P*c2)
fit=minimize_scalar(lambda lam:np.mean((radial(lam)-seen)**2),bounds=(.5,1.6),method='bounded',options={'xatol':1e-12})
assert fit.success
lam_rot=float(fit.x)
print('Building refined frozen field for vertical deformation',flush=True)
baryons=field.Baryons(True);extra=field.Completion(baryons,True)

def vertical_prediction(lam,q):
    mu=0. if q is None else 1/q
    additional=-field.force(extra,Rv,zv*mu)[:,1]*mu/conversion
    return lam*kv_b+lam**P*additional

def metrics(lam,q):
    vr=radial(lam);kv=vertical_prediction(lam,q)
    return dict(rotation_rms_kms=float(np.sqrt(np.mean((vr-seen)**2))),
        rotation_bias_kms=float(np.mean(vr-seen)),
        vertical_rms_surface_equivalent=float(np.sqrt(np.mean((kv-kv_seen)**2))),
        vertical_bias_surface_equivalent=float(np.mean(kv-kv_seen)),
        equal_observable_log_mse=float(.5*(np.mean(np.log(vr/seen)**2)+np.mean(np.log(kv/kv_seen)**2))))

def qfit(lam):
    result=minimize_scalar(lambda q:np.mean((vertical_prediction(lam,q)-kv_seen)**2),bounds=(.5,4.),method='bounded',options={'xatol':1e-10})
    assert result.success
    # Include edges: bounded minimization alone need not return an endpoint.
    candidates=[(.5,float(np.mean((vertical_prediction(lam,.5)-kv_seen)**2))),
                (4.,float(np.mean((vertical_prediction(lam,4.)-kv_seen)**2))),
                (float(result.x),float(result.fun))]
    return min(candidates,key=lambda x:x[1])[0]

joint=differential_evolution(lambda x:metrics(*x)['equal_observable_log_mse'],[(.5,1.6),(.5,4)],seed=9102026,tol=1e-8,popsize=10,polish=True)
assert joint.success
cases=[]
for name,lam,q in [('original',1.,1.),('shape_only',1.,qfit(1.)),('rotation_mass_original_shape',lam_rot,1.),
                   ('rotation_mass_then_vertical_shape',lam_rot,qfit(lam_rot)),('joint_compromise',float(joint.x[0]),float(joint.x[1]))]:
    cases.append(dict(case=name,lambda_mass=lam,q_potential=q,metrics=metrics(lam,q),
        q_at_search_boundary=bool(q<=.50001 or q>=3.99999)))

# Audit the encountered q=4 boundary via mu=1/q, including the cylindrical limit.
def q_from_mu(mu):return None if mu==0 else 1/float(mu)
extended=differential_evolution(lambda x:metrics(x[0],q_from_mu(x[1]))['equal_observable_log_mse'],[(.5,1.6),(0.,2.)],seed=9102027,tol=1e-9,popsize=10,polish=True)
assert extended.success
limit_fit=minimize_scalar(lambda lam:metrics(lam,None)['equal_observable_log_mse'],bounds=(.5,1.6),method='bounded')
if limit_fit.fun <= extended.fun:
    extended_lam=float(limit_fit.x);extended_q=None
else:
    extended_lam=float(extended.x[0]);extended_q=q_from_mu(extended.x[1])
cases.append(dict(case='extended_joint_with_cylindrical_limit',lambda_mass=extended_lam,q_potential=extended_q,
    metrics=metrics(extended_lam,extended_q),q_at_search_boundary=extended_q is None,
    note='q=null means infinite vertical stretch, a limiting diagnostic, not a finite galaxy source.'))
mu_grid=np.linspace(0,2,161);lambda_grid=np.linspace(.5,1.6,111)
grid_best=np.inf
for mu in mu_grid:
    added=(vertical_prediction(1.,q_from_mu(mu))-kv_b)
    vp=lambda_grid[:,None]*kv_b+lambda_grid[:,None]**P*added
    rp=np.sqrt(lambda_grid[:,None]*b2+lambda_grid[:,None]**P*c2)
    grid_scores=.5*(np.mean(np.log(rp/seen)**2,axis=1)+np.mean(np.log(vp/kv_seen)**2,axis=1))
    grid_best=min(grid_best,float(min(grid_scores)))
assert cases[-1]['metrics']['equal_observable_log_mse']<=grid_best+1e-8

# For an everywhere attractive added vertical force, ordinary matter gives a floor.
floor_residual=np.maximum(lam_rot*kv_b-kv_seen,0)
floor=dict(lambda_mass=lam_rot,rows_ordinary_already_above_observed=int(np.sum(floor_residual>0)),
    unavoidable_vertical_rms_lower_bound=float(np.sqrt(np.mean(floor_residual**2))),
    note='Conditional on fixed ordinary-matter shape, normalization, and published inferred vertical values. Allows arbitrary nonnegative extra pull independently at every row; hence a lower bound, not a realizable fitted model.')

# Field-source diagnostic on a finite cylindrical grid, not a proof of global positivity.
RR,ZZ=np.meshgrid(np.linspace(.5,30,61),np.linspace(0,10,41),indexing='ij')
R=RR.ravel();z=ZZ.ravel()
def deformed_force(R,z,q):
    a=field.force(extra,R,z/q);a[:,1]/=q
    return a
def source(q,h):
    a=deformed_force(R,z,q)
    derivative_R=(deformed_force(R+h,z,q)[:,0]-deformed_force(R-h,z,q)[:,0])/(2*h)
    derivative_z=(deformed_force(R,z+h,q)[:,1]-deformed_force(R,z-h,q)[:,1])/(2*h)
    return -(derivative_R+a[:,0]/R+derivative_z)/(4*np.pi*field.G)
source_rows=[]
for q in sorted(set(round(c['q_potential'],7) for c in cases if c['q_potential'] is not None)):
    coarse=source(q,.004);fine=source(q,.002)
    # A negative cell must exceed the step-difference estimate by 10 times.
    stable_negative=(fine<0)&(-fine>10*abs(fine-coarse))
    i=int(np.argmin(fine))
    source_rows.append(dict(q_potential=q,grid_points=len(R),negative_points=int((fine<0).sum()),
        negatives_exceeding_10_step_differences=int(stable_negative.sum()),
        minimum_effective_density_Msun_kpc3=float(fine[i]),minimum_at_R_z_kpc=[float(R[i]),float(z[i])],
        step_change_at_minimum_Msun_kpc3=float(abs(fine[i]-coarse[i])),
        note='Equivalent Newtonian additional source at lambda=1; sign unchanged by positive lambda^p. Signed effective density is not automatically forbidden in modified field equations, but is incompatible with interpreting that source everywhere as positive deposited mass.'))
predictions=[]
for case in cases:
    lam,q=case['lambda_mass'],case['q_potential']
    for i,r in enumerate(stars):predictions.append(dict(case=case['case'],observable='rotation_training',R_kpc=r['R_mean_kpc'],observed=float(seen[i]),predicted=float(radial(lam)[i])))
    for i in range(len(Rv)):predictions.append(dict(case=case['case'],observable='vertical_exposed',R_kpc=float(Rv[i]),observed=float(kv_seen[i]),predicted=float(vertical_prediction(lam,q)[i])))
original=cases[0]['metrics']
assert abs(original['rotation_rms_kms']-18.792085091038476)<1e-9
previous=json.loads((SOURCE/'results-refined.json').read_text())['scores']['vertical/conservative_completion']['rmse']
assert abs(original['vertical_rms_surface_equivalent']-previous)<1e-8
result=dict(classification='Post-exposure hypothetical shape extension using known coordinate deformation; not a derived companion mechanism or fresh holdout',
    formula='Phi_total(R,z)=lambda Phi_b(R,z)+lambda^p Phi_c_original(R,z/q)',
    lambda_scope='All ordinary matter scales together with fixed geometry; not a stellar-only mass-to-light change or observationally justified mass prior.',
    parameters_p=float(P),bounds=dict(lambda_mass=[.5,1.6],q_potential=[.5,4.]),
    extended_bounds=dict(lambda_mass=[.5,1.6],inverse_q=[0.,2.],zero_inverse_q_is_cylindrical_limit=True),
    cases=cases,vertical_nonnegative_addition_floor=floor,source_diagnostic=source_rows,
    optimization_check=dict(extended_grid_points=len(mu_grid)*len(lambda_grid),grid_minimum_log_score=grid_best,extended_optimizer_no_worse_than_grid=True),
    original_scores_reproduced=True,held_out_outcomes_opened=False,
    input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [SOURCE/'run.py',SOURCE/'results-refined.json',SOURCE/'predictions-refined.json',HERE.parent/'cepheid-common-frame/training-bins.json']},
    limitations=['No complete observational likelihood or component mass uncertainty.', 'Vertical estimates are previously exposed model-dependent Bovy-Rix inferences, not new raw-star force measurements.', 'Changing q leaves equatorial extra force unchanged, but changes off-plane geometry.', 'No companion transport/capture law predicts q or lambda.', 'Positive effective density, stable orbital support and finite outer boundary are separate requirements.', 'No relativistic lensing or photon-redshift completion.'])
for name,obj in [('results',result),('predictions',predictions)]:
    (HERE/f'{name}.json').write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2),flush=True)

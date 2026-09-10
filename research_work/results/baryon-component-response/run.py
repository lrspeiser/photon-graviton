"""Separate ordinary components and re-solve the frozen nonlinear field."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, minimize
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=HERE.parent/'conservative-field-completion'
spec=importlib.util.spec_from_file_location('frozen',SOURCE/'run.py')
f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
GROUPS=['stellar_disks','gas_disks','central_stars']
PARAMS=[(1.332e9,2.,.3,2.7,'exp'),(8.97e8,2.8,.9,2.7,'exp'),(5.81e7,7.,.085,4.,'sech2'),(2.68e9,1.5,.045,12.,'sech2')]

class DiskPart(f.Expansion):
    def __init__(self,refined,indices):
        r=np.geomspace(1e-5,500,3072 if refined else 1536)
        ells=np.arange(0,257 if refined else 129,2)
        mu,w=roots_legendre(1024 if refined else 512)
        leg,_=f.basis(ells,mu);angular=leg*(w[:,None]*(2*ells+1)[None,:]/2)
        coeff=np.empty((len(r),len(ells)))
        for start in range(0,len(r),32):
            rr=r[start:start+32,None];R=rr*np.sqrt(1-mu**2);z=rr*mu;rho=np.zeros_like(R)
            for i in indices:
                sigma,rd,h,hole,kind=PARAMS[i];surface=sigma*np.exp(-hole/R-R/rd)
                if kind=='exp':rho+=surface*np.exp(-abs(z)/h)/(2*h)
                else:rho+=surface*np.exp(-2*np.logaddexp(z/(2*h),-z/(2*h))+2*np.log(2))/(4*h)
            coeff[start:start+len(rr)]=rho@angular
        inn,out=f.scaled_integrals(r,coeff*r[:,None],coeff*r[:,None],ells)
        fac=4*np.pi*f.G/(2*ells+1)
        super().__init__(r,ells,-fac*(inn+out),fac*((ells+1)*inn-ells*out)/r[:,None])
        self.mass=float(4*np.pi*inn[-1,0]*r[-1])

class Components:
    def __init__(self,refined):
        self.stellar=DiskPart(refined,[0,1]);self.gas=DiskPart(refined,[2,3])
        cache=ROOT/'research_work/data-cache/bar-field'
        self.bar=f.CachedAxisymmetric(cache/'bar-L64.npz');self.nuclei=f.CachedAxisymmetric(cache/'nuclei-L16.npz')
        self.cache={}
    def parts(self,r,mu):
        r,mu=np.broadcast_arrays(np.atleast_1d(r),np.atleast_1d(mu))
        key=hashlib.sha256(r.tobytes()+mu.tobytes()).digest()
        if key not in self.cache:
            central=np.array(self.bar.evaluate(r,mu))+np.array(self.nuclei.evaluate(r,mu))
            q=r*r+.001**2
            bh=np.array([-f.G*4.1e6/np.sqrt(q),f.G*4.1e6*r/q**1.5,np.zeros_like(r)])
            self.cache[key]=np.array([self.stellar.evaluate(r,mu),self.gas.evaluate(r,mu),central,bh])
        return self.cache[key]

class Model:
    def __init__(self,components,scales):self.components=components;self.scales=np.r_[scales,1.]
    def evaluate(self,r,mu):return np.einsum('i,ijk->jk',self.scales,self.components.parts(r,mu))

catalog=ROOT/'research_work/data-cache/cepheid-stars/cepheid-common-frame-train-selected.parquet'
stars=pd.read_parquet(catalog)
assert set(stars.role)=={'train'}
bins=json.loads((HERE.parent/'cepheid-common-frame/training-bins.json').read_text())
obsR=np.array([r['jeans_proxy_kms'] for r in bins])
old=json.loads((SOURCE/'predictions-refined.json').read_text())
vertical=[r for r in old if r['observable']=='vertical' and r['model']=='ordinary']
Rv=np.array([r['R_kpc'] for r in vertical]);obsZ=np.array([r['observed'] for r in vertical])
obs=np.r_[obsR,obsZ];weight=np.r_[np.full(12,1/np.sqrt(24)),np.full(43,1/np.sqrt(86))]/obs
trace=[]
def solve(parts,scales,refined,label):
    print(label,dict(zip(GROUPS,np.asarray(scales).tolist())),flush=True)
    b=Model(parts,scales);c=f.Completion(b,refined)
    R=stars.R_kpc.to_numpy()
    radial2=-R*(f.force(b,R,0)[:,0]+f.force(c,R,0)[:,0])
    vr=np.array([np.sqrt(radial2[stars.bin.to_numpy()==i].mean()) for i in range(12)])
    vz=-(f.force(b,Rv,1.1)[:,1]+f.force(c,Rv,1.1)[:,1])/(2*np.pi*f.G*1e6)
    pred=np.r_[vr,vz]
    trace.append(dict(label=label,scales=np.asarray(scales).tolist(),refined=refined,
        rotation_rms_kms=float(np.sqrt(np.mean((vr-obsR)**2))),vertical_rms_surface_equivalent=float(np.sqrt(np.mean((vz-obsZ)**2))),
        equal_observable_fractional_mse=float(np.sum((weight*(pred-obs))**2))))
    return pred

print('Building separately normalized coarse components',flush=True)
coarse=Components(False)
reference=f.Baryons(False)
probeR=np.geomspace(.2,30,101);probez=np.linspace(-2,2,101)
combined=f.force(Model(coarse,np.ones(3)),probeR,probez);original=f.force(reference,probeR,probez)
decomposition_error=float(np.max(abs(combined-original))/np.max(abs(original)))
assert decomposition_error<1e-9
baseline=solve(coarse,np.ones(3),False,'coarse_original')
jac=[];delta=.05
for i,name in enumerate(GROUPS):
    low=np.ones(3);high=np.ones(3);low[i]-=delta;high[i]+=delta
    pl=solve(coarse,low,False,name+'_minus_5_percent')
    ph=solve(coarse,high,False,name+'_plus_5_percent')
    jac.append((ph-pl)/(2*delta))
jac=np.array(jac).T
fit=lsq_linear(weight[:,None]*jac,weight*(obs-baseline),bounds=(-.3,.3),tol=1e-12)
assert fit.success
scales=1+fit.x
linear_prediction=baseline+jac@fit.x
candidate=solve(coarse,scales,False,'coarse_linear_response_candidate')
# A second explicit balancing rule prevents improvement in one observable from
# hiding degradation in the other. This still uses the same local response.
baseline_rms=np.array([np.sqrt(np.mean((baseline[:12]-obs[:12])**2)),np.sqrt(np.mean((baseline[12:]-obs[12:])**2))])
def ratios(delta):
    residual=baseline+jac@delta-obs
    return np.array([np.sqrt(np.mean(residual[:12]**2)),np.sqrt(np.mean(residual[12:]**2))])/baseline_rms
balanced=minimize(lambda x:x[3],np.r_[np.zeros(3),1.],method='SLSQP',bounds=[(-.3,.3)]*3+[(0,2)],
    constraints=[{'type':'ineq','fun':lambda x:x[3]-ratios(x[:3])}],options={'ftol':1e-11,'maxiter':300})
assert balanced.success
balanced_scales=1+balanced.x[:3]
balanced_coarse=solve(coarse,balanced_scales,False,'coarse_balanced_candidate')
print('Building refined components and verifying candidate with full field solve',flush=True)
fine=Components(True)
fine_original=solve(fine,np.ones(3),True,'refined_original')
fine_candidate=solve(fine,scales,True,'refined_component_candidate')
fine_balanced=solve(fine,balanced_scales,True,'refined_balanced_candidate')
expected=json.loads((HERE.parent/'cepheid-common-frame/results.json').read_text())['training_proxy_scores']['completion']['rms_kms']
assert abs(next(x for x in trace if x['label']=='refined_original')['rotation_rms_kms']-expected)<1e-7
expected_z=json.loads((SOURCE/'results-refined.json').read_text())['scores']['vertical/conservative_completion']['rmse']
assert abs(next(x for x in trace if x['label']=='refined_original')['vertical_rms_surface_equivalent']-expected_z)<1e-7
responses=[]
for i,name in enumerate(GROUPS):
    responses.append(dict(group=name,median_rotation_change_for_10_percent_kms=float(np.median(.1*jac[:12,i])),
        median_vertical_change_for_10_percent_surface_equivalent=float(np.median(.1*jac[12:,i]))))
rows=[]
for i in range(len(obs)):
    rows.append(dict(observable='rotation_training' if i<12 else 'vertical_exposed',R_kpc=float(bins[i]['R_mean_kpc'] if i<12 else Rv[i-12]),
        observed=float(obs[i]),original=float(fine_original[i]),candidate=float(fine_candidate[i]),balanced=float(fine_balanced[i]),
        derivative_per_unit_component_scale={name:float(jac[i,j]) for j,name in enumerate(GROUPS)}))
result=dict(classification='Post-exposure local component response with two bounded fitting criteria, verified by full nonlinear field solves; not a global optimum, observational mass posterior or fresh validation',
    component_order=GROUPS,candidate_scales=scales.tolist(),exploratory_scale_bounds=[.7,1.3],
    balanced_candidate_scales=balanced_scales.tolist(),balanced_objective='Minimize the larger fractional RMS relative to original RMS across the two observables, using local response; full field checked afterwards.',
    bound_note='Illustrative +/-30 percent sensitivity box; not confidence intervals or observed allowed mass ranges.',
    fixed_items=['All spatial component shapes','Central black-hole mass','Empirical A, p and a_star','Original conservative field equation; no q deformation','Training/holdout assignment'],
    responses=responses,trace=trace,
    checks=dict(component_force_sum_relative_error=decomposition_error,
        max_linearized_vs_full_coarse_rotation_kms=float(np.max(abs(linear_prediction[:12]-candidate[:12]))),
        max_linearized_vs_full_coarse_vertical_units=float(np.max(abs(linear_prediction[12:]-candidate[12:]))),
        max_candidate_coarse_refined_rotation_kms=float(np.max(abs(candidate[:12]-fine_candidate[:12]))),
        max_candidate_coarse_refined_vertical_units=float(np.max(abs(candidate[12:]-fine_candidate[12:]))),
        original_scores_reproduced=True),
    baseline_disk_masses_Msun=dict(stellar=fine.stellar.mass,gas=fine.gas.mass),
    input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in [SOURCE/'run.py',SOURCE/'results-refined.json',catalog,HERE.parent/'cepheid-common-frame/training-bins.json']},
    held_out_outcomes_opened=False,
    limitations=['Only one local response step; not a proof no other component combination can work.', 'Stellar thin/thick disks co-scale, and gas components co-scale; scale lengths/heights fixed.', 'Central stellar bar and nuclei co-scale; bar azimuth averaged.', 'Published vertical inference remains model dependent with its original frame/selection assumptions.', 'No complete likelihood or priors; no full stellar orbit/population fit.', 'No photon production, capture, lensing or energy-supply validation.'])
for name,obj in [('results',result),('predictions',rows)]:
    (HERE/f'{name}.json').write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2),flush=True)

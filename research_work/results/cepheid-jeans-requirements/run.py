"""Post-exposure Jeans requirements and exploratory training-star tilt check."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
SOURCE=HERE.parent/'cepheid-common-frame'
spec=importlib.util.spec_from_file_location('cepheid_reduction',SOURCE/'run.py')
reduction=importlib.util.module_from_spec(spec);spec.loader.exec_module(reduction)
bins=json.loads((SOURCE/'training-bins.json').read_text())
raw=ROOT/'research_work/data-cache/cepheid-stars/cepheid-common-frame-train-selected.parquet'
t=pd.read_parquet(raw)
split=pd.read_parquet(raw.parent/'cepheid-common-frame-split.parquet')
assert set(t.source_id)<=set(split.loc[split.role.eq('train'),'source_id'])
records=[]
for row in bins:
    R=row['R_mean_kpc'];r2=row['corrected_radial_second_moment'];p2=row['corrected_azimuthal_second_moment']
    for model in ['ordinary','completion']:
        v2=row[model+'_vc_kms']**2
        # Known steady axisymmetric radial Jeans equation, with T=R/nu*d_z(nu*q).
        log_pressure_slope=(p2-v2)/r2-1
        density_gradient=log_pressure_slope/R+2/27.3
        needed_tilt=row['jeans_proxy_kms']**2-v2
        reconstructed=p2-r2*(1+log_pressure_slope)
        assert abs(reconstructed-v2)<1e-9
        records.append(dict(bin=row['bin'],model=model,R_kpc=R,n=row['n'],
            minimum_vc_for_nonincreasing_radial_pressure_kms=float(np.sqrt(max(p2-r2,0))),
            model_vc_kms=row[model+'_vc_kms'],
            required_dln_density_times_radial_second_moment_dlnR=log_pressure_slope,
            required_dln_density_dR_if_Rsigma_27_3_per_kpc=density_gradient,
            local_density_efold_growth_kpc=1/density_gradient if density_gradient>0 else None,
            required_tilt_term_kms2=needed_tilt,
            required_midplane_d_vRvz_dz_kms2_per_kpc=needed_tilt/R))

# Propagate the same declared measurement-error scenario into cross moments.
c=reduction.coords(t)
jac=np.stack([reduction.coords(t,delta={k:1.})[:,3:]-c[:,3:] for k in ['pmra','pmdec','rv']],axis=2)
sig=t[['pmra_error','pmdec_error','radial_velocity_error']].to_numpy()
cov=np.zeros((len(t),3,3))
for j in range(3):cov[:,j,j]=sig[:,j]**2
cov[:,0,1]=cov[:,1,0]=sig[:,0]*sig[:,1]*t.pmra_pmdec_corr.to_numpy()
outcov=np.einsum('nij,njk,nlk->nil',jac,cov,jac)
jd=(reduction.coords(t,scale=1.0001)[:,3:]-reduction.coords(t,scale=.9999)[:,3:])/.0002
outcov+=.07**2*np.einsum('ni,nj->nij',jd,jd)
assert np.max(abs(outcov.diagonal(axis1=1,axis2=2)-reduction.variances(t)))<1e-8
q=c[:,3]*c[:,5]-outcov[:,0,2]
z=c[:,2];rng=np.random.default_rng(6120910)
tilt=[]
for lo,hi in [(6,9),(9,12),(12,15),(15,18)]:
    mask=(c[:,0]>=lo)&(c[:,0]<(hi if hi<18 else hi+1e-10))
    zz=z[mask];qq=q[mask];N=len(zz)
    slope=float(zz@qq/(zz@zz))
    leave_one_out=((zz@qq)-zz*qq)/((zz@zz)-zz**2)
    free=np.linalg.lstsq(np.c_[np.ones(N),zz],qq,rcond=None)[0]
    ix=rng.integers(0,N,size=(2048,N));zb=zz[ix];qb=qq[ix]
    slopes=np.sum(zb*qb,axis=1)/np.sum(zb**2,axis=1)
    relevant=[r for r in records if r['model']=='completion' and lo<=r['R_kpc']<hi]
    required=float(np.average([r['required_midplane_d_vRvz_dz_kms2_per_kpc'] for r in relevant],weights=[r['n'] for r in relevant]))
    tilt.append(dict(R_low_kpc=lo,R_high_kpc=hi,n=N,above=int((zz>0).sum()),below=int((zz<0).sum()),
        median_abs_z_kpc=float(np.median(abs(zz))),max_abs_z_kpc=float(max(abs(zz))),
        reflection_linear_slope_kms2_per_kpc=slope,
        leave_one_out_slope_range_kms2_per_kpc=[float(leave_one_out.min()),float(leave_one_out.max())],
        maximum_single_star_slope_leverage=float(max(zz**2)/(zz@zz)),
        bootstrap_slope_percentiles_2_5_16_50_84_97_5=np.percentile(slopes,[2.5,16,50,84,97.5]).tolist(),
        free_intercept_kms2=float(free[0]),free_intercept_slope_kms2_per_kpc=float(free[1]),
        required_slope_for_completion_approx=required,
        qualification='Exploratory regression of corrected vR*vz on z across a finite radial/height range, not a selection-corrected midplane derivative or a formal significance. Bootstrap resamples stars only.'))
primary=[r for r in records if r['model']=='completion']
result=dict(classification='Known Jeans algebra applied to exposed training data; conditional diagnostic, not a new gravity law or holdout',
    training_rows=len(t),bins=len(bins),
    completion_below_zero_tilt_nonincreasing_pressure_floor_bins=sum(r['model_vc_kms']<r['minimum_vc_for_nonincreasing_radial_pressure_kms'] for r in primary),
    required_completion_log_pressure_slope_range=[min(r['required_dln_density_times_radial_second_moment_dlnR'] for r in primary),max(r['required_dln_density_times_radial_second_moment_dlnR'] for r in primary)],
    required_completion_density_efold_kpc_range=[min(r['local_density_efold_growth_kpc'] for r in primary),max(r['local_density_efold_growth_kpc'] for r in primary)],
    required_completion_midplane_tilt_slope_range=[min(r['required_midplane_d_vRvz_dz_kms2_per_kpc'] for r in primary),max(r['required_midplane_d_vRvz_dz_kms2_per_kpc'] for r in primary)],
    new_physics_parameters_fitted=False,validation_or_test_outcomes_opened=False,
    input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE/'training-bins.json',SOURCE/'run.py',raw]},
    limitations=['Pressure floor assumes steady axisymmetry and zero tilt; non-axisymmetry, time dependence, and biased distance/velocity inference can invalidate it.', 'Finite bins summarize midplane quantities approximately.', 'Density gradient is conditional on adopted radial second-moment scale; raw survey counts cannot establish physical density.', 'Tilt fit ignores position error, selection, radial variation within a bin and vertical density gradients away from the plane.', 'No complete astrophysical covariance or population likelihood.'])
for name,obj in [('results',result),('requirements',records),('tilt-training',tilt)]:
    (HERE/f'{name}.json').write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2));print(json.dumps(tilt,indent=2))

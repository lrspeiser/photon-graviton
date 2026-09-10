"""Shared absorption coefficient in a declared isotropic boundary bath."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np
from scipy.special import roots_legendre, logsumexp
from scipy.optimize import minimize_scalar

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache/baryonic-shielding'
CACHE.mkdir(parents=True,exist_ok=True)
CODE=HERE.parent/'baryon-attached-deposits/run.py'
spec=importlib.util.spec_from_file_location('attached',CODE)
attached=importlib.util.module_from_spec(spec);spec.loader.exec_module(attached)
SOURCE=HERE.parent/'baryon-attached-deposits/results.json'
d=json.loads(SOURCE.read_text())
for p,h in d['source_hashes'].items():
    assert attached.digest(ROOT/p)==h
rows=d['rows'];params=d['disk_parameters']
xyz=np.array([[r['R_kpc']*np.cos(r['phi_rad']),r['R_kpc']*np.sin(r['phi_rad']),r['z_kpc']] for r in rows])
target=np.array([r['required_extra_per_baryon'] for r in rows])

def density(x):
    R=np.hypot(x[:,0],x[:,1]);z=x[:,2]
    value=attached.bar(x)+attached.nuclei(x)
    for sigma,rd,h,hole,kind in params:
        surface=sigma*np.exp(-hole/R-R/rd)
        vertical=np.exp(-abs(z)/h)/(2*h) if kind=='exp' else np.exp(-2*np.logaddexp(z/(2*h),-z/(2*h))+2*np.log(2))/(4*h)
        value+=surface*vertical
    value+=3*4.1e6*.001**2/(4*np.pi*(np.sum(x*x,axis=1)+.001**2)**2.5)
    return value

def columns(nmu,nphi,ns,boundary):
    mu,w=roots_legendre(nmu);phi=np.arange(nphi)*2*np.pi/nphi
    m,p=np.meshgrid(mu,phi,indexing='ij')
    directions=np.c_[np.sqrt(1-m.ravel()**2)*np.cos(p.ravel()),np.sqrt(1-m.ravel()**2)*np.sin(p.ravel()),m.ravel()]
    weights=np.repeat(w,nphi)/(2*nphi)
    nodes,sw=roots_legendre(ns)
    result=[]
    for x in xyz:
        dot=directions@x
        length=dot+np.sqrt(dot*dot+boundary**2-x@x)
        # Split at the disk midplane cusp before Gauss integration.
        crossing=np.clip(x[2]/directions[:,2],0,length)
        col=np.zeros(len(directions))
        for left,right in [(np.zeros_like(length),crossing),(crossing,length)]:
            s=left[:,None]+(right-left)[:,None]*(nodes+1)/2
            points=x-s[:,:,None]*directions[:,None,:]
            values=density(points.reshape(-1,3)).reshape(len(directions),ns)
            col+=(values@sw)*(right-left)/2
        assert np.isfinite(col).all() and (col>=0).all()
        result.append(col)
    return np.array(result),weights

def fit(col,w):
    logtarget=np.log(target)
    def evaluate(logk):
        logs=logsumexp(-np.exp(logk)*col+np.log(w),axis=1)
        amplitude=float(np.mean(logtarget-logs))
        residual=logs+amplitude-logtarget
        return float(np.mean(residual**2)),amplitude,logs+amplitude
    grid=np.linspace(np.log(1e-12),np.log(1e-5),81)
    scores=np.array([evaluate(k)[0] for k in grid]);j=int(scores.argmin())
    opt=minimize_scalar(lambda k:evaluate(k)[0],bounds=(grid[max(0,j-1)],grid[min(80,j+1)]),method='bounded')
    assert opt.success
    score,amp,pred=evaluate(opt.x)
    return dict(kappa_kpc2_per_Msun=float(np.exp(opt.x)),log_amplitude=amp,
                rms_log_residual=float(np.sqrt(score)),geometric_scatter_factor=float(np.exp(np.sqrt(score))),
                grid_minimum_at_boundary=j in [0,80],log_prediction=pred.tolist())

def main():
    hashes={str(p.relative_to(ROOT)):attached.digest(p) for p in [Path(__file__),CODE,SOURCE]}
    np.testing.assert_allclose(density(xyz),[r['declared_baryon_density_Msun_kpc3'] for r in rows],rtol=1e-13)
    results={};cache={}
    for label,nmu,nphi,ns,boundary in [('coarse',8,16,64,30.),('fine',16,32,128,30.),('outer60',16,32,128,60.)]:
        col,w=columns(nmu,nphi,ns,boundary)
        path=CACHE/f'{label}.npz';np.savez_compressed(path,column_Msun_kpc2=col,weights=w)
        assert abs(w.sum()-1)<1e-14
        np.testing.assert_allclose(logsumexp(np.log(w)[None,:]+np.zeros_like(col),axis=1),0,atol=1e-14)
        result=fit(col,w);result.update(angular_directions=len(w),radial_nodes_per_segment=ns,boundary_kpc=boundary,cache_sha256=attached.digest(path))
        results[label]=result;cache[label]=(col,w)
        print(label,result['kappa_kpc2_per_Msun'],result['geometric_scatter_factor'],flush=True)
    base=results['coarse'];k=base['kappa_kpc2_per_Msun'];amp=base['log_amplitude']
    frozen={}
    for label,(col,w) in cache.items():
        pred=logsumexp(-k*col+np.log(w),axis=1)+amp
        delta=pred-np.array(base['log_prediction'])
        frozen[label]=dict(max_absolute_log_change=float(abs(delta).max()),median_absolute_log_change=float(np.median(abs(delta))),
                           rms_log_residual=float(np.sqrt(np.mean((pred-np.log(target))**2))))
    for p,h in hashes.items():assert attached.digest(ROOT/p)==h
    result=dict(scope='Synthetic source-shape test, not an observed-star fit or a derived companion capture mechanism.',
        postulate='beta=kappa*rho_b; isotropic uniform external bath; straight rays; fixed baryon geometry; equal duration; attached retained deposits.',
        formula='q_model=C*mean_direction exp(-kappa*integral_ray rho_b ds)',
        samples=len(rows),constant_loading_rms_log=float(np.std(np.log(target))),
        fits=results,frozen_coarse_parameter_checks=frozen,source_hashes=hashes,
        holdouts_opened=False,source_energy_budget_evaluated=False,local_emission_included=False,
        capture_momentum_and_baryon_evolution_included=False)
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

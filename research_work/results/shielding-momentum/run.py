"""Momentum implied by the existing absorption geometry; no new fit."""
from pathlib import Path
import importlib.util
import json
import numpy as np
from scipy.special import roots_legendre, logsumexp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
DRIVER=HERE.parent/'full-bar-orbits/run.py'
spec=importlib.util.spec_from_file_location('momentum_orbit',DRIVER)
orbit=importlib.util.module_from_spec(spec);spec.loader.exec_module(orbit)
SHIELD=HERE.parent/'baryonic-shielding/results.json'
SOURCE=HERE.parent/'baryon-attached-deposits/results.json'
import hashlib
def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    shield=json.loads(SHIELD.read_text());source=json.loads(SOURCE.read_text())
    for p,h in shield['source_hashes'].items():assert digest(ROOT/p)==h
    field=orbit.Field('full')
    files=[Path(__file__),DRIVER,SHIELD,SOURCE,*field.paths]
    hashes={str(p.relative_to(ROOT)):digest(p) for p in files}
    rows=source['rows']
    xyz=np.array([[r['R_kpc']*np.cos(r['phi_rad']),r['R_kpc']*np.sin(r['phi_rad']),r['z_kpc']] for r in rows])
    acc=np.concatenate([field.evaluate(x)[1] for x in np.array_split(xyz,16)])
    gravity=np.linalg.norm(acc,axis=1)
    assert (gravity>0).all()
    k=shield['fits']['fine']['kappa_kpc2_per_Msun']
    C=np.exp(shield['fits']['fine']['log_amplitude'])
    outputs={}
    for label,nmu,nphi in [('coarse',8,16),('fine',16,32),('outer60',16,32)]:
        path=ROOT/f'research_work/data-cache/baryonic-shielding/{label}.npz'
        assert digest(path)==shield['fits'][label]['cache_sha256']
        hashes[str(path.relative_to(ROOT))]=digest(path)
        with np.load(path) as f:col,w=f['column_Msun_kpc2'],f['weights']
        mu,gw=roots_legendre(nmu);phi=np.arange(nphi)*2*np.pi/nphi
        m,p=np.meshgrid(mu,phi,indexing='ij')
        directions=np.c_[np.sqrt(1-m.ravel()**2)*np.cos(p.ravel()),np.sqrt(1-m.ravel()**2)*np.sin(p.ravel()),m.ravel()]
        np.testing.assert_allclose(w,np.repeat(gw,nphi)/(2*nphi),rtol=1e-13)
        np.testing.assert_allclose(w@directions,0,atol=1e-14)
        logits=np.log(w)[None,:]-k*col
        total=logsumexp(logits,axis=1)
        delivered=np.exp(logits-total[:,None])
        xi=delivered@directions
        norm=np.linalg.norm(xi,axis=1)
        assert (norm<=1+1e-12).all()
        q=C*np.exp(total)
        # Constant loading over T, cold Newtonian attached-source interpretation:
        # P/(M_current c^2)=q/[(1+q)T]. q is the MODEL loading, not the target.
        # a_capture = c*xi*q/[(1+q)T] at leading order for a stationary receiver.
        # Acceleration units (km/s)^2/kpc yield T in kpc/(km/s).
        T=299792.458*norm*q/(1+q)/gravity
        gyr=T*3.085677581491367e16/(365.25*86400*1e9)
        alignment=np.sum(xi*acc,axis=1)/np.maximum(norm*gravity,1e-300)
        output=[]
        for i,r in enumerate(rows):
            output.append(dict(R_kpc=r['R_kpc'],z_kpc=r['z_kpc'],phi_rad=r['phi_rad'],
                predicted_loading=float(q[i]),momentum_per_energy_over_c=xi[i].tolist(),
                anisotropy_magnitude=float(norm[i]),alignment_with_gravity=float(alignment[i]),
                gravity_kms2_per_kpc=float(gravity[i]),critical_loading_duration_Gyr=float(gyr[i]),
                duration_for_force_below_10_percent_Gyr=float(10*gyr[i])))
        outputs[label]=output
    fine=outputs['fine']
    xi=np.array([r['momentum_per_energy_over_c'] for r in fine])
    comparisons={}
    for label in ['coarse','outer60']:
        other=np.array([r['momentum_per_energy_over_c'] for r in outputs[label]])
        comparisons[label]=dict(max_absolute_anisotropy_vector_change=float(np.linalg.norm(other-xi,axis=1).max()))
    # Independent optically thin uniform-sphere control: xi/beta -> -x/3.
    # The even square-root part of the boundary length has zero vector moment.
    test_x=np.array([2.,1.,.5]);dot=directions@test_x
    length=dot+np.sqrt(dot*dot+30**2-test_x@test_x)
    control=[]
    for beta in [1e-5,1e-6,1e-7]:
        prob=w*np.exp(-beta*length);vector=prob@directions/prob.sum()
        control.append(float(np.linalg.norm(vector/beta+test_x/3)))
    assert control[-1]<1e-5 and control[-1]<control[0]
    for p,h in hashes.items():assert digest(ROOT/p)==h
    result=dict(scope='Conditional leading-order momentum consequence of the previously fitted shielding model.',
        assumptions=['Companion momentum is E/c in its travel direction.',
                     'Complete absorption, no compensating outgoing momentum.',
                     'Stationary-receiver leading-order force; motion-dependent drag omitted.',
                     'Inertial added mass equals effective Newtonian deposited source.',
                     'Constant loading over unknown T in the fixed snapshot geometry.'],
        loading_time_is_inferred=False,cosmic_age_assumed=False,holdouts_opened=False,
        source_energy_budget_evaluated=False,rows=outputs,quadrature_comparisons=comparisons,
        fitted_kappa_fixed=k,normalization_fixed=float(C),source_hashes=hashes,
        uniform_sphere_linear_limit_errors=control,
        reference_gravity='Frozen full-bar target field, not gravity recomputed from the mismatching shielding source.')
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print('anisotropy',min(r['anisotropy_magnitude'] for r in fine),max(r['anisotropy_magnitude'] for r in fine))
    print('critical Gyr',min(r['critical_loading_duration_Gyr'] for r in fine),max(r['critical_loading_duration_Gyr'] for r in fine))
    print(comparisons)

if __name__=='__main__':main()

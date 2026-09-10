"""Orbit readiness for the conservative full-bar empirical response."""
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FULL=HERE.parent/'full-bar-completion'
spec=importlib.util.spec_from_file_location('full_bar_completion',FULL/'run.py')
full=importlib.util.module_from_spec(spec);spec.loader.exec_module(full)
from field import point_potential
CACHE=ROOT/'research_work/data-cache/full-bar-orbits'
CACHE.mkdir(parents=True,exist_ok=True)
SOURCE=HERE.parent/'rotating-bar-orbits/results.json'
Y0=np.array(json.loads(SOURCE.read_text())['initial_conditions'])
TIMES=np.linspace(0,.25,501)
OMEGA=37.5


def cross(v):
    return np.c_[-v[:,1],v[:,0],np.zeros(len(v))]


def rotate(v,angle):
    c,s=np.cos(angle),np.sin(angle)
    out=v.copy();out[:,0]=c*v[:,0]-s*v[:,1];out[:,1]=s*v[:,0]+c*v[:,1]
    return out


class Field:
    def __init__(self,kind):
        _,self.disk_path=full.disk_model()
        f=np.load(self.disk_path)
        disk=full.FastMultipole(f['r'],f['phi']*np.sqrt(4*np.pi/(2*f['ell']+1)),f['ell'],np.zeros_like(f['ell']))
        barpath=ROOT/'research_work/data-cache/bar-field/bar-L64.npz'
        nucpath=ROOT/'research_work/data-cache/bar-field/nuclei-L16.npz'
        nuc=np.load(nucpath);m0=nuc['m']==0
        nuclei=full.FastMultipole(nuc['r'],nuc['coeff'][:,m0],nuc['l'][m0],nuc['m'][m0])
        self.parts=[disk,full.load_field(barpath),nuclei]
        self.paths=[self.disk_path,barpath,nucpath]
        if kind!='ordinary':
            ap=full.CACHE/'axisymmetric-reference-80-L128-N512.npz'
            self.parts.append(full.load_field(ap));self.paths.append(ap)
        if kind in ['full','lower_field']:
            cp=full.CACHE/('finer.npz' if kind=='full' else 'fine.npz')
            self.parts.append(full.load_field(cp));self.paths.append(cp)
        self.bounds=[float('inf'),0.,0.]

    def evaluate(self,x):
        rad=np.linalg.norm(x,axis=1);R=np.hypot(x[:,0],x[:,1])
        self.bounds[0]=min(self.bounds[0],float(rad.min()))
        self.bounds[1]=max(self.bounds[1],float(rad.max()))
        self.bounds[2]=max(self.bounds[2],float(abs(x[:,2]).max()))
        # Each component checks its own interpolation domain. No extrapolation.
        p,a=point_potential(x)
        for field in self.parts:
            q,b=field.evaluate(x);p+=q;a+=b
        return p,a


def integrate(field,tol,inertial=False,times=TIMES):
    def rhs(t,y):
        s=y.reshape(-1,6);x,p=s[:,:3],s[:,3:]
        if inertial:
            _,a=field.evaluate(rotate(x,-OMEGA*t))
            return np.c_[p,rotate(a,OMEGA*t)].ravel()
        _,a=field.evaluate(x)
        return np.c_[p-OMEGA*cross(x),a-OMEGA*cross(p)].ravel()
    sol=solve_ivp(rhs,[times[0],times[-1]],Y0.ravel(),t_eval=times,
                  method='DOP853',rtol=tol,atol=tol*.01)
    assert sol.success
    return sol.y.T.reshape(len(times),len(Y0),6),sol.nfev


def summary(field,states):
    x,p=states[:,:,:3],states[:,:,3:]
    potentials=np.array([field.evaluate(q)[0] for q in x])
    L=x[:,:,0]*p[:,:,1]-x[:,:,1]*p[:,:,0]
    J=.5*np.sum(p*p,axis=2)+potentials-OMEGA*L
    R=np.hypot(x[:,:,0],x[:,:,1]);z=x[:,:,2]
    return [dict(probe=j,R_min_kpc=float(R[:,j].min()),R_max_kpc=float(R[:,j].max()),
                 max_abs_height_kpc=float(abs(z[:,j]).max()),
                 mean_rotation_kms=float(np.mean(L[:,j]/R[:,j])),
                 sampled_plane_fraction=float(np.mean(abs(z[:,j])<.2)),
                 sampled_offplane_fraction=float(np.mean((abs(z[:,j])>.5)&(abs(z[:,j])<1.5))),
                 Jacobi_drift_over_220_squared=float(np.max(abs(J[:,j]-J[0,j]))/220**2)) for j in range(len(Y0))]


def main():
    results={};trajectories={};inputs=[Path(__file__),SOURCE,FULL/'run.py']
    for kind in ['ordinary','axisymmetric_extra','full']:
        field=Field(kind);inputs+=field.paths
        attempts=[];last=None;retained=None
        for tol in [2e-9,2e-11,2e-13]:
            start=time.monotonic();states,nfev=integrate(field,tol)
            item=dict(rtol=tol,nfev=nfev,seconds=time.monotonic()-start)
            if last is not None:
                dx=float(np.max(np.linalg.norm(states[:,:,:3]-last[:,:,:3],axis=2)))
                dv=float(np.max(np.linalg.norm(states[:,:,3:]-last[:,:,3:],axis=2)))
                item.update(position_difference_kpc=dx,velocity_difference_kms=dv,
                            integration_refinement_pass=bool(dx<1e-4 and dv<.01))
            attempts.append(item)
            print(kind,json.dumps(item),flush=True)
            retained=states
            if item.get('integration_refinement_pass',False):break
            last=states
        rows=summary(field,retained)
        np.savez_compressed(CACHE/f'{kind}.npz',times=TIMES,trajectory=retained)
        trajectories[kind]=retained
        results[kind]=dict(attempts=attempts,probes=rows,all_evaluation_domain_bounds=field.bounds,
                           Jacobi_gate_pass=max(r['Jacobi_drift_over_220_squared'] for r in rows)<1e-5,
                           integration_refinement_pass=attempts[-1].get('integration_refinement_pass',False))
        (HERE/'partial-results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8',newline='\n')
    field=Field('full');lower=Field('lower_field');inputs+=lower.paths
    points=trajectories['full'][:,:,:3].reshape(-1,3)
    # Compare extra force only, normalizing away neither a dominant ordinary
    # field nor an additive potential gauge.
    old=Field('ordinary')
    errors=[]
    for i in range(0,len(points),32):
        q=points[i:i+32];base=old.evaluate(q)[1]
        fine=field.evaluate(q)[1]-base;coarse=lower.evaluate(q)[1]-base
        errors.extend((np.linalg.norm(fine-coarse,axis=1)/np.linalg.norm(fine,axis=1)).tolist())
    # An independent inertial-frame integration tests rotating-frame signs.
    times=np.linspace(0,.05,101)
    rotating,_=integrate(field,2e-13,times=times)
    inertial,_=integrate(field,2e-13,inertial=True,times=times)
    converted=np.array([np.c_[rotate(s[:,:3],-OMEGA*t),rotate(s[:,3:],-OMEGA*t)] for t,s in zip(times,inertial)])
    frame_dx=float(np.max(np.linalg.norm(converted[:,:,:3]-rotating[:,:,:3],axis=2)))
    frame_dv=float(np.max(np.linalg.norm(converted[:,:,3:]-rotating[:,:,3:],axis=2)))
    result=dict(scope='Six synthetic orbit probes; not a fitted stellar population or data likelihood.',
        known_equations='H_J=|p|²/2+Phi-Omega*Lz; xdot=p-Omega cross x; pdot=-grad(Phi)-Omega cross p',
        initial_conditions=Y0.tolist(),pattern_speed_kms_per_kpc=OMEGA,duration_kpc_per_kms=.25,
        models=results,actual_trajectory_extra_force_refinement_max=float(max(errors)),
        actual_trajectory_extra_force_refinement_median=float(np.median(errors)),
        actual_trajectory_field_gate_pass=bool(max(errors)<.01),
        frame_check_duration_kpc_per_kms=.05,inertial_rotating_position_difference_kpc=frame_dx,
        inertial_rotating_velocity_difference_kms=frame_dv,frame_check_pass=bool(frame_dx<1e-4 and frame_dv<.01),
        input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in set(inputs)},
        holdouts_opened=False,stellar_population_fit=False,photon_origin_derived=False)
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Trajectory field refinement max',max(errors),'frame differences',frame_dx,frame_dv,flush=True)


if __name__=='__main__':main()

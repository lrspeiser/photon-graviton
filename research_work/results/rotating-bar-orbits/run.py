"""Bar-frame orbits under a conservative potential; no observed stars fitted."""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
FOUNDATION=HERE.parent/'bar-field-foundation'
sys.path.insert(0,str(FOUNDATION))
from field import Multipole,disks,disk_evaluate,point_potential,halo,companion
from fast_multipole import FastMultipole

CACHE=ROOT/'research_work/data-cache/rotating-bar-orbits'
CACHE.mkdir(parents=True,exist_ok=True)
FIELD_CACHE=ROOT/'research_work/data-cache/bar-field'
PROTOCOL=json.loads((HERE/'protocol.json').read_text())
OMEGA=PROTOCOL['bar_pattern_speed_kms_per_kpc']
TIMES=np.linspace(0,PROTOCOL['duration_kpc_per_kms'],501)


def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def digest(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def rotate(v,angle):
    v=np.asarray(v);out=v.copy();c=np.cos(angle);s=np.sin(angle)
    out[...,0]=c*v[...,0]-s*v[...,1]
    out[...,1]=s*v[...,0]+c*v[...,1]
    return out


def cross_z(v):
    return np.stack([-v[...,1],v[...,0],np.zeros_like(v[...,0])],axis=-1)


class DiskGrid:
    def __init__(self,nr,nz,potential):
        # Fingerprint source implementation and grid parameters, never silently
        # reuse a table built from changed disk code.
        source_hash=digest(FOUNDATION/'field.py')
        self.path=CACHE/f'disk-N24-Rmin0.001-{nr}-{nz}-{source_hash[:12]}.npz'
        if not self.path.exists():
            lr=np.linspace(np.log(.001),np.log(40.),nr)
            zz=np.linspace(0,np.arcsinh(10/.05),nz)
            vals=np.empty((nr,nz));started=time.monotonic()
            for i,r in enumerate(np.exp(lr)):
                xyz=np.c_[np.full(nz,r),np.zeros(nz),.05*np.sinh(zz)]
                # Only potential samples are stored; all accelerations below
                # are derivatives of that same spline.
                vals[i]=[potential(r/8,z/8,use_physical=False)*220**2 for z in xyz[:,2]]
                if i%32==0:print('disk grid',nr,nz,'row',i,'seconds',round(time.monotonic()-started,1),flush=True)
            assert np.isfinite(vals).all()
            np.savez_compressed(self.path,lr=lr,zz=zz,values=vals,source_hash=source_hash)
        f=np.load(self.path)
        assert str(f['source_hash'])==source_hash
        self.lr=f['lr'];zp=f['zz'];v=f['values']
        self.zz=np.r_[-zp[:0:-1],zp]
        self.spline=RectBivariateSpline(self.lr,self.zz,np.c_[v[:, :0:-1],v],kx=3,ky=3,s=0)

    def evaluate(self,xyz):
        xyz=np.atleast_2d(xyz);x,y,z=xyz.T;R=np.hypot(x,y)
        lr=np.log(R);zz=np.arcsinh(z/.05)
        if np.any((lr<self.lr[0])|(lr>self.lr[-1])|(abs(zz)>self.zz[-1])):
            raise ValueError(f'Orbit leaves disk domain: R range {R.min(),R.max()}, max abs z {abs(z).max()}')
        pot=self.spline.ev(lr,zz)
        ar=-self.spline.ev(lr,zz,dx=1)/R
        az=-self.spline.ev(lr,zz,dy=1)/np.sqrt(z*z+.05**2)
        return pot,np.c_[ar*x/R,ar*y/R,az]


class Field:
    def __init__(self,disk,model='ordinary_matter',order=40):
        self.disk=disk;self.model=model
        self.bar=FastMultipole.load(FIELD_CACHE/f'bar-L{order}.npz')
        self.nuclei=FastMultipole.load(FIELD_CACHE/'nuclei-L16.npz')

    def evaluate(self,x):
        p,a=self.bar.evaluate(x)
        for pair in [self.nuclei.evaluate(x),self.disk.evaluate(x),point_potential(x)]:
            p+=pair[0];a+=pair[1]
        if self.model=='halo_comparison':
            q,b=halo(x);p+=q;a+=b
        elif self.model.startswith('companion_'):
            q,b=companion(x,self.model.split('_')[1]);p+=q;a+=b
        return p,a


def integrate(field,y0,tol,times=TIMES,inertial=False):
    def rhs(t,y):
        state=y.reshape(-1,6);x=state[:,:3];p=state[:,3:]
        if inertial:
            _,a=field.evaluate(rotate(x,-OMEGA*t))
            return np.c_[p,rotate(a,OMEGA*t)].ravel()
        _,a=field.evaluate(x)
        return np.c_[p-OMEGA*cross_z(x),a-OMEGA*cross_z(p)].ravel()
    out=solve_ivp(rhs,(times[0],times[-1]),y0.ravel(),method='DOP853',
                  rtol=tol,atol=tol*.01,t_eval=times)
    assert out.success
    return out.y.T.reshape(len(times),-1,6),out.nfev


def summarize(field,trajectory):
    pots=np.array([field.evaluate(s[:,:3])[0] for s in trajectory])
    x=trajectory[:,:,:3];p=trajectory[:,:,3:]
    L=x[:,:,0]*p[:,:,1]-x[:,:,1]*p[:,:,0]
    J=.5*np.sum(p*p,axis=2)+pots-OMEGA*L
    R=np.hypot(x[:,:,0],x[:,:,1]);z=x[:,:,2]
    vphi=L/R
    records=[]
    for j in range(x.shape[1]):
        records.append(dict(probe=j,R_min_kpc=float(R[:,j].min()),R_max_kpc=float(R[:,j].max()),
                            max_abs_z_kpc=float(abs(z[:,j]).max()),mean_vphi_kms=float(vphi[:,j].mean()),
                            sampled_plane_fraction=float(np.mean(abs(z[:,j])<.2)),
                            sampled_offplane_fraction=float(np.mean((abs(z[:,j])>.5)&(abs(z[:,j])<1.5))),
                            max_Jacobi_drift_over_220_squared=float(np.max(abs(J[:,j]-J[0,j]))/220**2)))
    return records


def main():
    started=time.monotonic();d=disks(24)
    coarse=DiskGrid(128,65,d);fine=DiskGrid(256,129,d)
    rng=np.random.default_rng(73051)
    R=np.exp(rng.uniform(np.log(.2),np.log(25),80))
    z=rng.uniform(-2,2,80)
    points=np.c_[R,np.zeros(len(R)),z]
    pp,aa=disk_evaluate(d,points)
    checks={}
    for name,g in [('coarse',coarse),('fine',fine)]:
        p,a=g.evaluate(points)
        rel=np.linalg.norm(a-aa,axis=1)/np.linalg.norm(aa,axis=1)
        checks[name]=dict(max_force_relative_error=float(rel.max()),median_force_relative_error=float(np.median(rel)),
                          worst_point_kpc=points[rel.argmax()].tolist())
    save('disk-interpolation-checks.json',checks)
    print('disk interpolation checks',checks,flush=True)
    assert checks['fine']['max_force_relative_error']<PROTOCOL['acceptance_checks']['disk_interpolation_force_error_fraction_of_reference_disk_force']
    baseline=Field(fine)
    x=np.array([[r*np.cos(np.pi/6),r*np.sin(np.pi/6),z] for r in [1,3,8] for z in [.1,.8]])
    _,a=baseline.evaluate(x);R=np.hypot(x[:,0],x[:,1]);er=x[:,:2]/R[:,None]
    ar=np.sum(a[:,:2]*er,axis=1);vc=np.sqrt(-R*ar)
    p=np.c_[20*er[:,0]-.8*vc*er[:,1],20*er[:,1]+.8*vc*er[:,0],np.full(len(x),20)]
    y0=np.c_[x,p]
    results={};trajectories={}
    for model in PROTOCOL['models']:
        f=Field(fine,model)
        low,nlow=integrate(f,y0,2e-9)
        high,nhigh=integrate(f,y0,2e-11)
        records=summarize(f,high)
        dx=float(np.max(np.linalg.norm(high[:,:,:3]-low[:,:,:3],axis=2)))
        dv=float(np.max(np.linalg.norm(high[:,:,3:]-low[:,:,3:],axis=2)))
        first_refinement=dict(rtol_pair=[2e-9,2e-11],position_difference_kpc=dx,velocity_difference_kms=dv)
        used_tolerance=2e-11
        if dx>1e-4 or dv>.01:
            refined,nref=integrate(f,y0,2e-13)
            dx=float(np.max(np.linalg.norm(refined[:,:,:3]-high[:,:,:3],axis=2)))
            dv=float(np.max(np.linalg.norm(refined[:,:,3:]-high[:,:,3:],axis=2)))
            high=refined;nhigh=nref;used_tolerance=2e-13
            records=summarize(f,high)
        results[model]=dict(probes=records,tolerance_position_difference_kpc=dx,
                            tolerance_velocity_difference_kms=dv,function_evaluations=[nlow,nhigh],
                            original_refinement=first_refinement,retained_rtol=used_tolerance)
        trajectories[model]=high
        np.savez_compressed(CACHE/(model+'-trajectory.npz'),times=TIMES,trajectory=high)
        (CACHE/'partial-results.json').write_text(json.dumps(results,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
        print('completed',model,'seconds',round(time.monotonic()-started,1),'dx',dx,'dv',dv,flush=True)
    # Separate inertial-coordinate dynamics, transformed back for comparison.
    times=TIMES[TIMES<=.05]
    inertial,neval=integrate(baseline,y0,2e-13,times,inertial=True)
    rotating=trajectories['ordinary_matter'][:len(times)]
    inertial[:,:,:3]=rotate(inertial[:,:,:3],-OMEGA*times[:,None])
    inertial[:,:,3:]=rotate(inertial[:,:,3:],-OMEGA*times[:,None])
    frame_error=float(np.max(np.linalg.norm(inertial[:,:,:3]-rotating[:,:,:3],axis=2)))
    # Existing potential-order uncertainty must remain distinct from integrator error.
    lower=Field(fine,order=24)
    lower_traj,_=integrate(lower,y0,2e-13)
    lower_summary=summarize(lower,lower_traj)
    results['bar_order_comparison']=dict(L24_probes=lower_summary,
        max_trajectory_position_difference_kpc=float(np.max(np.linalg.norm(lower_traj[:,:,:3]-trajectories['ordinary_matter'][:,:,:3],axis=2))),
        status='Potential truncation sensitivity, not observational uncertainty or convergence proof')
    results['verification']=dict(inertial_rotating_max_position_difference_kpc=frame_error,
        elapsed_seconds=time.monotonic()-started,observed_stars_fitted=False,heldout_scores_evaluated=False)
    results['initial_conditions']=y0.tolist()
    results['cache_provenance']=[dict(file=str(path.relative_to(ROOT)),sha256=digest(path))
         for path in [fine.path,coarse.path,FIELD_CACHE/'bar-L24.npz',FIELD_CACHE/'bar-L40.npz',FIELD_CACHE/'nuclei-L16.npz']]
    # Preserve outputs even if a declared numerical gate fails.
    save('results.json',results)
    np.savez_compressed(CACHE/'trajectories.npz',times=TIMES,**trajectories,L24=lower_traj)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for j,ax in enumerate(axes.ravel()):
        for name in ['ordinary_matter','companion_equatorial','companion_caps']:
            t=trajectories[name][:,j]
            ax.plot(np.hypot(t[:,0],t[:,1]),t[:,2],lw=.8,label=name.replace('_',' '))
        ax.set(xlabel='Radius from rotation axis [kpc]',ylabel='Height above disk [kpc]',title=f'Synthetic orbit {j+1}')
    axes[0,0].legend(fontsize=7)
    fig.suptitle('Illustrative rotating-bar trajectories — not fitted observed stars')
    fig.savefig(HERE/'orbit-comparison.png',dpi=150);plt.close(fig)
    gates=PROTOCOL['acceptance_checks']
    assert frame_error<gates['max_inertial_rotating_coordinate_difference_kpc']
    for model in PROTOCOL['models']:
        r=results[model]
        assert r['tolerance_position_difference_kpc']<gates['max_tolerance_refinement_position_difference_kpc']
        assert r['tolerance_velocity_difference_kms']<gates['max_tolerance_refinement_velocity_difference_kms']
        assert max(p['max_Jacobi_drift_over_220_squared'] for p in r['probes'])<gates['max_Jacobi_drift_divided_by_220_squared']
    print('Declared orbit integration gates passed.',flush=True)


if __name__=='__main__':main()

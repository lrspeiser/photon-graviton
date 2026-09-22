"""Independent JR-5 numerical checks and frozen-coefficient source sensitivities."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.special import ellipk,ellipe
from core import AxisModel,gas_density,transfer,baseline,G,C
from data_screen import prepare,predict_one,score


def direct_equatorial_force(R,alpha,beta,kind,n=280):
    """Independent spherical-volume integral with analytic azimuth kernel.
    No Legendre multipoles. Gauss-Legendre log-radius intervals split at R.
    """
    x,w=np.polynomial.legendre.leggauss(n)
    rr=[];ww=[]
    for lo,hi in [(np.log(1e-5),np.log(R)),(np.log(R),np.log(1e4))]:
        r=np.exp((lo+hi)/2+(hi-lo)*x/2);rr.extend(r);ww.extend(w*(hi-lo)/2*r)
    r=np.array(rr);wr=np.array(ww)
    x,w=np.polynomial.legendre.leggauss(n);mu=(x+1)/2;wm=w/2
    edges=np.r_[0.,r];mid=(edges[:-1]+edges[1:])/2
    f=transfer(edges,gas_density(mid,mu,1e9,2.,.1,kind),alpha,beta)
    rho0=baseline(r,1,1,20,2)[1];delta=rho0[:,None]*(2*f-1)
    s=r[:,None]*np.sqrt(1-mu[None,:]**2);z=r[:,None]*mu[None,:]
    den2=(R+s)**2+z*z;den=np.sqrt(den2);u=4*R*s/den2
    one_minus=((R-s)**2+z*z)/den2
    K=ellipk(u);E=ellipe(u)
    kp=np.where(u>1e-7,(E/one_minus-K)/(2*np.maximum(u,1e-300)),np.pi/8*(1+9*u/8))
    du=4*s/den2-8*R*s*(R+s)/den2**2
    kernel=(R+s)*K/den**3-kp*du/den
    delta_g=8*G*np.sum(wr[:,None]*r[:,None]**2*wm[None,:]*delta*kernel)
    return float(baseline(np.array([R]),1,1,20,2)[0][0]+delta_g)


def projection_check(model,b=2.):
    # Direct cylinder mass by angular quadrature of the untruncated density.
    # Linear interpolation of the saved angular flux, not its multipole fit.
    x,w=np.polynomial.legendre.leggauss(128);r=model.r
    mum=np.sqrt(np.maximum(0.,1-(b/r)**2));integ=[]
    for i,lower in enumerate(mum):
        mu=lower+(1-lower)*(x+1)/2
        f=np.interp(mu,model.mu,model.f[i])
        integ.append(np.sum(w*f)*(1-lower)) # factor 2*f times half weights
    mass=4*np.pi*np.trapezoid(r*r*model.rho0*np.array(integ),r)
    direct=4*G*mass/(C*C*b)
    return dict(direct_cylinder_bend=float(direct),multipole_ray_bend=model.bend(b,'face'),
        relative_difference=float(direct/model.bend(b,'face')-1))


def main():
    p=argparse.ArgumentParser(__doc__);p.add_argument('--jr1',type=Path,required=True);p.add_argument('--results',type=Path,required=True);a=p.parse_args()
    target=a.results/'audit.json'
    if target.exists():raise FileExistsError(target)
    summary=json.loads((a.results/'data-summary.json').read_text());sel=summary['selection']
    profiles=json.loads((a.results/'gas-profiles.json').read_text())
    coarse,_=prepare(a.jr1,profiles=profiles);fine,_=prepare(a.jr1,nr=512,nmu=64,lmax=24,profiles=profiles)
    refinement=[]
    for gc,gf in zip(coarse,fine):
        c=predict_one(gc,sel['alpha'],sel['beta'],sel['quadratic']);f=predict_one(gf,sel['alpha'],sel['beta'],sel['quadratic'])
        refinement.append(dict(name=gc['name'],max_velocity_relative=float(np.max(abs(c/f-1))),max_velocity_kms=float(np.max(abs(c-f)))))
    finer_score=score(fine,sel['alpha'],sel['beta'],sel['quadratic'])
    (a.results/'fine-selected-full.json').write_text(json.dumps(finer_score,indent=2)+'\n')
    sensitivity=[]
    for h in [.05,.2]:
        gg,_=prepare(a.jr1,nr=512,nmu=64,lmax=24,height=h,profiles=profiles)
        res=score(gg,sel['alpha'],sel['beta'],sel['quadratic'])
        (a.results/f'height-{h}-full.json').write_text(json.dumps(res,indent=2)+'\n')
        sensitivity.append(dict(height_over_gas_scale=h,groups=res['groups']))
    direct=[]
    for kind in ['disk','ring']:
        model=AxisModel(1,1,20,2,nr=1536,nmu=192,lmax=60).solve(1e9,2,.3,.3,kind=kind)
        forces=[]
        for n in [160,320]:forces.append(direct_equatorial_force(2,.3,.3,kind,n))
        direct.append(dict(kind=kind,force_direct=forces,force_multipole=float(model.radial(np.array([2.]))[0]),
            direct_self_refinement=forces[0]/forces[1]-1,
            relative_force_difference=forces[-1]/float(model.radial(np.array([2.]))[0])-1,
            lens_projection=projection_check(model)))
    # Independent readback calculations; no optimizer calls.
    checks=[]
    for fn in ['null-full.json','selected-full.json','active-full.json','monopole-control.json','fine-selected-full.json']:
        record=json.loads((a.results/fn).read_text())
        for row in record['rows']:
            pred=np.array(row['predicted_kms']);y=np.array(row['observed_kms']);er=np.array(row['error_kms'])
            checks.extend([abs(np.sqrt(np.mean((pred-y)**2))-row['RMS_kms'])<1e-10,
             abs(np.sqrt(np.mean(((pred-y)/y)**2))-row['fractional_RMS'])<1e-12,
             abs(np.sum(((pred-y)/er)**2)-row['chi2'])<1e-7])
    stats=dict(checks=len(checks),failed=int(np.count_nonzero(np.logical_not(checks))))
    hashes={str(p.relative_to(a.jr1)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
        a.jr1/'evidence/R10_spheroid_stellar_population.json',
        a.jr1/'repository/temporal_candidate_audit/data/Rotmod_LTG.zip',
        a.jr1/'repository/temporal_candidate_audit/data/SPARC_Lelli2016c.mrt']}
    result=dict(refinement=refinement,maximum_velocity_relative=max(x['max_velocity_relative'] for x in refinement),
        above_one_percent=sum(x['max_velocity_relative']>.01 for x in refinement),fine_groups=finer_score['groups'],
        thickness_sensitivity=sensitivity,independent_toy_forces=direct,readback=stats,input_sha256=hashes,
        interpretation='Arithmetic and finite-resolution checks, not confirmation of physical model or source reconstructions')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()

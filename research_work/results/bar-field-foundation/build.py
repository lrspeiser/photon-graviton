"""Numerical foundation check before any stellar holdout fit."""
from pathlib import Path
import json,time
import numpy as np
from field import Multipole,bar,nuclei,point_potential,halo,companion,disks,disk_evaluate,G,H

CACHE=H.parents[2]/'research_work/data-cache/bar-field';CACHE.mkdir(parents=True,exist_ok=True)
save=lambda name,obj:(H/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

def cache(name,density,**kwargs):
    p=CACHE/(name+'.npz')
    if p.exists():return Multipole.load(p)
    t=time.time();v=Multipole.build(density,**kwargs);v.save(p)
    print(name,'built',round(time.time()-t,2),'s; mass',v.integrated_mass,flush=True);return v

if __name__=='__main__':
    # A known spherical Plummer model gives an independent analytic answer.
    mass=1e10;soft=1.
    pl=cache('plummer-check',lambda xyz:3*mass/(4*np.pi*soft**3)*(1+np.sum(xyz*xyz,axis=1)/soft**2)**(-2.5),lmax=0,nr=1200,ntheta=20,nphi=16)
    R=np.geomspace(.05,30,60);points=np.c_[R,R*.3,R*.4]
    pp,aa=pl.evaluate(points);exactp,exacta=point_potential(points,mass,soft)
    errors={'plummer_max_potential_fraction':float(np.max(abs(pp/exactp-1))),
            'plummer_max_force_fraction':float(np.max(np.linalg.norm(aa-exacta,axis=1)/np.linalg.norm(exacta,axis=1)))}
    assert errors['plummer_max_force_fraction']<.002
    low=cache('bar-L16',bar,lmax=16,nr=600,ntheta=64,nphi=96)
    high=cache('bar-L24',bar,lmax=24,nr=900,ntheta=96,nphi=144)
    refined=cache('bar-L40',bar,lmax=40,nr=1200,ntheta=144,nphi=192)
    nucleus=cache('nuclei-L16',nuclei,lmax=16,nr=900,ntheta=72,nphi=32)
    # Positions chosen before evaluating stellar velocity likelihoods.
    xyz=[]
    for R in [.5,1,2,3,5,8,12,20]:
        for z in [.1,.5,1.1]:
            for phi in [0,np.pi/4,np.pi/2]:xyz.append([R*np.cos(phi),R*np.sin(phi),z])
    xyz=np.array(xyz)
    pl,al=low.evaluate(xyz);ph,ah=high.evaluate(xyz)
    rel=np.linalg.norm(ah-al,axis=1)/np.maximum(np.linalg.norm(ah,axis=1),1)
    errors['bar_L16_to_L24_max_force_fraction']=float(rel.max())
    errors['bar_L16_to_L24_median_force_fraction']=float(np.median(rel))
    errors['bar_mass_low_Msun']=low.integrated_mass;errors['bar_mass_high_Msun']=high.integrated_mass
    pr,ar=refined.evaluate(xyz)
    refined_error=np.linalg.norm(ar-ah,axis=1)/np.maximum(np.linalg.norm(ar,axis=1),1)
    errors['bar_L24_to_L40_max_force_fraction']=float(refined_error.max())
    errors['bar_L24_to_L40_median_force_fraction']=float(np.median(refined_error))
    errors['bar_mass_refined_Msun']=refined.integrated_mass
    errors['nuclear_mass_Msun']=nucleus.integrated_mass
    print('Bar convergence',errors,flush=True)
    p0=xyz[[2,17,36]]
    _,a=high.evaluate(p0);numeric=[];step=1e-5
    for j in range(3):
        shift=np.eye(3)[j]*step
        numeric.append(-(high.evaluate(p0+shift)[0]-high.evaluate(p0-shift)[0])/(2*step))
    errors['bar_gradient_max_abs_error_kms2_per_kpc']=float(np.max(abs(np.array(numeric).T-a)))
    assert errors['bar_gradient_max_abs_error_kms2_per_kpc']<.02
    t=time.time();disk=disks(24);print('Disks built in',time.time()-t,flush=True)
    dp,da=disk_evaluate(disk,xyz)
    npot,na=nucleus.evaluate(xyz);bp,ba=point_potential(xyz)
    baselinep=pr+dp+npot+bp;baselinea=ar+da+na+ba
    errors['bar_resolution_change_fraction_of_total_force']=float(np.max(np.linalg.norm(ar-ah,axis=1)/np.linalg.norm(baselinea,axis=1)))
    hp,ha=halo(xyz)
    rows=[]
    for name in ['ordinary_matter','halo_comparison','companion_equatorial','companion_caps','companion_shell']:
        p=baselinep.copy();a=baselinea.copy()
        if name=='halo_comparison':p+=hp;a+=ha
        elif name.startswith('companion_'):
            cp,ca=companion(xyz,name.split('_')[1]);p+=cp;a+=ca
        for i,pos in enumerate(xyz):
            rows.append({'model':name,'x_kpc':float(pos[0]),'y_kpc':float(pos[1]),'z_kpc':float(pos[2]),
                         'potential_kms2':float(p[i]),'acceleration_kms2_per_kpc':a[i].tolist()})
    save('field-predictions.json',rows)
    save('numerical-checks.json',errors)
    print(json.dumps(errors,indent=2),flush=True)

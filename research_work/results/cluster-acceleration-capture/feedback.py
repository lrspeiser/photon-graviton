"""Conditional immobile-storage feedback, including its gravity in capture."""
from pathlib import Path
import json
import numpy as np
from scipy.special import roots_legendre
OUT=Path(__file__).resolve().parent

def evolve(n,step,outer=1000.,angular=48):
    edges=np.r_[0.,np.geomspace(.001,outer,n)]
    r=(.5*(edges[1:]**3+edges[:-1]**3))**(1/3)
    volume=4*np.pi/3*np.diff(edges**3)
    mu,w=roots_legendre(angular)
    impact2=r[:,None]**2*(1-mu[None,:]**2);z=r[:,None]*mu[None,:]
    root=np.sqrt(np.maximum(0.,edges[None,None,:]**2-impact2[:,:,None]))
    inside=np.clip(z[:,:,None]+root,0,2*root)
    lengths=np.diff(inside,axis=2)
    mass_b=r**3/(1+r*r)**1.5
    # Independent full-chord capture integration, split at shell boundaries.
    bg,wg=roots_legendre(12)
    b=(.5*(edges[:-1,None]+edges[1:,None])+.5*np.diff(edges)[:,None]*bg).ravel()
    bw=(.5*np.diff(edges)[:,None]*wg).ravel()
    chord_lengths=2*np.diff(np.sqrt(np.maximum(0.,edges[None,:]**2-b[:,None]**2)),axis=1)
    def cross_section(mass):
        enclosed=np.cumsum(mass)-.5*mass
        opacity=((mass_b+enclosed)/(r*r))**2
        return float(np.sum(2*np.pi*b*bw*(-np.expm1(-np.einsum('ij,j->i',chord_lengths,opacity)))))
    def rate(mass):
        enclosed=np.cumsum(mass)-.5*mass
        g=(mass_b+enclosed)/(r*r)
        opacity=g*g # A=1, G=Mb=a=1; common coefficient unchanged.
        tau=np.einsum('ijk,k->ij',lengths,opacity,optimize=False)
        q=opacity*(np.exp(-tau)@w)/2
        mass_rate=q*volume
        return mass_rate,mass_rate.sum()
    zero=np.zeros(n);initial,total0=rate(zero)
    # Evolve by total added mass m; dt/dm=1/total capture power.
    mass=zero.copy();t=0.;snap=[]
    targets=[.1,.5,1.]
    for target in targets:
        while mass.sum()<target-1e-12:
            dm=min(step,target-mass.sum())
            one,p1=rate(mass);trial=mass+dm*one/p1
            two,p2=rate(trial)
            mass+=dm*.5*(one/p1+two/p2)
            t+=dm*.5*(1/p1+1/p2)
        current,power=rate(mass)
        intercepted=cross_section(mass)
        assert abs(power/intercepted-1)<.02
        cumulative=np.cumsum(mass)
        median=float(np.interp(target/2,np.r_[0,cumulative],edges))
        snap.append(dict(deposit_mass_over_baryonic_mass=target,dimensionless_exposure=t,capture_power=power,
           capture_power_over_initial=power/total0,half_stored_mass_radius_over_a=median,
           stored_mass_inside_a=float(np.interp(1.,edges,np.r_[0,cumulative])),ledger_error=float(mass.sum()-target),
           independent_chord_capture=intercepted,volume_vs_chord_relative=power/intercepted-1))
    return dict(radial_shells=n,mass_step=step,outer_boundary_over_a=outer,angular_nodes=angular,initial_capture_power=total0,snapshots=snap)

runs=[evolve(128,.02),evolve(256,.02),evolve(256,.01)]
reference=json.loads((OUT/'results.json').read_text())['cases'][1]['effective_area_over_a2']
assert abs(runs[-1]['initial_capture_power']/reference-1)<.02
for run in runs:
    assert all(abs(x['ledger_error'])<1e-10 for x in run['snapshots'])
    assert all(x['capture_power']>run['initial_capture_power'] for x in run['snapshots'])
differences=[]
for a,b,c in zip(runs[0]['snapshots'],runs[1]['snapshots'],runs[2]['snapshots']):
    differences.append(dict(mass=c['deposit_mass_over_baryonic_mass'],radial_power_relative=b['capture_power']/a['capture_power']-1,
        step_power_relative=c['capture_power']/b['capture_power']-1,step_exposure_relative=c['dimensionless_exposure']/b['dimensionless_exposure']-1))
assert max(abs(v['radial_power_relative']) for v in differences)<.03
assert max(abs(v['step_power_relative']) for v in differences)<.002
boundary=evolve(256,.01,outer=2000.,angular=96)
result=dict(scope='Conditional immobile stored mass feedback in spherical Newtonian geometry, A=1, imposed constant incoming companion bath. No support stresses, orbital motion, field-energy budget or global bath depletion. Not a complete conserving relativistic evolution.',
    physical_time_unit='Mb*c/(u_infinity*a^2); dimensionless time = physical time / this unit',runs=runs,refinement=differences,
    combined_boundary_angular_check=boundary,
    combined_boundary_angular_final_power_relative=boundary['snapshots'][-1]['capture_power']/runs[-1]['snapshots'][-1]['capture_power']-1)
(OUT/'feedback-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(dict(finest=runs[-1],refinement=differences),indent=2),flush=True)

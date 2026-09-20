"""CC-2W preregistered linear exterior and lensing calculation."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import subprocess
import numpy as np
from scipy.integrate import quad
from scipy.special import k1
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_audit import HERE,ROOT,save,digest


def main():
    out=HERE/'weak-field-v1';out.mkdir(exist_ok=False)
    sources=[HERE/'weak_field.py',HERE/'weak-field-protocol.md',HERE/'run_audit.py']
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                 hashes={p.relative_to(ROOT).as_posix():digest(p) for p in sources}))
    c=.08**2*6/(4*np.pi);radii=np.geomspace(1,100,101)
    rows=[];lenses=[]
    for omega in (0,.02,.2,1):
        potential=lambda r:-c*np.exp(-omega*r)/r
        acceleration=lambda r:c*np.exp(-omega*r)*(1+omega*r)/r**2
        velocity=lambda r:np.sqrt(r*acceleration(r))
        for r in radii:
            x=omega*r;h=r*1e-4
            force=acceleration(r);speed=velocity(r)
            slope=-.5-x*x/(2*(1+x))
            epicycle=c*np.exp(-x)*(1+x-x*x)/r**3
            force_fd=(potential(r+h)-potential(r-h))/(2*h)
            slope_fd=(np.log(velocity(r+h))-np.log(velocity(r-h)))/(np.log(r+h)-np.log(r-h))
            epicycle_fd=((r+h)**3*acceleration(r+h)-(r-h)**3*acceleration(r-h))/(2*h*r**3)
            errors=[abs(force-force_fd),abs(slope-slope_fd)/max(1,abs(slope)),abs(epicycle-epicycle_fd)]
            rows.append(dict(omega=omega,r=r,acceleration=force,speed=speed,log_speed_slope=slope,
                             epicyclic_frequency_squared=epicycle,flat_target=bool(abs(slope)<=.1),
                             stable_circular_orbit=bool(epicycle>0),numerical_errors=errors,numerical_pass=bool(max(errors)<=1e-7)))
        for impact in np.geomspace(1,10,12):
            def integrand(theta):
                co=np.cos(theta)
                return 2*c/impact*np.exp(-omega*impact/co)*(co+omega*impact)
            bend,error=quad(integrand,0,np.pi/2,epsabs=1e-14,epsrel=1e-11)
            reference=2*c/impact if omega==0 else 2*c*omega*k1(omega*impact)
            relative=abs(bend-reference)/abs(reference)
            lenses.append(dict(omega=omega,impact=impact,bend=bend,reference=reference,quadrature_error=error,
                               relative_error=relative,numerical_pass=bool(relative<=1e-7),
                               ratio_to_unscreened_spatial_curvature_benchmark=bend/(4*c/impact)))
    theta=np.arange(2048)*2*np.pi/2048
    positions=np.column_stack((np.cos(theta),np.sin(theta),np.zeros_like(theta)))
    currents=.2*6/2048*np.column_stack((-np.sin(theta),np.cos(theta),np.zeros_like(theta)))
    ring=[]
    for r in (10,20,40,80):
        distance=np.linalg.norm(np.array([r,0,0])-positions,axis=1)
        field=-.5*.08*np.sum(currents/distance[:,None],axis=0)/(4*np.pi)
        ring.append(dict(r=r,field=field,amplitude=np.linalg.norm(field)))
    slopes=np.diff(np.log([r['amplitude'] for r in ring]))/np.diff(np.log([r['r'] for r in ring]))
    summary=dict(orbit_cases=len(rows),orbit_numerical_passes=sum(r['numerical_pass'] for r in rows),
                 flat_target_passes=sum(r['flat_target'] for r in rows),
                 unstable_orbit_counts={str(w):sum(not r['stable_circular_orbit'] for r in rows if r['omega']==w) for w in (0,.02,.2,1)},
                 lens_cases=len(lenses),lens_numerical_passes=sum(r['numerical_pass'] for r in lenses),
                 maximum_lens_reference_error=max(r['relative_error'] for r in lenses),
                 ring_total_current=currents.sum(axis=0),ring_amplitude_slopes=slopes,
                 massless_scalar_lensing_ratio=.5)
    save(out/'orbits.json',rows);save(out/'lensing.json',lenses);save(out/'ring.json',ring);save(out/'summary.json',summary)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for omega in (0,.02,.2,1):
        values=[r for r in rows if r['omega']==omega]
        axes[0].loglog(radii,[r['speed'] for r in values],label=f'omega={omega:g}')
        lens=[r for r in lenses if r['omega']==omega]
        axes[1].semilogx([r['impact'] for r in lens],[r['ratio_to_unscreened_spatial_curvature_benchmark'] for r in lens],label=f'omega={omega:g}')
    axes[0].set(xlabel='Radius (model units)',ylabel='Circular speed',title='Linear exterior: no flat outer branch',ylim=(1e-6,.1))
    axes[0].legend();axes[1].axhline(1,color='black',ls='--',label='Spatial-curvature benchmark')
    axes[1].set(xlabel='Impact parameter',ylabel='Candidate deflection / (4 C / b)',title='Static lapse-only light response',ylim=(0,1.1));axes[1].legend(fontsize=8)
    fig.suptitle('CC-2W structural screen — no observed data fitted')
    fig.savefig(HERE/'weak-field.png',dpi=170);plt.close(fig)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(encoding='utf8'))


if __name__=='__main__':main()

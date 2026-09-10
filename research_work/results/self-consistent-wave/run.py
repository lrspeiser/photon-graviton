"""Stationary Schrodinger-Poisson source in a fixed Hernquist ordinary galaxy."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_bvp, cumulative_trapezoid, simpson, quad
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent

def solve(eta,f,L=60.,tol=2e-7,previous=None,eps=1e-3):
    x=np.geomspace(eps,L,1500)
    width=min(max(eta**(1/3),eta)/(1+f/3),L/8)
    u=x*np.exp(-x/width)
    u*=np.sqrt(f/simpson(u*u,x=x))
    M=cumulative_trapezoid(u*u,x,initial=0)
    potential=-f/L-(cumulative_trapezoid(M/x**2,x,initial=0)[-1]-cumulative_trapezoid(M/x**2,x,initial=0))
    y=np.array([u,u*(1/x-1/width),potential,M])
    energy=-.3-f/(2*width)
    if previous is not None:
        y=previous.sol(x);y[:2]*=np.sqrt(f/previous.f);y[2:]*=f/previous.f
        outside=x>previous.x[-1]
        y[0:2,outside]=0; y[2,outside]=-f/x[outside]; y[3,outside]=f
        energy=previous.p[0]
    def rhs(x,y,e):
        return np.array([y[1],2/eta*(-1/(1+x)+y[2]-e[0])*y[0],y[3]/x**2,y[0]**2])
    def bc(a,b,e):
        return np.array([a[0]-eps*a[1],a[3]-eps*a[0]**2/3,b[0],b[3]-f,b[2]+f/L])
    sol=solve_bvp(rhs,bc,x,y,p=[energy],tol=tol,max_nodes=45000)
    sol.f=f
    if not sol.success: raise RuntimeError(f'{sol.message}; max residual={max(sol.rms_residuals)} at x={sol.x[np.argmax(sol.rms_residuals)]}')
    xx=np.geomspace(eps,L*.95,2000);uu=sol.sol(xx)[0]
    if np.min(uu)<-1e-7*np.max(np.abs(uu)) and np.max(uu)>1e-7*np.max(np.abs(uu)): raise RuntimeError('Nodal/excited solution')
    return sol

def diagnostics(sol,eta,f,L):
    x=np.geomspace(sol.x[0],L,12000);u,up,v,M=sol.sol(x)
    kinetic=eta/2*simpson((up-u/x)**2,x=x)
    self_energy=.5*simpson(u*u*v,x=x)
    bary_energy=simpson(-u*u/(1+x),x=x)
    external_virial=simpson(u*u*x/(1+x)**2,x=x)
    virial=2*kinetic+self_energy-external_virial
    half=brentq(lambda r:sol.sol(r)[3]-f/2,1e-3,L)
    norm=simpson(u*u,x=x)
    return {'eta':eta,'source_to_baryon_mass_ratio':f,'domain_over_a':L,
            'inner_radius_over_a':float(sol.x[0]),'eigenvalue':float(sol.p[0]),'half_mass_radius_over_a':half,
            'normalization_relative_error':abs(norm/f-1),
            'virial_relative_residual':abs(virial)/(2*kinetic+abs(self_energy)+external_virial),
            'kinetic_energy':kinetic,'self_gravitational_energy':self_energy,
            'external_potential_energy':bary_energy,
            'mass_fraction_beyond_0p8_domain':float(1-sol.sol(.8*L)[3]/f),
            'max_solver_rms_residual':float(np.max(sol.rms_residuals))}

def main():
    records=[]; profiles=[]

    for eta in [.3,1.,3.]:
        prev=None
        for f in sorted(set(np.geomspace(.03,3.,65).tolist()+[.1,1.,3.])):
            sol=solve(eta,f,previous=prev);prev=sol
            if f not in [.1,1.,3.]:continue
            base=diagnostics(sol,eta,f,60.)
            # Wider domain and tighter tolerance; interpolate prior solution only as a guess.
            refined=solve(eta,f,L=90.,tol=2e-9,previous=sol,eps=5e-4)
            fine=diagnostics(refined,eta,f,90.)
            radii=np.geomspace(.05,30.,150)
            m0=sol.sol(radii)[3];m1=refined.sol(radii)[3]
            # Scale by total source mass to avoid meaningless relative error in tiny central masses.
            err=float(np.max(np.abs(m0-m1))/f)
            fine['coarse_to_refined_max_enclosed_mass_difference_over_total']=err
            fine['coarse_to_refined_half_mass_radius_relative_change']=abs(fine['half_mass_radius_over_a']/base['half_mass_radius_over_a']-1)
            assert err<2e-5 and fine['virial_relative_residual']<2e-5, (err,fine)
            assert fine['normalization_relative_error']<1e-6
            lens=[]
            for b in [.5,1.,2.,5.]:
                def mass_at(r): return float(refined.sol(r)[3]) if r<90 else f
                bending=quad(lambda t:mass_at(b/np.cos(t))*np.cos(t)/b,0,np.pi/2,epsabs=1e-9)[0]
                projected=mass_at(b)+quad(lambda r:refined.sol(r)[0]**2*(b/r)**2/(1+np.sqrt(1-(b/r)**2)),b,90,epsabs=1e-9)[0]
                assert abs(bending*b/projected-1)<1e-6
                lens.append({'impact_over_a':b,'extra_bending_in_4GMb_over_ac2':bending,
                             'projected_mass_identity_relative_error':abs(bending*b/projected-1)})
            fine['conditional_weak_lensing']=lens
            records.append(fine)
            for r in radii:
                u,up,v,M=refined.sol(r)
                profiles.append({'eta':eta,'mass_ratio':f,'r_over_a':float(r),
                    'enclosed_source_mass_over_baryon_mass':float(M),
                    'density_in_Mb_over_4pi_a3':float((u/r)**2),
                    'extra_acceleration_in_GMb_over_a2':float(M/r**2),
                    'total_circular_speed_squared_in_GMb_over_a':float(r/(1+r)**2+M/r)})
            print(json.dumps(fine),flush=True)
    (HERE/'results.json').write_text(json.dumps({'classification':'Numerical stationary source equilibria; not capture, stability or observed validation',
        'equations':'Free single-state Schrodinger-Poisson in fixed Hernquist ordinary potential',
        'cases':records},indent=2)+'\n',newline='\n')
    (HERE/'profiles.json').write_text(json.dumps(profiles,indent=2)+'\n',newline='\n')


if __name__ == "__main__":
    main()

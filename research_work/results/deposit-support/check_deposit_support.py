"""Conditional support requirements for a photon-supplied reservoir in a galaxy."""
from pathlib import Path
import json,os
import numpy as np
from scipy.constants import c,parsec
from scipy.integrate import solve_ivp,quad
from scipy.special import erf

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'deposit-support'

def main():
    # Unit r0=vc=1; log potential used only from r=.1 to 1.
    def stop(t,y):return y[0]-.1
    stop.terminal=True;stop.direction=-1
    fall=solve_ivp(lambda t,y:[y[1],-1/y[0]],(0,2),[1.,0.],events=stop,
                   rtol=1e-11,atol=1e-13,max_step=.01)
    expected=np.sqrt(np.pi/2)*erf(np.sqrt(np.log(10)))
    assert fall.success and len(fall.t_events[0])==1
    t=float(fall.t_events[0][0])
    energy=fall.y[1]**2/2+np.log(fall.y[0])
    assert abs(t-expected)<1e-9 and np.max(abs(energy))<1e-8
    jeans=[]
    radii=np.geomspace(1,.01,201)
    for beta in [-1.,0.,.5,.9,1.]:
        exponent=2*(1-beta)
        analytic=-np.log(radii) if beta==1 else -np.expm1(exponent*np.log(radii))/exponent
        sol=solve_ivp(lambda r,y:[(-1+exponent*y[0])/r],(1,.01),[0.],
                       t_eval=radii,rtol=1e-11,atol=1e-13)
        error=float(np.max(abs(sol.y[0]-analytic)))
        assert sol.success and error<1e-8 and np.min(sol.y[0])>=-1e-12
        jeans.append({'anisotropy_beta':beta,'maximum_Jeans_solution_error':error,
                      'sigma_r_squared_over_vc_squared_at_r_over_R_01':float(analytic[-1]),
                      'scope':'Moment equation with zero radial stress at outer boundary; distribution positivity and stability not established.'})

    # Positive isotropic distribution f=exp[-(v^2/2+Phi)/s^2] in fixed log Phi.
    # Test its density/moments by quadrature and the collisionless equation directly.
    s2=.5
    velocity_integral=4*np.pi*quad(lambda v:v*v*np.exp(-v*v/(2*s2)),0,np.inf)[0]
    second=4*np.pi*quad(lambda v:v**4*np.exp(-v*v/(2*s2)),0,np.inf)[0]/velocity_integral
    assert abs(velocity_integral-(2*np.pi*s2)**1.5)<1e-10
    assert abs(second-3*s2)<1e-10
    rng=np.random.default_rng(9026);errors=[]
    for _ in range(100):
        x=rng.normal(size=3);v=rng.normal(size=3);radius=np.linalg.norm(x)
        gradient=x/radius**2
        f=np.exp(-(np.dot(v,v)/2+np.log(radius))/s2)
        grad_x_f=-f*gradient/s2;grad_v_f=-f*v/s2
        errors.append(float(abs(np.dot(v,grad_x_f)-np.dot(gradient,grad_v_f))))
    assert max(errors)<1e-10
    # Tangential orbital support: circular orbit in the same fixed potential.
    def orbit_rhs(t,y):
        pos=y[:2];return np.r_[y[2:],-pos/np.dot(pos,pos)]
    orbit=solve_ivp(orbit_rhs,(0,4*np.pi),[1.,0.,0.,1.],rtol=1e-11,atol=1e-13,max_step=.02)
    radius=np.linalg.norm(orbit.y[:2],axis=0)
    angular=orbit.y[0]*orbit.y[3]-orbit.y[1]*orbit.y[2]
    assert orbit.success and np.max(abs(radius-1))<1e-8 and np.max(abs(angular-1))<1e-8
    vc=200000.;r0=10e3*parsec;year=365.25*86400
    radiation_slope=4*(vc/c)**2
    values={'scope':'Fixed-potential mechanical benchmarks; no new material population, energy supply, full halo formation, stability or astronomical fit.',
            'cold_fall':{'radius_initial_kpc':10.,'radius_final_kpc':1.,'circular_speed_km_s':200.,
                         'fall_time_Myr':t*r0/vc/year/1e6,'dimensionless_ODE_time':t,
                         'dimensionless_analytic_time':float(expected),'max_specific_energy_error_in_vc2_units':float(np.max(abs(energy)))},
            'finite_boundary_Jeans':jeans,
            'positive_isotropic_distribution':{'sigma_r_over_vc':float(np.sqrt(s2)),
                         'max_collisionless_equation_residual':max(errors),'velocity_second_moment_over_vc2':float(second),
                         'kinetic_energy_over_rest_energy_at_200_km_s':float(3*vc*vc/(4*c*c)),
                         'limitations':'Scale-free equilibrium has singular center and infinite extent/mass if extrapolated. Not a finite, growing, self-consistent galaxy.'},
            'circular_support':{'max_radius_error':float(np.max(abs(radius-1))),
                                'max_specific_angular_momentum_error':float(np.max(abs(angular-1)))},
            'isotropic_radiation':{'assumptions':'Static isotropic radiation with P=u/3, ordinary stress conservation and specified lapse; not collisionless streaming or a confined solution.',
                                  'required_density_log_slope_at_200_km_s':float(radiation_slope),
                                  'energy_density_ratio_r10_to_r1':float(np.exp(-radiation_slope*np.log(10))),
                                  'target_r_minus_2_ratio_r10_to_r1':.01},
            'next_requirements':['Receiver or field identity and phase-space injection source',
                                 'Angular momentum, velocity dispersion or supporting field stresses',
                                 'Finite core, outer boundary and escape',
                                 'Time-dependent self-gravity with energy-source depletion',
                                 'Joint motion/lensing and stability against perturbations']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'deposit-support-results.json').write_text(json.dumps(values,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'cold_fall':values['cold_fall'],'radiation':values['isotropic_radiation']},indent=2))

if __name__=='__main__':main()

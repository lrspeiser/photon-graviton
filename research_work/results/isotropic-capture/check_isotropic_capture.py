"""Direction-weighted absorption by a moving receiver and an attenuating sphere."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.special import roots_legendre

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'isotropic-capture'

def moments(beta,slope):
    gamma=1/np.sqrt(1-beta*beta)
    power=gamma**slope*quad(lambda mu:(1-beta*mu)**(slope+1),-1,1,epsabs=1e-12)[0]/2
    force=gamma**slope*quad(lambda mu:mu*(1-beta*mu)**(slope+1),-1,1,epsabs=1e-12)[0]/2
    return power,force

def main():
    constant=[]
    for beta in [0.,.001,.3,.8]:
        power,force=moments(beta,0.)
        assert abs(power-1)<1e-12 and abs(force+beta/3)<1e-12
        gamma=1/np.sqrt(1-beta*beta)
        uprime=quad(lambda mu:(gamma*(1-beta*mu))**2,-1,1)[0]/2
        fluxprime=quad(lambda mu:gamma**2*(1-beta*mu)*(mu-beta),-1,1)[0]/2
        assert abs(uprime-gamma**2*(1+beta*beta/3))<1e-12
        assert abs(fluxprime+4*gamma**2*beta/3)<1e-12
        # Inverse Lorentz transform of the receiver-frame absorbed four-force,
        # followed by d proper time / d lab time, must equal lab angular integrals.
        assert abs(uprime+beta*fluxprime-power)<1e-12
        assert abs(fluxprime+beta*uprime-force)<1e-12
        constant.append({'beta':beta,'power_over_c_sigma_u':power,
                         'force_over_sigma_u':force,'receiver_energy_density_over_u':uprime,
                         'receiver_flux_over_c_u':fluxprime})
    evolution=[]
    for beta in [.001,.3,.8]:
        gamma=1/np.sqrt(1-beta*beta);e0=gamma;p0=gamma*beta
        def rhs(t,y):
            energy,momentum,mass_integral=y
            v=momentum/energy;mass=np.sqrt(energy*energy-momentum*momentum)
            return [1.,-v/3,energy/mass*(1+v*v/3)]
        times=np.linspace(0,e0,101)
        sol=solve_ivp(rhs,(0,e0),[e0,p0,1.],t_eval=times,rtol=1e-11,atol=1e-13)
        energies=e0+times;momenta=p0*(energies/e0)**(-1/3)
        masses=np.sqrt(energies**2-momenta**2)
        error=float(np.max(abs(sol.y-np.array([energies,momenta,masses]))))
        assert sol.success and error<1e-9
        evolution.append({'initial_beta':beta,'energy_doubling_final_beta':float(sol.y[1,-1]/sol.y[0,-1]),
                          'speed_ratio':float(sol.y[1,-1]/sol.y[0,-1]/beta),
                          'analytic_speed_ratio':float(2**(-4/3)),
                          'max_exact_solution_error':error})
    spectral=[]
    for slope in [-6.,-4.,-1.,0.,2.,6.]:
        beta=1e-4;power,force=moments(beta,slope)
        coefficient=1-force/(power*beta)
        expected=(slope+4)/3
        assert abs(coefficient-expected)<1e-6
        spectral.append({'cross_section_local_energy_power':slope,'small_speed_deceleration_coefficient':coefficient,
                         'analytic_coefficient':expected})
    cancellation=[]
    for beta in [.001,.3,.8]:
        power,force=moments(beta,-4.)
        assert abs(power-1)<1e-10 and abs(force-beta)<1e-10
        cancellation.append({'beta':beta,'power':power,'force':force,
                             'velocity_acceleration_numerator':float(force-beta*power)})
    sphere=[]
    mu,w=roots_legendre(512)
    for tau in [.01,1.,10.]:
        for radius in [0.,.5,.9,1.]:
            distance=radius*mu+np.sqrt(1-radius**2+radius**2*mu*mu)
            intensity=np.exp(-tau*distance)
            angular=float(intensity@w);flux=float((mu*intensity)@w)
            xi=flux/angular
            assert xi<=1e-12 and abs(xi)<=1
            if radius==0:assert abs(xi)<1e-14
            if radius==1:
                exact0=1-np.expm1(-2*tau)/(2*tau)
                exact1=-.5+(1-(1+2*tau)*np.exp(-2*tau))/(4*tau*tau)
                assert abs(xi-exact1/exact0)<1e-4
            sphere.append({'optical_depth':tau,'radius_over_R':radius,
                           'radial_delivered_momentum_over_energy_div_c':xi,
                           'energy_density_over_unattenuated_bath':angular/2})
    result={'scope':'Prescribed radiation and absorption cross sections. Derived direction weighting, not a physical material, complete capture law or self-consistent galaxy.',
            'constant_cross_section':constant,'free_receiver_exact_evolution':evolution,
            'frequency_selective_comparisons':spectral,'inverse_fourth_power_cancellation':cancellation,
            'attenuated_external_sphere':sphere,
            'adiabatic_mass_doubling_radius_ratios':{'frequency_independent_isotropic_bath_flat_speed':2**(-4/3),
                                                   'frequency_independent_isotropic_bath_Kepler':2**(-8/3),
                                                   'scope':'Slow, nearly circular loading, fixed dominant potential; not a full orbit integration.'},
            'unresolved':['Actual receiver spectrum and finite bandwidth',
                          'Companion illumination and sources',
                          'All emitted/reverse channels, storage capacity and thermal consistency',
                          'Coupled transport and orbit evolution with changing gravity',
                          'Formation, stability and joint lensing']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'isotropic-capture-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'orbital_ratios':result['adiabatic_mass_doubling_radius_ratios'],
                      'spectral_comparisons':spectral},indent=2))

if __name__=='__main__':main()

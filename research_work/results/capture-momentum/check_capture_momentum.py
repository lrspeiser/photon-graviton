"""Exact local absorption kinematics and conditional galactic loading constraints."""
from pathlib import Path
import json,os
import numpy as np
from scipy.constants import c,parsec

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'capture-momentum'

def absorb(mass,velocity,energy,direction,anisotropy):
    # c=1. The packet bundle has four-momentum (energy, xi*energy*n).
    # xi=0 is balanced delivered radiation, not automatically an isotropic bath.
    gamma=1/np.sqrt(1-np.dot(velocity,velocity))
    before=np.r_[gamma*mass,gamma*mass*velocity]
    delivered=np.r_[energy,anisotropy*energy*direction]
    after=before+delivered
    final_mass=np.sqrt(after[0]**2-np.dot(after[1:],after[1:]))
    return before,after,final_mass,after[1:]/after[0]

def main():
    rng=np.random.default_rng(9092026);errors=[]
    for _ in range(200):
        direction=rng.normal(size=3);direction/=np.linalg.norm(direction)
        velocity=rng.normal(size=3);velocity*=rng.uniform(0,.8)/np.linalg.norm(velocity)
        mass=rng.uniform(.5,2);energy=mass*10**rng.uniform(-5,1);xi=rng.uniform(-1,1)
        before,after,mf,vf=absorb(mass,velocity,energy,direction,xi)
        reconstructed=np.r_[mf,mf*vf]/np.sqrt(1-np.dot(vf,vf))
        errors.append(float(np.max(abs(reconstructed-after))))
        expected=mass**2+2*energy*(before[0]-xi*np.dot(before[1:],direction))+energy**2*(1-xi*xi)
        assert np.isclose(mf*mf,expected,rtol=1e-11)
    assert max(errors)<1e-9
    examples=[]
    velocity=np.array([0.,200000/c,0.]);direction=np.array([1.,0.,0.])
    for xi in [0.,1.]:
        for ratio in [1e-6,1e-3,1.]:
            before,after,mf,vf=absorb(1.,velocity,ratio,direction,xi)
            assert after[2]==before[2] # angular momentum at r=(1,0,0)
            rest=mf-1;kinetic=(after[0]-mf)-(before[0]-1)
            assert abs(rest+kinetic-ratio)<1e-12
            examples.append({'delivered_anisotropy':xi,'incident_energy_over_initial_rest_energy':ratio,
                             'final_rest_mass_over_initial':float(mf),'radial_speed_km_s':float(vf[0]*c/1000),
                             'tangential_speed_km_s':float(vf[1]*c/1000),
                             'rest_energy_gain_over_incident':float(rest/ratio),
                             'kinetic_energy_change_over_incident':float(kinetic/ratio)})
    # Infinitesimal four-momentum update vs independently differentiated formulas.
    velocity=np.array([.2,.3,0.]);gamma=1/np.sqrt(1-np.dot(velocity,velocity))
    derivative=[]
    for xi in [-1.,0.,.5,1.]:
        epsilon=1e-8
        _,_,mf,vf=absorb(1.,velocity,epsilon,direction,xi)
        predicted_v=(xi*direction-velocity)/gamma
        predicted_m=gamma*(1-xi*np.dot(velocity,direction))
        ev=float(np.max(abs((vf-velocity)/epsilon-predicted_v)))
        em=float(abs((mf-1)/epsilon-predicted_m))
        assert ev<1e-7 and em<1e-7
        derivative.append({'xi':xi,'velocity_derivative_error':ev,'mass_derivative_error':em})
    r=10000*parsec;vc=200000.;year=365.25*86400;duration=1e10*year
    growth=np.log(2)/duration
    pressure_timescale=c*r/vc**2
    ratio=growth*pressure_timescale
    result={'scope':'Local special-relativistic capture plus weak-field instantaneous force and adiabatic-orbit comparisons; no actual radiation field, capture rate, self-consistent growing galaxy or astronomical fit.',
            'random_cases':200,'maximum_four_momentum_reconstruction_error':max(errors),
            'packet_examples':examples,'differential_checks':derivative,
            'galactic_example':{'radius_kpc':10.,'circular_speed_km_s':200.,
                                'assumed_receiver_mass_doubling_time_Gyr':10.,
                                'outward_force_over_gravity_for_xi_1':float(ratio),
                                'absolute_xi_for_radiation_force_below_gravity':float(1/ratio),
                                'absolute_xi_for_radiation_force_below_tenth_gravity':float(.1/ratio),
                                'mass_growth_efold_time_for_equal_force_at_xi_1_Gyr':float(pressure_timescale/year/1e9),
                                'note':'Radius, total potential and fractional growth are frozen local diagnostics, not an evolved trajectory.'},
            'balanced_delivered_momentum_adiabatic_orbit':{'assumptions':'No radiation torque/force, constant external potential, slow loading, nearly circular orbit; delivered balance is stipulated, not derived from an isotropic moving-target bath.',
                                                         'mass_doubling_radius_ratio_flat_speed':.5,
                                                         'mass_doubling_radius_ratio_Kepler':.25},
            'next_requirements':['Net energy and momentum delivered after all emission channels',
                                 'Actual angular radiation distribution in the receiver frame',
                                 'Capture rates tied to receiver mass growth and finite capacity',
                                 'Angular momentum source and evolving supported distribution',
                                 'Joint source depletion, gravity and lensing']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'capture-momentum-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result['galactic_example'],indent=2))

if __name__=='__main__':main()

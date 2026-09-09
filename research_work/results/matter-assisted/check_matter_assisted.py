"""Recoil and a polarizability-mediated scalar companion: conditional EFT example."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'matter-assisted'

def dot(a,b): return a[0]*b[0]-a[1:]@b[1:]
def electric(k,e,u): return k*dot(e,u)-e*dot(k,u)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(2026090904)
    recoil=[]; ward=[]; angular=[]; disallowed=0
    u=np.array([1.,0,0,0])
    for mass in [10.,1000.,1e6]:
        for _ in range(100):
            energy=1.; final=rng.uniform(.05,.95)
            theta=rng.uniform(0,np.pi); phi=rng.uniform(0,2*np.pi)
            n=np.array([0.,0,1.]); np_=np.array([np.sin(theta),0,np.cos(theta)])
            qn=np.array([np.sin(phi),0,np.cos(phi)])
            delta=energy-final; transfer=energy*n-final*np_
            denominator=mass+delta-transfer@qn
            omega=(mass*delta-energy*final*(1-n@np_))/denominator
            if omega<=0:
                disallowed+=1
                continue
            target_p=transfer-omega*qn
            target_energy=np.sqrt(mass**2+target_p@target_p)
            # Rationalized expression avoids loss of small recoil energy.
            kinetic=(target_p@target_p)/(target_energy+mass)
            error=abs(delta-omega-kinetic)
            assert error<1e-12 and kinetic>=0
            recoil.append(dict(target_mass_over_initial_energy=mass,
                               companion_energy_fraction_of_photon_loss=float(omega/delta),
                               energy_residual=float(error)))
            k=np.r_[energy,energy*n]; kp=np.r_[final,final*np_]
            e=np.array([0.,1,0,0]); ep=np.r_[0.,np.cos(theta),0.,-np.sin(theta)]
            ward.extend([abs(dot(electric(k,k,u),electric(kp,ep,u))),
                         abs(dot(electric(k,e,u),electric(kp,kp,u)))])
    assert max(ward)<1e-12 and len(recoil)+disallowed==300
    # Sum explicit physical polarizations, independent of the scalar formula.
    for theta in np.linspace(0,np.pi,31):
        for phi in [0.,.3,1.]:
            n=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
            out_pol=[np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)]),
                     np.array([-np.sin(phi),np.cos(phi),0.])]
            value=sum((e@ep)**2 for e in [np.array([1.,0,0]),np.array([0.,1,0])] for ep in out_pol)/2
            expected=(1+np.cos(theta)**2)/2
            assert abs(value-expected)<1e-12
            angular.append(float(abs(value-expected)))
    # Direct integrated leading-heavy-target phase space for beta=1.
    # d sigma = E E'^3 (E-E') dE' dOmega_gamma dOmega_comp A / [8(2pi)^5].
    integrated=[]
    angular_integral=2*np.pi*quad(lambda mu:(1+mu*mu)/2,-1,1)[0]*4*np.pi
    for energy in [.5,1.,2.,4.]:
        sigma=angular_integral*quad(lambda ep:energy*ep**3*(energy-ep),0,energy)[0]/(8*(2*np.pi)**5)
        expected=energy**6/(480*np.pi**3)
        assert abs(sigma/expected-1)<1e-12
        integrated.append(dict(energy_in_reference_units=energy,sigma_for_unit_beta=float(sigma),
                               cross_section_ratio_to_reference=energy**6))
    norm=quad(lambda x:20*x**3*(1-x),0,1)[0]
    loss=quad(lambda x:20*x**3*(1-x)**2,0,1)[0]
    loss2=quad(lambda x:20*x**3*(1-x)**3,0,1)[0]
    assert abs(norm-1)<1e-12 and abs(loss-1/3)<1e-12 and abs(loss2-1/7)<1e-12
    # Heavy-target, dilute, independent forward channel. Calibrate only the
    # initial fractional drift at E_ref to alpha_ref, not an astronomical fit.
    # Dimensionless path t=alpha_ref R; jump rate is 3 E^6 in E_ref units.
    histories=[]; number=20000; length=.1
    for initial in [.5,1.,2.,4.]:
        energies=[]; counts=[]; max_error=0.
        for _ in range(number):
            energy=initial; companion=0.; position=0.; count=0
            while True:
                position+=rng.exponential(1/(3*energy**6))
                if position>length: break
                fraction=rng.beta(4,2)
                transferred=energy*(1-fraction)
                companion+=transferred; energy*=fraction; count+=1
            max_error=max(max_error,abs(energy+companion-initial))
            energies.append(energy); counts.append(count)
        energies=np.array(energies); counts=np.array(counts)
        no_event=float(np.mean(counts==0)); analytic=float(np.exp(-3*initial**6*length))
        se=np.sqrt(analytic*(1-analytic)/number)
        assert abs(no_event-analytic)<6*se+1/number and max_error<1e-12
        histories.append(dict(initial_energy_in_reference_units=initial,number_histories=number,
                              dimensionless_path=length,mean_final_energy=float(energies.mean()),
                              centroid_shift=initial/float(energies.mean())-1,
                              added_relative_energy_rms=float(energies.std()/energies.mean()),
                              unscattered_fraction=no_event,analytic_unscattered_fraction=analytic,
                              mean_event_count=float(counts.mean()),energy_ledger_max_error=max_error,
                              standard_error_mean_energy=float(energies.std()/np.sqrt(number))))
    straight=[]
    for mass in [10.,1000.,1e6]:
        delta=.1
        # photon undeflected, companion perpendicular: target carries recoil.
        omega=mass*delta/(mass+delta)
        kinetic=delta-omega
        straight.append(dict(target_mass_over_initial_energy=mass,photon_deflection=0,
                             companion_angle_radians=float(np.pi/2),companion_energy=omega,
                             recoil_energy=kinetic))
    result=dict(scope='One stipulated, leading-order polarizability EFT with scalar companions, not ordinary gravitons or a fitted complete theory',
                recoil_cases=len(recoil),attempted_recoil_cases=300,
                nonpositive_companion_energy_cases_disallowed=disallowed,
                max_recoil_energy_residual=max(r['energy_residual'] for r in recoil),
                max_ward_identity_residual=max(ward),max_polarization_sum_error=max(angular),
                leading_heavy_target_cross_sections=integrated,
                mean_event_loss=loss,second_event_loss_moment=loss2,
                photon_direction_mean_cosine=0.,straight_photon_recoil_examples=straight,
                illustrative_jump_histories=histories,
                checks=dict(recoil_conservation=True,gauge_invariance=True,polarization_sum=True,
                            heavy_target_phase_space=True,jump_moments=True,simulation_energy_and_survival=True),
                limitations=['Action coefficient and target abundance unmeasured; no absolute conversion rate predicted.',
                             'Heavy-target, long-wavelength polarizability limit only; higher operators and form factors omitted.',
                             'Scalar companion is a comparison identity, not a derived ordinary graviton.',
                             'Monte Carlo omits inverse reactions in a stipulated dilute escaping-companion limit.',
                             'No capture, stable deposit, gravitational response or event-time stretching derived.'])
    (OUT/'matter-assisted-results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],cross_section_ratios=integrated,histories=histories),indent=2))

if __name__=='__main__': main()

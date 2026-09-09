"""Finite-mode electromagnetic waves coupled to a dynamical matched index field."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'matched-wave'

def main():
    ks=np.array([.3,1.,3.,7.]);q0=np.array([.3,.2,.1,.08]);p0=np.array([.1,-.1,.2,-.1])
    actions=(p0*p0+ks*ks*q0*q0)/(2*ks);J=float(np.dot(ks,actions))
    runs=[]
    for inertia,ndot,end in [(10.,.02,20.),(1000.,.02,20.),(10.,5.,2.)]:
        # n, canonical field momentum, optical time, q_k, p_k.
        initial=np.r_[1.,inertia*ndot,0.,q0,p0]
        def rhs(t,y):
            n,field_p,eta=y[:3];q=y[3:7];p=y[7:11]
            radiation=np.sum(p*p+ks*ks*q*q)/(2*n)
            return np.r_[field_p/inertia,radiation/n,1/n,p/n,-ks*ks*q/n]
        sol=solve_ivp(rhs,(0,end),initial,method='DOP853',rtol=1e-11,atol=1e-13,
                      t_eval=np.linspace(0,end,1001))
        assert sol.success
        n,field_p,eta=sol.y[:3];q=sol.y[3:7];p=sol.y[7:11]
        current_actions=(p*p+ks[:,None]**2*q*q)/(2*ks[:,None])
        action_error=float(np.max(abs(current_actions-actions[:,None])))
        radiation=np.sum(p*p+ks[:,None]**2*q*q,axis=0)/(2*n)
        field=field_p*field_p/(2*inertia)
        energy_error=float(np.max(abs(field+radiation-field[0]-radiation[0])))
        phase=ks[:,None]*eta
        exact_q=q0[:,None]*np.cos(phase)+(p0/ks)[:,None]*np.sin(phase)
        exact_p=p0[:,None]*np.cos(phase)-(ks*q0)[:,None]*np.sin(phase)
        mode_error=float(max(np.max(abs(q-exact_q)),np.max(abs(p-exact_p))))
        assert action_error<1e-9 and energy_error<1e-8 and mode_error<1e-8
        reduced=solve_ivp(lambda t,y:[y[1],J/(inertia*y[0]**2)],(0,end),[1.,ndot],
                          method='DOP853',rtol=1e-12,atol=1e-14,t_eval=sol.t)
        assert reduced.success and np.max(abs(reduced.y[0]-n))<1e-8
        runs.append({'field_inertia':inertia,'initial_index_velocity':ndot,'end_time':end,
                     'initial_max_adiabaticity_parameter_abs_ndot_over_k':float(ndot/ks.min()),
                     'final_index':float(n[-1]),'all_mode_final_frequency_ratio':float(1/n[-1]),
                     'radiation_energy_lost':float(radiation[0]-radiation[-1]),
                     'field_energy_gained':float(field[-1]-field[0]),
                     'max_mode_action_error':action_error,'max_total_energy_error':energy_error,
                     'max_exact_optical_time_wave_error':mode_error,
                     'max_wave_vs_reduced_index_error':float(np.max(abs(reduced.y[0]-n)))})
    clocks=[]
    for ratio in [1.1,2.,10.]:
        # Minimal epsilon=mu=n, fixed charges/masses: nu_A~n^-2, length~n.
        emitted_atomic_frequency=1. # normalize emission epoch n_e=1
        received_line=emitted_atomic_frequency/ratio
        present_atomic_frequency=1/(ratio*ratio)
        observed_stretch=present_atomic_frequency/received_line
        assert np.isclose(observed_stretch,1/ratio)
        clocks.append({'n_observed_over_n_emitted':ratio,'propagation_frequency_ratio':1/ratio,
                       'same_species_spectroscopic_one_plus_z':observed_stretch,
                       'scope':'Leading Coulomb atom with fixed masses/charges; not a precision atomic theory or an adopted changing-rod model.'})
    result={'scope':'Positive classical finite-mode Hamiltonian, exact matched-wave evolution and leading charged-matter consistency diagnostic. No complete relativistic matter/gravity theory, cosmic history or graviton deposit.',
            'mode_wavenumbers':ks.tolist(),'wave_field_runs':runs,'minimal_atomic_response':clocks,
            'exact_invariant':'I_k=(p_k^2+k^2 q_k^2)/(2k); H=P_n^2/(2K)+sum k I_k/n, n>0.',
            'new_evidence':'The wave action gives exact homogeneous mode evolution and conserved receiving-field energy without an adiabatic approximation; it reproduces the prior ray model dynamically.',
            'known_limits_retained':['The minimal electromagnetic coupling changes electrostatics and atomic standards',
                                    'Fixed atoms cannot simply be stipulated independently of the same action',
                                    'Changing rods is not an adopted operationally nonexpanding solution',
                                    'The receiving n field is not an established graviton population or permanent well deposit',
                                    'Preferred-frame completion, spatial stability, dispersion and quantum vacuum require further work']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'matched-wave-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'runs':runs,'atomic_response':clocks},indent=2))

if __name__=='__main__':main()

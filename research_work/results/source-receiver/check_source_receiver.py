"""Join a specified converter spectrum to a receiver in a one-zone kinetic model."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad,solve_ivp

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'source-receiver'

def response(w,omega,damping_ratio):
    gamma=damping_ratio*omega
    denominator=(omega*omega-w*w)**2+gamma*gamma*w*w
    line=gamma*gamma*omega*omega/denominator # unity at omega, not necessarily maximum
    slope=-(4*w**4+2*(gamma*gamma-2*omega*omega)*w*w)/denominator
    return denominator,line,slope

def source(w,photon_energy,damping):
    cutoff=1. if np.isinf(photon_energy) else max(0.,1-w/photon_energy)**3
    return w*w*cutoff/response(w,1.,damping)[0]

def integrate(function,upper,receiver=1.):
    cuts=sorted(set([0.,*([x for x in [.5,1.,2.,receiver*.5,receiver,receiver*2] if x<upper]),upper]))
    return sum(quad(function,a,b,epsabs=1e-10,epsrel=2e-9,limit=200)[0] for a,b in zip(cuts[:-1],cuts[1:]))

def evaluate(photon_energy,damping,receiver,loading):
    total=integrate(lambda w:source(w,photon_energy,damping),photon_energy,receiver)
    def weight(w):
        line=response(w,receiver,damping)[1]
        return source(w,photon_energy,damping)*line/(1+loading*line)
    weighted=integrate(weight,photon_energy,receiver)
    drag=integrate(lambda w:weight(w)*(4+response(w,receiver,damping)[2])/3,photon_energy,receiver)/weighted
    captured=loading*weighted/total
    escape=integrate(lambda w:source(w,photon_energy,damping)/(1+loading*response(w,receiver,damping)[1]),photon_energy,receiver)/total
    assert abs(captured+escape-1)<1e-8
    assert -1e-12<=captured<=1+1e-12
    return {'photon_energy_over_emitter_Omega':photon_energy,'emitter_damping_ratio':damping,
            'receiver_Omega_over_emitter_Omega':receiver,'receiver_damping_ratio':damping,
            'capture_rate_at_receiver_resonance_times_escape_time':loading,
            'steady_received_power_fraction':captured,'steady_escape_power_fraction':escape,
            'normalized_absorption_drag_coefficient':drag,
            'receiving_to_incident_power_ledger_error':abs(captured+escape-1)}

def main():
    matched=[]
    for damping in [.03,.1,.5]:
        for energy in [10.,100.,np.inf]:
            row=evaluate(energy,damping,1.,0.)
            # Matched emitter/receiver response allows an independent exact identity.
            def weighted(w):return source(w,energy,damping)/response(w,1.,damping)[0]
            normal=integrate(weighted,energy)
            correction=0 if np.isinf(energy) else .5*integrate(lambda w:weighted(w)*w/(energy-w),energy)/normal
            expected=5/6+correction
            assert abs(row['normalized_absorption_drag_coefficient']-expected)<1e-8
            row['matched_analytic_coefficient']=expected
            row['identity_error']=abs(row['normalized_absorption_drag_coefficient']-expected)
            if np.isinf(energy):row['photon_energy_over_emitter_Omega']='asymptotic reference, not finite source'
            matched.append(row)
    scans=[evaluate(100.,.1,receiver,loading) for receiver in [.5,.8,1.,1.2,2.]
           for loading in [.01,1.,100.,10000.]]
    # Frequency-resolved energy ledger tested dynamically rather than only at steady state.
    kinetic=[]
    for w in [.2,1.,3.,30.]:
        production=source(w,100.,.1);capture=10*response(w,1.,.1)[1];escape=1.
        rate=capture+escape;end=10/rate
        def rhs(t,y):
            reservoir,received,escaped=y
            return [production-rate*reservoir,capture*reservoir,escape*reservoir]
        sol=solve_ivp(rhs,(0,end),[0.,0.,0.],rtol=1e-11,atol=1e-13)
        exact=production/rate*(-np.expm1(-rate*sol.t))
        error=float(np.max(abs(sol.y[0]-exact)))
        ledger=float(np.max(abs(sol.y.sum(axis=0)-production*sol.t)))
        assert sol.success and error<1e-8 and ledger<1e-8
        kinetic.append({'companion_energy_over_emitter_Omega':w,'capture_rate_times_escape_time':capture,
                        'max_reservoir_error':error,'max_energy_ledger_error':ledger})
    result={'scope':'Conditional dilute/incoherent source-receiver spectral closure and one-zone energy kinetics. The converter retains its known color/angle failures. No complete microscopic theory, geometry, permanent deposit, stellar supply normalization or observational validation.',
            'source_spectrum':'S(w) proportional to w^2 (1-w/E)^3 / [(1-w^2)^2+d_em^2*w^2], 0<w<E; derived from the earlier point-oscillator conversion differential rate.',
            'kinetic_equation':'du_w/dt=S_w-(Gamma_abs(w)+1/t_escape)u_w; receiving power Gamma_abs(w)u_w and escape power u_w/t_escape are explicit.',
            'optically_weak_matched_receiver':matched,'detuning_and_loading_scan':scans,
            'time_dependent_energy_checks':kinetic,
            'force_branch_note':'C is normalized to receiving-channel absorption power. In the constant-branch scalar oscillator, reradiated companion scattering has the same spectral shape and adds force; permanent storage efficiency and receiver identity remain separate.',
            'unresolved':['Absolute converter rate, stellar photon spectrum and finite source history',
                          'Spatial/angular transport and field-mediated emitter-receiver correlations',
                          'Real receiver, capacity, heating and inverse transitions',
                          'Full momentum balance, orbital support and evolving gravity',
                          'Achromatic redshift, sharp images, timing and joint lensing']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'source-receiver-results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'matched_example':matched[4],'loading_examples':[x for x in scans if x['receiver_Omega_over_emitter_Omega']==1.]},indent=2))

if __name__=='__main__':main()

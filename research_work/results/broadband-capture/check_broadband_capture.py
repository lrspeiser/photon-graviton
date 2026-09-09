"""Broadband drag criterion and an explicitly normalized scalar oscillator receiver."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'broadband-capture'

def shape(energy,damping):
    denominator=(1-energy*energy)**2+damping*damping*energy*energy
    cross=1/denominator
    slope=-(4*energy**4+2*(damping*damping-2)*energy*energy)/denominator
    return cross,slope

def thermal_integral(function):
    return sum(quad(function,a,b,epsabs=1e-9,epsrel=2e-10,limit=200)[0]
               for a,b in [(0,.5),(.5,1),(1,2),(2,np.inf)])

def main():
    ledger=[]
    # g^2=4*pi*gamma_rad; all frequencies and rates in Omega=1 units.
    for gr,gs in [(.01,.09),(.1,.9),(.1,0.)]:
        g2=4*np.pi*gr
        for omega in [.1,1.,10.]:
            den=(1-omega*omega)**2+(gr+gs)**2*omega*omega
            amplitude=np.sqrt(g2)/complex(1-omega*omega,-(gr+gs)*omega)
            flux=omega*omega/2 # unit incident scalar amplitude
            radiation=gr*omega*omega*abs(amplitude)**2/2
            receiving=gs*omega*omega*abs(amplitude)**2/2
            input_work=np.sqrt(g2)*omega*amplitude.imag/2
            assert abs(input_work-radiation-receiving)<1e-10
            assert np.isclose(receiving/flux,g2*gs/den,rtol=1e-12,atol=1e-14)
            forward_amplitude=g2/(4*np.pi*complex(1-omega*omega,-(gr+gs)*omega))
            extinction=4*np.pi*forward_amplitude.imag/omega
            assert abs(extinction-(radiation+receiving)/flux)<1e-10
            ledger.append({'gamma_rad':gr,'gamma_receiving':gs,'omega_over_Omega':omega,
                           'radiated_power':radiation,'receiver_power':receiving,
                           'input_work':input_work,'absorption_cross_section':receiving/flux})
    thermal=[]
    for damping in [.1,1.]:
        for temperature in [.1,1.,10.]:
            def spectrum(e):
                x=e/temperature
                return 0. if x>700 else e**3/np.expm1(x)
            def power(e):return spectrum(e)*shape(e,damping)[0]
            def direct(e):
                cross,slope=shape(e,damping)
                return spectrum(e)*cross*(4+slope)/3
            def positive(e):
                x=e/temperature
                return power(e)*x/(-np.expm1(-x))/3
            weight=thermal_integral(power)
            a=thermal_integral(direct)/weight;b=thermal_integral(positive)/weight
            assert abs(a-b)<1e-7 and b>1/3
            thermal.append({'gamma_over_Omega':damping,'temperature_over_Omega':temperature,
                            'drag_coefficient_direct':a,'drag_coefficient_positive_form':b,
                            'identity_error':abs(a-b)})
    monochromatic=[]
    damping=.1;zero=np.sqrt(2/(2-damping*damping))
    for energy in [zero,2.,10.,100.]:
        cross,slope=shape(energy,damping)
        monochromatic.append({'energy_over_Omega':energy,'cross_section_shape':cross,
                              'small_speed_drag_coefficient':(4+slope)/3})
    assert abs(monochromatic[0]['small_speed_drag_coefficient'])<1e-10
    bands=[]
    for center in [zero,10.]:
        for width in [.001,.05,.2,.5]:
            # Gaussian in log energy. Integrate in standard-normal variable;
            # split at the oscillator resonance to resolve a narrow line.
            cuts=sorted(set([-10.,10.,0.,*([np.log(1/center)/width] if abs(np.log(1/center)/width)<10 else [])]))
            def term(z,moment):
                e=center*np.exp(width*z);cross,slope=shape(e,damping)
                return np.exp(-z*z/2)*cross*(1 if moment==0 else (4+slope)/3)
            sums=[sum(quad(lambda z:term(z,m),a,b,epsabs=1e-10,epsrel=1e-9,limit=200)[0]
                      for a,b in zip(cuts[:-1],cuts[1:])) for m in [0,1]]
            bands.append({'center_energy_over_Omega':center,'log_energy_standard_deviation':width,
                          'band_averaged_drag_coefficient':sums[1]/sums[0]})
    # An unregularized epsilon^-4 cross section in a thermal bath has divergent
    # power at low frequency. Show the missing cutoff cannot be ignored.
    divergent=[]
    for cutoff in [.01,.001,.0001]:
        value=quad(lambda e:1/(e*np.expm1(e)),cutoff,1,epsabs=1e-8)[0]
        divergent.append({'lower_energy_cutoff_over_temperature':cutoff,
                          'partial_absorbed_power_in_arbitrary_units':value})
    assert all(b['partial_absorbed_power_in_arbitrary_units']>5*a['partial_absorbed_power_in_arbitrary_units']
               for a,b in zip(divergent[:-1],divergent[1:]))
    result={'scope':'Passive scalar oscillator and low-speed broadband radiation-force diagnostics. No ordinary-graviton receiver, permanent store, coupled galaxy or astronomical validation.',
            'driven_oscillator_ledger':ledger,'thermal_drag':thermal,
            'monochromatic_oscillator_response':monochromatic,'finite_band_response':bands,
            'unregularized_inverse_fourth_power_thermal_problem':divergent,
            'general_criterion':'C=[integral u(epsilon) sigma(epsilon) (4+dlnsigma/dlnepsilon) depsilon]/[3 integral u sigma depsilon]. For vanishing endpoint term this equals integral sigma(3u-epsilon du/depsilon)/(3 integral u sigma).',
            'thermal_scope':'For isotropic Planck spectral energy density, passive nonnegative sigma, finite integrals and vanishing endpoint term, C>1/3. Does not assume companions are thermal or invoke a Big Bang.',
            'unresolved':['Physical receiver, continuum and ultraviolet cutoff',
                          'Actual companion spectrum from the same conversion mechanism',
                          'Permanent storage, capacity and receiver heating',
                          'Scattering and all reverse/emitted channels',
                          'Coupled transport, orbit evolution, gravity and lensing']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'broadband-capture-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'thermal_drag':thermal,'finite_band_response':bands},indent=2))

if __name__=='__main__':main()

"""Causal oscillator response, energy ledger and the inherited spectral kernel."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad,solve_ivp

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'oscillator-response'

def moments(r,d):
    # u=companion frequency / oscillator frequency; r=photon E / that frequency.
    # Log-frequency integration resolves both a narrow resonance and a long tail.
    ans=[]
    for power in range(3):
        def ordinary(u): return (1-u/r)**3*u**(power+1)/((1-u*u)**2+d*d*u*u)
        a=quad(ordinary,0,.5,epsabs=1e-11,epsrel=1e-10)[0]
        def logarithmic(v):
            u=np.exp(v)
            return ordinary(u)*u
        boundaries=sorted(set([np.log(.5),-d,0.,d,np.log(r)]))
        value=a+sum(quad(logarithmic,lo,hi,epsabs=1e-10,epsrel=1e-10,limit=200)[0]
                    for lo,hi in zip(boundaries[:-1],boundaries[1:]))
        ans.append(value)
    return ans

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    response=[]
    for gamma in [.02,.1,.5]:
        def rhs(t,y):
            q,v,radiated=y
            return [v,-q-gamma*v,gamma*v*v]
        sol=solve_ivp(rhs,(0,12/gamma),[0.,1.,0.],rtol=2e-10,atol=1e-12)
        energy=(sol.y[0]**2+sol.y[1]**2)/2
        residual=float(np.max(abs(energy+sol.y[2]-.5)))
        wd=np.sqrt(1-gamma**2/4)
        exact=np.exp(-gamma*sol.t/2)*np.sin(wd*sol.t)/wd
        error=float(np.max(abs(sol.y[0]-exact)))
        assert sol.success and residual<1e-8 and error<1e-8
        # Static susceptibility Kramers-Kronig sum rule for the ideal response.
        integral=sum(quad(lambda w:gamma/((1-w*w)**2+gamma*gamma*w*w),lo,hi,
                          epsabs=1e-11,epsrel=1e-10)[0] for lo,hi in [(0,1),(1,np.inf)])
        static=2*integral/np.pi
        assert abs(static-1)<1e-9
        poles=np.roots([1.,gamma,1.]); assert np.max(poles.real)<0
        response.append(dict(gamma_over_omega=gamma,energy_ledger_max_error=residual,
                             analytic_impulse_max_error=error,static_sum_rule=static,
                             final_stored_energy=float(energy[-1]),radiated_energy=float(sol.y[2,-1])))
    # A forced oscillator makes the source work explicit, rather than assigning
    # a dissipative filter an unaccounted energy sink.
    def drive(t): return np.sin(t)*np.exp(-(t-10)**2/8)
    def forced(t,y):
        q,v,rad,work=y
        return [v,drive(t)-q-.1*v,.1*v*v,drive(t)*v]
    sol=solve_ivp(forced,(0,200),[0.,0.,0.,0.],rtol=2e-10,atol=1e-12)
    residual=float(np.max(abs((sol.y[0]**2+sol.y[1]**2)/2+sol.y[2]-sol.y[3])))
    assert sol.success and residual<1e-8
    spectra=[]; colors=[]
    for d in [.01,.1]:
        for r in [10.,100.,1000.,1e6]:
            m=moments(r,d)
            spectra.append(dict(gamma_over_omega=d,photon_energy_over_resonance=r,
                                mean_fractional_transfer=m[1]/m[0]/r,
                                second_fractional_transfer=m[2]/m[0]/r**2,
                                mean_companion_energy_over_resonance=m[1]/m[0],
                                companion_energy_rms_over_resonance=float(np.sqrt(m[2]/m[0]-(m[1]/m[0])**2))))
        reference=moments(100,d)[1]
        for ratio in [.5,1.,2.,4.]:
            r=100*ratio
            value=ratio**3*moments(r,d)[1]/reference
            if ratio>=1: assert value>=ratio**3-1e-10
            colors.append(dict(gamma_over_omega=d,photon_energy_ratio=ratio,
                               initial_fractional_loss_ratio=value))
    # The point-target slope bound applies beyond the Lorentzian. Independently
    # integrate positive trial frequency weights: alpha(E)=int (E-w)^3 W(w) dw.
    bound=[]
    for name,weight in [('exponential',lambda w:np.exp(-w)),
                        ('two_bands',lambda w:np.exp(-((w-.3)/.1)**2)+.3*np.exp(-((w-2)/.2)**2)),
                        ('broad',lambda w:1/(1+w*w))]:
        for e in [1.,3.,10.]:
            a=quad(lambda w:(e-w)**3*weight(w),0,e,epsabs=1e-11)[0]
            derivative=3*quad(lambda w:(e-w)**2*weight(w),0,e,epsabs=1e-11)[0]
            slope=e*derivative/a
            assert slope>=3-1e-10
            bound.append(dict(weight=name,photon_energy=e,logarithmic_slope=slope))
    result=dict(scope='Point-target oscillator-mediated scalar comparison; no combined collective array or ordinary graviton law',
                response_checks=response,forced_energy_ledger_max_error=residual,
                forced_source_work=float(sol.y[3,-1]),forced_radiated_energy=float(sol.y[2,-1]),
                transfer_spectra=spectra,color_dependence=colors,positive_weight_slope_checks=bound,
                checks=dict(impulse_response=True,stable_response_poles=True,static_dispersion_sum_rule=True,
                            stored_radiated_and_input_work_balance=True,positive_weight_color_slope_bound=True),
                limitations=['Oscillator frequency and couplings are unspecified physical parameters, not fitted discoveries.',
                             'Radiative scalar damping approximation assumes negligible competing channels and renormalized positive oscillator frequency.',
                             'Two-photon decay from the light coupling and incoming-companion fluctuations cannot be discarded in a complete model.',
                             'No jointly derived extended-target self-energy, stable deposit, event dilation or gravity/lensing response.'])
    (OUT/'oscillator-response-results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],spectra=spectra,color_dependence=colors),indent=2))

if __name__=='__main__': main()

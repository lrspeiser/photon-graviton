"""Inverse source requirements for a stationary Planck bath under energy drift.

Dimensionless x=h_P nu/(k_B T), time in units of the conversion time.
This derives the required source, not an astronomical source model.
"""
from pathlib import Path
import json, os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'background-replenishment'

def number(x):
    return x*x/np.expm1(x) if x < 700 else 0.0

def slope_factor(x):
    return x/(-np.expm1(-x))

def source(x, a):
    return number(x)*(a+slope_factor(x)-3)

def main():
    n=quad(number,0,100,epsabs=1e-11)[0]
    u=quad(lambda x:x*number(x),0,100,epsabs=1e-11)[0]
    crossing=brentq(lambda x:slope_factor(x)-3,1,5,xtol=1e-14)
    rows=[]
    for a in [0.,1.,2.,3.,10.]:
        jn=quad(lambda x:source(x,a),0,100,epsabs=1e-10)[0]
        ju=quad(lambda x:x*source(x,a),0,100,epsabs=1e-10)[0]
        assert np.isclose(jn,a*n,atol=1e-10,rtol=1e-10)
        assert np.isclose(ju,(a+1)*u,atol=1e-10,rtol=1e-10)
        # Independent differential check of the inverse source equation.
        residuals=[]
        for x in np.geomspace(.001,30,120):
            dx=x*1e-4
            f=lambda y:y*number(y)
            derivative=(-f(x+2*dx)+8*f(x+dx)-8*f(x-dx)+f(x-2*dx))/(12*dx)
            residuals.append(abs(derivative+source(x,a)-a*number(x))/number(x))
        assert max(residuals)<1e-7
        rows.append({'removal_rate_over_conversion_rate':a,
                     'low_frequency_source_over_conversion_rate_times_number_limit':a-2,
                     'source_nonnegative_at_all_positive_frequencies':a>=2,
                     'injected_number_over_bath_number_per_conversion_time':jn/n,
                     'injected_energy_over_bath_energy_per_conversion_time':ju/u,
                     'max_scaled_differential_residual':max(residuals)})
    # Integrate the transport solution along characteristics, independently of
    # the local derivative test. The prescribed positive source must retain N.
    characteristic_errors=[]
    for a in [2.,3.,10.]:
        for t in [.1,.5,1.]:
            for x in np.geomspace(.001,10,20):
                evolved=np.exp((1-a)*t)*number(x*np.exp(t))+quad(
                    lambda s:np.exp((1-a)*s)*source(x*np.exp(s),a),0,t,
                    epsabs=1e-12,epsrel=1e-11)[0]
                characteristic_errors.append(abs(evolved/number(x)-1))
    assert max(characteristic_errors)<1e-9
    result={'scope':'Exact inverse source diagnostic for a fixed-volume, constant-speed, stationary Planck bath with frequency-independent conversion and removal rates. No source population, capture mechanism, FIRAS fit or nonexpanding history is established.',
            'conversion_rate_symbol':'H_c; not the expansion rate or Planck constant',
            'minimum_gray_removal_rate_over_conversion_rate':2,
            'zero_removal_negative_source_below_dimensionless_frequency':crossing,
            'dimensionless_planck_number_integral':n,
            'dimensionless_planck_energy_integral':u,
            'mean_bath_photon_energy_over_kT':u/n,
            'minimum_gray_case_mean_injected_photon_energy_over_kT':1.5*u/n,
            'gray_removal_scenarios':rows,
            'max_characteristic_stationarity_relative_error':max(characteristic_errors),
            'energy_ledger':'P_source=(H_c+kappa) U; P_conversion=H_c U; P_removed=kappa U. Removed energy must be received or exported. A stationary photon bath does not imply stationary companions, deposits, fuel or heat.',
            'not_adopted':['engineered source spectrum','gray photon removal','permanent storage','thermal companion population']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'background-replenishment-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()

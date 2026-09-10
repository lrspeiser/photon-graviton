"""Conditional path-resolved transfer/capture benchmark, not a time-field derivation."""
import json
import math
from pathlib import Path
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
records=[]
# Photon speed 1, conversion x=0..1; companions move independently at speed 0.5.
# Capture only x=2..3, with dimensionless integrated optical depth K.
for A in (0.,.1,math.log(2)):
    for K in (0.,1.,5.):
        for shape in ('uniform','increasing'):
            alpha=(lambda x:A) if shape=='uniform' else (lambda x:2*A*x)
            depth=(lambda x:A*x) if shape=='uniform' else (lambda x:A*x*x)
            produced=quad(lambda x:alpha(x)*math.exp(-depth(x)),0,1,epsabs=1e-13)[0]
            photons=math.exp(-A)
            deposited=quad(lambda x:alpha(x)*math.exp(-depth(x))*(-math.expm1(-K)),0,1,epsabs=1e-13)[0]
            escaped=quad(lambda x:alpha(x)*math.exp(-depth(x))*math.exp(-K),0,1,epsabs=1e-13)[0]
            assert abs(produced-(-math.expm1(-A)))<1e-12
            assert abs(photons+deposited+escaped-1)<1e-12
            assert abs(deposited-(-math.expm1(-A))*(-math.expm1(-K)))<1e-12
            # Every signal traverses fixed length 4 at fixed speed 1.
            launch=[0.,.2,.4];arrival=[t+4 for t in launch]
            timing=(arrival[-1]-arrival[0])/(launch[-1]-launch[0])
            assert abs(timing-1)<1e-12
            records.append({'integrated_conversion_A':A,'integrated_capture_K':K,'conversion_shape':shape,
                'photon_energy_fraction_at_detector':photons,'deposited_energy_fraction':deposited,
                'escaped_companion_energy_fraction':escaped,'wavelength_stretch':math.exp(A),
                'whole_event_stretch':timing,'arrival_times':arrival,
                'total_accounted_energy_fraction':photons+deposited+escaped})
out={'status':'conditional spatial bookkeeping and timing counterexample; no microscopic coupling, time-field source or gravity response',
     'formula_provenance':'standard continuity/characteristic solutions applied to stipulated transfer/capture laws; not unique formulas',
     'cases':len(records),'companion_speed':.5,'photon_speed':1.,
     'assumptions':['fixed nonexpanding spatial coordinates and fixed speeds',
                    'stationary frequency-independent fractional energy-transfer rate',
                    'photon number retained; no scattering, delay or reverse transfer',
                    'companions keep energy during free flight and lose it only into deposits',
                    'conversion and capture regions disjoint; all companions traverse capture region',
                    'permanent deposited reservoir; recoil and a gravitational law unmodeled'],
     'records':records}
(HERE/'spatial-transfer-results.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps({'cases':len(records),'max_energy_error':max(abs(r['total_accounted_energy_fraction']-1) for r in records),
                  'half_energy_case':[r for r in records if r['integrated_conversion_A']==math.log(2) and r['integrated_capture_K']==1 and r['conversion_shape']=='uniform']},indent=2))

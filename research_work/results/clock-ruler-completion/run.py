"""Observable standards in a conditional homogeneous universal metric family."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
A=.1

def main():
    sources=[Path(__file__),HERE.parent/'directional-wave-reservoir/run.py',
             HERE.parent/'atomic-line-response/report.md',HERE.parent/'matter-clock-closure/report.md']
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    rows=[]
    for exponent in [0.,.5,1.,1.5,2.,2.5]:
        for distance in [1.,3.,5.]:
            S=np.exp(A*distance)
            for birth in [0.,.5,1.]:
                def arrival(t):
                    return brentq(lambda end:quad(lambda q:1/(1+A*q),t,end,epsabs=1e-12)[0]-distance,t,t+20,xtol=1e-13)
                first,last=arrival(birth),arrival(birth+.2)
                ne,no=1+A*birth,1+A*first
                lapse=lambda t:(1+A*t)**(exponent-2)
                emit_interval=quad(lapse,birth,birth+.2,epsabs=1e-13)[0]
                receive_interval=quad(lapse,first,last,epsabs=1e-13)[0]
                spectral=S*(no/ne)**(exponent-2)
                exact=S**(exponent-1)
                assert abs(spectral-exact)<1e-12 and abs(receive_interval/emit_interval-exact)<1e-11
                # Atomic coordinate length scales as n^(1-exponent); clocks as n^(exponent-2).
                n=no;clock=n**(exponent-2);rod=n**(1-exponent);speed=1/n
                measured=speed/(rod*clock)
                assert abs(measured-1)<1e-12
                # Hamiltonian small-momentum kinetic mass from the metric completion.
                spatial=n**(exponent-1);N=clock;p=.03
                H_metric=N*np.sqrt(1+p*p/spatial**2)
                H_disp=np.sqrt((p/n)**2+n**(2*(exponent-2)))
                assert abs(H_metric-H_disp)<1e-12
                rows.append(dict(inertial_mass_exponent=exponent,distance=distance,birth_time=birth,
                    reference_frequency_stretch=S,measured_one_plus_redshift=float(spectral),
                    measured_event_interval_ratio=float(receive_interval/emit_interval),
                    measured_light_speed=measured,spatial_scale_ratio=float((no/ne)**(exponent-1)),
                    lapse_ratio=float((no/ne)**(exponent-2))))
    examples=[]
    for r in [0.,1.,2.]:
        S=1.4
        examples.append(dict(exponent=r,index_ratio=S,lapse_ratio=S**(r-2),
            coordinate_atomic_rod_ratio=S**(1-r),physical_spatial_scale_ratio=S**(r-1),
            measured_one_plus_redshift=S**(r-1),measured_c_ratio=1.))
    result=dict(scope='Conditional homogeneous clock/ruler completion, not a general theorem against nonexpanding mechanisms.',
        assumed_family='epsilon=mu=n; fixed dimensionless atomic coupling; common matter/photon limiting speed 1/n; coordinate inertial mass n^r.',
        derived_metric='ds^2=-n^(2r-4) dt^2+n^(2r-2) dx^2',
        rows=rows,examples=examples,source_hashes=hashes,
        matter_coupling_derived_from_wave_hamiltonian=False,einstein_equations_solved=False,
        nonuniform_endpoints_tested=False,holdouts_opened=False)
    for p in sources:assert hashlib.sha256(p.read_bytes()).hexdigest()==hashes[str(p.relative_to(ROOT))]
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(examples,indent=2))

if __name__=='__main__':main()

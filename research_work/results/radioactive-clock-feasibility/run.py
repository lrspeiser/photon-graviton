"""Conditional radioactive-tail clock with an illustrative deposition law."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
TAU=111.3 # illustrative Co-56-like e-folding benchmark, not a new measurement
rows=[]
for S in [1.,1.5,2.2]:
    for b in [0.,1.]:
        A=S**b
        for t0 in [35.,70.]:
            for t in [100.,200.,400.]:
                fpos=.035;x=(t0/t)**2;dep=fpos+(1-fpos)*(-np.expm1(-x))
                derivative=-(1-fpos)*np.exp(-x)*2*x/t
                observed_rate=(1/TAU-derivative/dep)/A
                def flux(observer_time):
                    source_time=observer_time/A
                    return np.exp(-source_time/TAU)*(fpos+(1-fpos)*(-np.expm1(-(t0/source_time)**2)))/S**(1+b)
                eps=.001*A
                numerical=-(np.log(flux(A*t+eps))-np.log(flux(A*t-eps)))/(2*eps)
                assert abs(numerical/observed_rate-1)<1e-8
                naive=1/(TAU*observed_rate);assert naive<=A+1e-12
                rows.append({'S':S,'b':b,'true_event_stretch':A,'source_age_days':t,'t0_days':t0,'deposition_fraction':float(dep),'observed_log_decline_per_day':float(observed_rate),'naive_stretch_if_escape_ignored':float(naive),'naive_to_true_stretch':float(naive/A),'finite_difference_relative_error':float(abs(numerical/observed_rate-1))})
(HERE/'results.json').write_text(json.dumps({'scope':'Conditional source-physics stress test; no observational fit or proven independent clock.','tau_days_benchmark':TAU,'illustrative_retained_particle_energy_fraction':.035,'rows':rows,'maximum_derivative_error':max(v['finite_difference_relative_error'] for v in rows)},indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'checks':len(rows),'minimum_naive_to_true_stretch':min(v['naive_to_true_stretch'] for v in rows),'maximum_naive_to_true_stretch':max(v['naive_to_true_stretch'] for v in rows),'max_derivative_error':max(v['finite_difference_relative_error'] for v in rows)}))

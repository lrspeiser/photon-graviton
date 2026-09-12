"""Check nested width grids with identical shape quadrature."""
from pathlib import Path
import json, sys
import numpy as np
from scipy.special import logsumexp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'timing-population'))
from estimator import fit_population,log_width_weights
protocol=json.loads((HERE/'protocol.json').read_text())
gate=json.loads((HERE/'width-check-protocol.json').read_text())
with np.load(HERE/'refined-event-likelihoods.npz') as a:
    x=a['log_width'];logs=a['event_log_likelihoods'];z=a['redshift']
    coarse=fit_population(x[::2],logs[:,::2],z,protocol)
    fine=fit_population(x,logs,z,protocol)
    cc=logs[:,::2]-logs[:,::2].max(axis=1)[:,None]
    fc=logs-logs.max(axis=1)[:,None]
    interp=np.array([np.interp(x,x[::2],row) for row in cc])
    post=fc+log_width_weights(x,z,[fine['a'],fine['b'],np.log(fine['sigma'])])
    post=np.exp(post-logsumexp(post,axis=1)[:,None])
    change=float(np.mean(np.sum(post*abs(interp-fc),axis=1)))
    bd=abs(coarse['b']-fine['b'])
    passed=bd<=gate['maximum_b_change'] and change<=gate['maximum_mean_event_log_change'] and all(r['optimizer_success'] and r['interior_solution'] and r['maximum_population_mass_outside_width_bounds']<=protocol['numerical_gate']['maximum_continuous_population_probability_outside_width_bounds'] for r in [coarse,fine])
    result={'coarse':coarse,'fine':fine,'absolute_b_change':bd,'mean_event_log_change':change,'nested_grid_check_pass':bool(passed),'scope':gate['scope'],'limitation':gate['limitation']}
    (HERE/'width-grid-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'absolute_b_change':bd,'mean_event_log_change':change,'nested_grid_check_pass':bool(passed)}))

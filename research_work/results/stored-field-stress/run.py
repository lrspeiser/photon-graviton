"""Stress consistency checks; no observational fit or held-out score."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import quad

here = Path(__file__).resolve().parent
source = here.parent/'deposit-boundary-admissibility/configurations.json'
fit = here.parent/'joint-galaxy-audit/results.json'
rows = json.loads(source.read_text())
p = json.loads(fit.read_text())['sparc']['parameters']['p']
assert len(rows) == 99 and all(r['role']=='training' for r in rows)
errors = []
for row in rows:
    c = row['c_dimensionless']
    xt = row['cutoff_over_Re']*1.8153
    for x in np.geomspace(xt*.001,xt*.999,12):
        g = c/(1+x)**(2*p)
        lap = g*(2/x-2*p/(1+x))
        assert lap > 0
        step = x*1e-5
        def mass(t): return t*t*c/(1+t)**(2*p)
        numerical = (mass(x+step)-mass(x-step))/(2*step*x*x)
        errors.append(abs(numerical/lap-1))
assert max(errors)<1e-7
# Q=sin(t)^2/2 and V=cos(t)^2/2 for m=F=1.
mean_u = quad(lambda t:(np.sin(t)**2+np.cos(t)**2)/2,0,2*np.pi)[0]/(2*np.pi)
mean_p = quad(lambda t:(np.sin(t)**2-np.cos(t)**2)/2,0,2*np.pi)[0]/(2*np.pi)
assert abs(mean_u-.5)<1e-12 and abs(mean_p)<1e-12
result = {'classification':'Conditional canonical scalar stress audit; no galaxy validation',
          'training_configurations':len(rows),'positive_source_probes':len(errors),
          'max_laplacian_relative_error':max(errors),
          'static_nonnegative_potential_scalar_matches_required_source':False,
          'harmonic_local_example':{'mean_energy':mean_u,'mean_pressure':mean_p,
                                    'finite_spatial_equilibrium_established':False},
          'input_sha256':{str(f.relative_to(here.parent)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [source,fit]}}
(here/'results.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))

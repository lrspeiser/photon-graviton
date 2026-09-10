"""Dimensionless exterior checks, not fitted cutoff radii."""
from pathlib import Path
import json
import numpy as np
H=Path(__file__).resolve().parent
p=json.loads((H/'source-density-results.json').read_text())['parameters']['p']
q=2*p;x=np.array([1.,2.,10.,100.])
# Acceleration unit g_t at a hypothetical outer transition r_t; mass unit
# g_t*r_t^2/G. Interior scaling is exterior point-baryon asymptotics only.
rows=[]
for s in x:
    gpower=s**(-q);gcut=gpower*np.exp(-(s-1));gfinite=s**(-2)
    rows.append(dict(r_over_rt=s,power_g_over_gt=gpower,power_M_over_Mt=s**(2-q),
        naive_exponential_g_over_gt=gcut,naive_exponential_density_sign_factor=2-q-s,
        source_truncated_g_over_gt=gfinite,source_truncated_M_over_Mt=1.))
# Known spherical source truncation retains enclosed mass; force remains
# continuous, exterior density zero, and the exterior potential integral finite.
assert all(r['source_truncated_g_over_gt']*r['r_over_rt']**2==1 for r in rows)
out={'p':p,'hypothetical_rt_value_assigned':False,'rows':rows,
     'naive_acceleration_exponential_becomes_negative_density_at_r_over_rt':2-q,
     'interpretation':'A positive-source outer continuation is mathematically available. Its radius, smoothness, source support and capture cause are not predicted or fitted.'}
(H/'exterior-check.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))

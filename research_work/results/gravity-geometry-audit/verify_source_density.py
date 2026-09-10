"""Quantify source-reconstruction error; preserve negative numerical artifacts."""
from pathlib import Path
import json
import numpy as np
H=Path(__file__).resolve().parent
rows=json.loads((H/'source-density-rows.json').read_text())
def sample(order,h):return [r for r in rows if r['bar_order']==order and r['h_kpc']==h]
a=sample(64,.001);b=sample(64,.002);c=sample(40,.001)
assert len(a)==len(b)==len(c)==63
assert [(r['R_kpc'],r['z_kpc']) for r in a]==[(r['R_kpc'],r['z_kpc']) for r in b]==[(r['R_kpc'],r['z_kpc']) for r in c]
out={}
for label,key in [('spherical','spherical_density_Msun_pc3'),('shape','shape_density_Msun_pc3')]:
    vals=np.array([r[key] for r in a]);step=np.array([abs(x[key]-y[key]) for x,y in zip(a,b)]);order=np.array([abs(x[key]-y[key]) for x,y in zip(a,c)])
    out[label]={'negative_probes_in_any_run':sum(r[key]<0 for r in rows),
                'max_step_change_Msun_pc3':float(step.max()),'max_bar_order_change_Msun_pc3':float(order.max()),
                'minimum_density_minus_step_and_order_changes_Msun_pc3':float(np.min(vals-step-order))}
out['chain_rule_vs_acceleration_divergence_max_difference_Msun_pc3']=max(abs(r['shape_density_Msun_pc3']-r['shape_divergence_density_Msun_pc3']) for r in a)
out['declared_baryon_minimum_Msun_pc3']=min(r['declared_baryon_density_Msun_pc3'] for r in a)
out['maximum_baryonic_density_reconstruction_error_Msun_pc3']=max(abs(r['baryon_density_Msun_pc3']-r['declared_baryon_density_Msun_pc3']) for r in a)
out['maximum_shape_change_when_using_declared_density_Msun_pc3']=max(abs(r['shape_density_Msun_pc3']-r['shape_with_declared_source_Msun_pc3']) for r in a)
out['shape_source_substitution_minimum_Msun_pc3']=min(r['shape_with_declared_source_Msun_pc3'] for r in a)
out['negative_reconstructed_baryon_examples']=[r for r in a if r['baryon_density_Msun_pc3']<0]
out['limits']='Finite probe and resolution diagnostics, not a global positivity theorem or rigorous error bound. Declared-source substitution is a sensitivity check, not the exact Laplacian of the cached potential.'
assert out['chain_rule_vs_acceleration_divergence_max_difference_Msun_pc3']<1e-5
assert out['declared_baryon_minimum_Msun_pc3']>0
assert out['shape_source_substitution_minimum_Msun_pc3']>0
(H/'source-density-verification.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))

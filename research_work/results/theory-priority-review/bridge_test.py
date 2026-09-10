"""Check whether the drafted gravity fit independently predicts photon loss.

This evaluates an exact parameter degeneracy on all retained SPARC rows.
No new fit, held-out selection, or physical coupling is introduced.
"""
from pathlib import Path
import hashlib,json
import numpy as np
H=Path(__file__).resolve().parent;J=H.parent/'joint-galaxy-audit'
save=lambda name,obj:(H/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
r=json.loads((J/'results.json').read_text())
data=json.loads((J/'galaxy-rotation-predictions.json').read_text())
b=[d for d in data if d['model']=='baryons'];g=[d for d in data if d['model']=='power']
assert [(d['galaxy'],d['R_kpc'],d['split']) for d in b]==[(d['galaxy'],d['R_kpc'],d['split']) for d in g]
KPC=3.085677581491367e19
R=np.array([d['R_kpc'] for d in b]);vb=np.array([d['predicted_kms'] for d in b])
gb=vb**2*1e6/(R*KPC)
par=r['sparc']['parameters'];A=par['A'];p=par['p'];astar=par['a_star_m_s2']
alpha=r['redshift']['parameters']['alpha_per_Mpc']
def predict(A,astar):return np.sqrt(vb**2+A*astar*(gb/astar)**p*R*KPC/1e6)
v=predict(A,astar)
assert np.max(abs(v-np.array([d['predicted_kms'] for d in g])))<1e-9
scans=[]
for factor in [.01,.1,1,10,100]:
    adjustedA=A*factor**(p-1)
    vv=predict(adjustedA,astar*factor)
    scans.append({'alpha_multiplier':factor,'alpha_per_Mpc':alpha*factor,'adjusted_A':adjustedA,
       'max_rotation_prediction_change_kms':float(np.max(abs(vv-v))),
       'predicted_z_at_100_Mly':float(np.expm1(alpha*factor*30.660139))})
extra=A*astar*(gb/astar)**p
base=extra*R*KPC/(2e6*v*v)
# Derivatives of ln v with respect to ln A, ln alpha and p.
Jac=np.c_[base,(1-p)*base,base*np.log(gb/astar)]
singular=np.linalg.svd(Jac,compute_uv=False)
assert singular[-1]/singular[0]<1e-12
assert max(s['max_rotation_prediction_change_kms'] for s in scans)<1e-9
out={'rows':len(b),'galaxies':len(set(d['galaxy'] for d in b)),
    'baseline_parameters':par,'scan':scans,'rotation_Jacobian_singular_values':singular.tolist(),
    'rotation_Jacobian_rank_at_relative_tolerance_1e_10':int(np.sum(singular>singular[0]*1e-10)),
    'interpretation':'Rotation constrains A*(c^2*alpha)^(1-p), not alpha independently while A is free. Redshift data can constrain alpha separately; the combined datasets are not claimed to have the same rank deficiency.',
    'physical_gap':'No independently derived coupling/capture/storage law currently fixes A. A successful gravity fit does not establish photon conversion as its source.',
    'inputs_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [J/'results.json',J/'galaxy-rotation-predictions.json']}}
save('bridge-identifiability.json',out)
print(json.dumps(out,indent=2))

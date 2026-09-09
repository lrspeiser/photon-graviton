from pathlib import Path
import json,itertools,math
P=Path(__file__).resolve().parent
# Leading nonrelativistic Coulomb/Pauli toy: eps=n^a, mu=n^(2-a), both masses=n^b, e and h fixed.
rows=[]
for a,b in itertools.product([0.,1.,2.],[0.,1.,2.]):
 s=b-2*a;l=a-b
 cav=b-a-1;hf_ratio=2-2*a
 rows.append(dict(a=a,b=b,optical_exponent=s,rod_exponent=l,cavity_to_optical_exponent=cav-s,hyperfine_to_optical_exponent=hf_ratio,observable_redshift_exponent=1+s))
protected=[r for r in rows if r['cavity_to_optical_exponent']==0 and r['hyperfine_to_optical_exponent']==0]
for r in protected:assert r['observable_redshift_exponent']==-r['rod_exponent']
# Exact operational identity, independent of leading atomic toy: n*L*nu_atom constant.
errors=[]
for ne,no,Le,Lo in itertools.product([.3,.5,1.],[1.,1.2],[1.,2.],[.5,1.]):
 nue=1/(ne*Le);nuo=1/(no*Lo)
 observed=(no/ne)*(nuo/nue)
 errors.append(abs(observed-Le/Lo))
gamma=7.7315e-11;p=.26390611655880836
out=dict(status='Analytic identities and synthetic checks, not new observational fits',leading_atomic_toy=rows,
 protected_leading_cases=protected,maximum_identity_error=max(errors),
 present_fixed_atom_fixed_rod_cavity_fractional_drift_per_year=-gamma,
 redshift_drift_per_year={str(z):gamma*((1+z)**(1-p)-(1+z)) for z in [.1,1.,2.,4.]},
 limitations=['Leading Coulomb/Pauli scaling does not include relativistic, nuclear, radiative or solid-state corrections.',
 'Exact identity assumes homogeneous nondispersive adiabatic propagation, identical cavity geometry with fixed mode, and rods tracking the chosen cavity length standard.',
 'A nonzero free-wave redshift with stable cavity/atom ratio then requires distances to change in those rod units for fixed coordinate separation.',
 'No general no-go theorem for dispersive, nonlocal, environment-dependent, multimode or nonadiabatic theories.'])
(P/'checks.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

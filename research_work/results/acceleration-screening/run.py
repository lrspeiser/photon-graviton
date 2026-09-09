from pathlib import Path
import sympy as s
import json,hashlib
HERE=Path(__file__).resolve().parent
y,B=s.symbols('Y B',nonnegative=True)
f=y/(1+y)
mu=1+B*s.diff(f,y)
long=s.simplify(mu+2*y*s.diff(mu,y))
assert s.simplify(mu-(1+B/(1+y)**2))==0
assert s.simplify(long-(1+B*(1-3*y)/(1+y)**3))==0
h=(long-1)/B
assert s.simplify(s.diff(h,y)-6*(y-1)/(1+y)**4)==0
cases=[dict(B=b,min_longitudinal=float(long.subs({B:b,y:1})),strictly_elliptic_all_accelerations=b<4) for b in [0,1,3,4,5,10]]
out=dict(protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),switch=str(f),gravity_coefficient=str(mu),longitudinal_eigenvalue=str(long),minimum_longitudinal='1-B/4 at Y=1 for B>0',spherical_force_direction='mu>=1 implies a<=a_Newton for fixed enclosed source mass',cases=cases,checks_passed=True)
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps(out,indent=2))

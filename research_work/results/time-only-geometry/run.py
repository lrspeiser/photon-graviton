from pathlib import Path
import sympy as s
import json,hashlib
HERE=Path(__file__).resolve().parent
t,x,y,z=s.symbols('t x y z',real=True);coords=[t,x,y,z]
q=s.Function('q')(t,x)
g=s.diag(-q*q,1,1,1);inv=g.inv()
Gamma=[[[s.simplify(sum(inv[a,d]*(s.diff(g[d,c],coords[b])+s.diff(g[d,b],coords[c])-s.diff(g[b,c],coords[d])) for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
Ricci=s.zeros(4)
for a in range(4):
 for b in range(4):
  Ricci[a,b]=s.simplify(sum(s.diff(Gamma[c][a][b],coords[c])-s.diff(Gamma[c][a][c],coords[b])+sum(Gamma[c][c][d]*Gamma[d][a][b]-Gamma[c][b][d]*Gamma[d][a][c] for d in range(4)) for c in range(4)))
R=s.simplify(sum(inv[a,b]*Ricci[a,b] for a in range(4) for b in range(4)))
Einstein=s.simplify(Ricci-g*R/2)
expected=s.diag(0,0,s.diff(q,x,2)/q,s.diff(q,x,2)/q)
assert s.simplify(Einstein-expected)==s.zeros(4)
profile=1-s.Rational(1,10)*s.exp(-x*x)
transverse=s.diff(profile,x,2)/profile
values={str(v):float(transverse.subs(x,v)) for v in [0,2]}
assert values['0']>0 and values['2']<0
out=dict(metric=str(g),ricci=str(Ricci),scalar_curvature=str(R),einstein_tensor=str(Einstein),normal_energy_density_in_GR='0',example_transverse_stress_in_units_c4_over_8piG=values,checks_passed=True,protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest())
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
print(json.dumps(out,indent=2))

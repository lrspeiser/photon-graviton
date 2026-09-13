"""Independent Jacobi initial-value check of lens-source distance factors."""
import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
HERE=Path(__file__).resolve().parent
fit=json.loads((HERE.parent/'brightness-distance-consistency/regular-area-results.json').read_text())
result=json.loads((HERE/'regular-optics-results.json').read_text())
f,q=sp.symbols('f q',positive=True)
H=-(1-f)*sp.log(1-f)*sp.sqrt(1+f/(1+q*f))
fun=sp.lambdify((f,q),sp.diff(H,f,2)/H,'numpy')
h=sp.lambdify((f,q),H,'numpy');rows=[]
for item in result['rows']:
    if item['model']!='attenuated':continue
    g=item['geometry'];zl=g['z_lens'];zs=g['z_source'];fl=zl/(1+zl);fs=zs/(1+zs)
    sol=solve_ivp(lambda x,y:[y[1],fun(x,fit['q'])*y[0]],(fl,fs),[0,1+zl],rtol=1e-10,atol=1e-12)
    assert sol.success
    ratio=float(sol.y[0,-1]/h(fs,fit['q']));error=abs(ratio/g['Dls_over_Ds']-1)
    assert error<1e-8
    rows.append(dict(Name=item['Name'],quadrature_ratio=g['Dls_over_Ds'],Jacobi_ODE_ratio=ratio,relative_error=error))
(HERE/'regular-optics-check.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
print('Maximum relative discrepancy:',max(r['relative_error'] for r in rows))

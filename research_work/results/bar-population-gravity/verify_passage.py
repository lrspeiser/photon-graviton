"""Independent adaptive integration of the local straight-passage formulas."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
rows=[]
L=1.;v=2.
for b in (.1,.01,.001):
    potential=-quad(lambda t:1/np.sqrt(b*b+(v*t)**2),0,L/v,epsabs=1e-11,epsrel=1e-11)[0]*v/L
    force=quad(lambda t:b/(b*b+(v*t)**2)**1.5,0,L/v,epsabs=1e-10,epsrel=1e-11)[0]*v/L
    p=-np.arcsinh(L/b)/L;a=1/(b*np.sqrt(b*b+L*L))
    ep=abs(potential-p)/abs(p);ea=abs(force-a)/abs(a)
    assert max(ep,ea)<1e-10
    rows.append(dict(b=b,potential=potential,force=force,potential_relative_error=ep,force_relative_error=ea))
Path(__file__).with_name('passage-checks.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf8',newline='\n')
print(rows)

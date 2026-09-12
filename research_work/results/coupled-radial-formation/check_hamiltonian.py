"""Independent finite-difference check of the regularized shell force."""
from pathlib import Path
import json
import numpy as np

rng=np.random.default_rng(619)
x=np.linspace(.1,3,24)*rng.choice([-1,1],24)
m=rng.uniform(.001,.01,len(x))
results=[]
for epsilon in (.05,.025,.00625):
    def U(x):
        radii=abs(x)
        kernel=1/np.sqrt(np.maximum(radii[:,None],radii[None,:])**2+epsilon**2)
        return -.5*float(m@kernel@m)-float(m@(1/np.sqrt(1+x*x)))
    order=np.argsort(abs(x));r=abs(x[order]);enclosed=np.cumsum(m[order])-.5*m[order]
    a=np.empty(len(x));a[order]=-np.sign(x[order])*enclosed*r/(r*r+epsilon*epsilon)**1.5
    a-=x/(1+x*x)**1.5
    numerical=[]
    for i in range(len(x)):
        h=1e-5;plus=x.copy();minus=x.copy();plus[i]+=h;minus[i]-=h
        numerical.append(-(U(plus)-U(minus))/(2*h*m[i]))
    error=float(max(abs(np.array(numerical)-a)))
    results.append(dict(epsilon=epsilon,max_absolute_force_error=error,passes=error<1e-7))
path=Path(__file__).resolve().parent/'hamiltonian-check.json'
path.write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
assert all(x['passes'] for x in results)
print(results)

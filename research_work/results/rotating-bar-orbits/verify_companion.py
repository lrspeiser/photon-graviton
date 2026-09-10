"""Verify cached geometry evaluation against the prior force implementation."""
import time
import numpy as np
from run import save
from field import companion
from cached_companion import CachedCompanion

rng=np.random.default_rng(3081);points=rng.uniform(-4,4,(48,3));rows=[]
for name in ['equatorial','caps','shell']:
    fast=CachedCompanion(name)
    p,a=companion(points,name);q,b=fast.evaluate(points)
    ep=float(np.max(abs(p-q)/np.maximum(abs(p),1)))
    ea=float(np.max(np.linalg.norm(a-b,axis=1)/np.maximum(np.linalg.norm(a,axis=1),1)))
    assert ep<1e-13 and ea<1e-13
    timing={}
    for kind,func in [('original',lambda:companion(points[:6],name)),('cached',lambda:fast.evaluate(points[:6]))]:
        t=time.perf_counter()
        for i in range(80):func()
        timing[kind]=(time.perf_counter()-t)/80
    rows.append(dict(geometry=name,max_potential_fraction_difference=ep,max_force_fraction_difference=ea,
                     seconds_per_six_points=timing))
save('companion-evaluation-check.json',rows)
print(rows)

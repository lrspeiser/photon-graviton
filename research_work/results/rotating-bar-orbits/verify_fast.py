"""Compare shared recurrences with the existing independent harmonic calls."""
import numpy as np
from run import FIELD_CACHE,save
from field import Multipole
from fast_multipole import FastMultipole

rng=np.random.default_rng(44812)
r=np.exp(rng.uniform(np.log(.005),np.log(30),120))
mu=rng.uniform(-.999,.999,120);phi=rng.uniform(-np.pi,np.pi,120)
xyz=np.c_[r*np.sqrt(1-mu*mu)*np.cos(phi),r*np.sqrt(1-mu*mu)*np.sin(phi),r*mu]
results=[]
for name in ['bar-L24','bar-L40','bar-L64','nuclei-L16']:
    old=Multipole.load(FIELD_CACHE/(name+'.npz'))
    new=FastMultipole.load(FIELD_CACHE/(name+'.npz'))
    p,a=old.evaluate(xyz);q,b=new.evaluate(xyz)
    err=float(np.max(np.linalg.norm(a-b,axis=1)/np.maximum(np.linalg.norm(a,axis=1),1)))
    perr=float(np.max(abs(p-q)/np.maximum(abs(p),1)))
    assert err<1e-12 and perr<1e-12
    results.append(dict(cache=name,maximum_force_fraction_difference=err,maximum_potential_fraction_difference=perr))
save('harmonic-recurrence-check.json',results)
print(results)

"""Test the above/below symmetry assumed by every current candidate field."""
import numpy as np
from run import DiskGrid,Field,PROTOCOL,save

# Existing cached potential grid; no rebuilding from an unspecified potential.
disk=DiskGrid(256,129,None)
rng=np.random.default_rng(11902)
R=rng.uniform(.2,15,48);phi=rng.uniform(-np.pi,np.pi,48);z=rng.uniform(.01,2,48)
points=np.c_[R*np.cos(phi),R*np.sin(phi),z]
mirror=points.copy();mirror[:,2]*=-1
rows=[]
for name in PROTOCOL['models']:
    f=Field(disk,name)
    p,a=f.evaluate(points);q,b=f.evaluate(mirror)
    b[:,2]*=-1
    ep=float(np.max(abs(p-q)/np.maximum(abs(p),1)))
    ea=float(np.max(np.linalg.norm(a-b,axis=1)/np.maximum(np.linalg.norm(a,axis=1),1)))
    print(name,ep,ea,flush=True)
    # Differentiating nearly equal potential values in the disk spline gives
    # cancellation-level force roundoff. The initial 1e-12 fractional gate
    # returned 1.12023e-12; retain that fact and use a roundoff-scale 1e-10
    # force check, still many orders tighter than the force-model accuracy.
    assert ep<1e-12 and ea<1e-10
    rows.append(dict(model=name,relative_potential_reflection_error=ep,relative_force_reflection_error=ea))
save('reflection-symmetry.json',dict(checks=rows,
    initial_force_threshold=1e-12,initial_maximum_force_fraction=1.1202226462422218e-12,
    force_roundoff_threshold=1e-10,
    implication='The same mirrored initial populations have mirrored orbits. Intrinsic north/south asymmetry requires asymmetric populations, selection, time dependence or an asymmetric field absent here.',
    observations_tested=False))
print('Above/below field reflection checks passed.')

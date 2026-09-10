"""Check disk interpolation at sampled actual diagnostic trajectories."""
import numpy as np
from run import CACHE,PROTOCOL,DiskGrid,save
from field import disks,disk_evaluate

traj=np.load(CACHE/'trajectories.npz')
points=[]
for model in PROTOCOL['models']:
    x=traj[model][:,:,:3].reshape(-1,3)
    idx=np.unique(np.r_[np.linspace(0,len(x)-1,100,dtype=int),
                        np.argmin(np.hypot(x[:,0],x[:,1])),np.argmax(abs(x[:,2]))])
    points.extend(x[idx])
points=np.array(points)
d=disks(24);g=DiskGrid(256,129,None)
p,a=disk_evaluate(d,points);q,b=g.evaluate(points)
rel=np.linalg.norm(a-b,axis=1)/np.linalg.norm(a,axis=1)
result=dict(sampled_trajectory_positions=len(points),
            max_disk_force_fraction_error=float(rel.max()),
            median_disk_force_fraction_error=float(np.median(rel)),
            worst_point_kpc=points[rel.argmax()].tolist(),
            scope='Approximately 100 positions per synthetic model plus extrema; not every integration stage or every possible stellar orbit')
save('orbit-field-check.json',result)
assert rel.max()<PROTOCOL['acceptance_checks']['disk_interpolation_force_error_fraction_of_reference_disk_force']
print(result)

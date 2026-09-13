"""Full directional catalog volume/aperture integration with scrambled Sobol checks."""
from pathlib import Path
import ast,json,hashlib
import numpy as np
from scipy.stats import qmc

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]

def definitions(path,stop):
    tree=ast.parse(path.read_text());body=[]
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id==stop for t in node.targets):break
        body.append(node)
    ns={'__file__':str(path)}
    exec(compile(ast.Module(body=body,type_ignores=[]),str(path),'exec'),ns)
    return ns

point=definitions(ROOT/'research_work/results/cluster-projected-deposition/run.py','axis')
finite=definitions(ROOT/'research_work/results/cluster-m87-finite-emitter/project.py','coarse')
radial=finite['calculate'](600)
from scipy.interpolate import PchipInterpolator
interp=PchipInterpolator(np.log(radial[0]),np.log(radial[1]))
L=float(point['lum'][np.argmin(np.linalg.norm(point['positions'],axis=1))]);a=point['a'];k=point['k']
pilot=json.loads((ROOT/'research_work/results/cluster-catalog-pilot/results.json').read_text())['runs'][0]
expected=pilot['internal_power_Lsun']+pilot['external_power_Lsun']-finite['prior']['point_emitter_power_Lsun']+finite['prior']['extended_emitter_power_Lsun']
apertures=np.array([.1,.3,.5,.9,1.])
runs=[]
for seed in [4701,4702,4703,4704]:
    unit=qmc.Sobol(3,scramble=True,seed=seed).random_base2(17)
    r=unit[:,0]**(1/3);z=2*unit[:,1]-1;phi=2*np.pi*unit[:,2]
    xyz=r[:,None]*np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
    values=np.concatenate([point['density'](chunk) for chunk in np.array_split(xyz,8)])
    old=k*L*a*np.exp(-a*r)*(-np.expm1(-(k-a)*r))/((k-a)*4*np.pi*r*r)
    new=np.exp(interp(np.log(np.maximum(r,1e-8))))
    values=values-old+new
    assert np.all(values>0) and np.all(np.isfinite(values))
    b=np.linalg.norm(xyz[:,:2],axis=1)
    for exponent in [15,17]:
        n=2**exponent; power=values[:n]*(4*np.pi/3)
        p3=[float(np.mean(power*(r[:n]<=t))) for t in apertures]
        p2=[float(np.mean(power*(b[:n]<=t))) for t in apertures]
        runs.append(dict(seed=seed,exponent=exponent,sphere_power_Lsun=p3,projected_power_Lsun=p2))
    print(json.dumps(dict(seed=seed,total=runs[-1]['sphere_power_Lsun'][-1],reference=expected)),flush=True)

fine=[v for v in runs if v['exponent']==17];coarse=[v for v in runs if v['exponent']==15]
p3=np.array([v['sphere_power_Lsun'] for v in fine]);p2=np.array([v['projected_power_Lsun'] for v in fine])
mean3=p3.mean(axis=0);mean2=p2.mean(axis=0)
assert np.all(mean2>=mean3)
assert np.all(np.diff(mean3)>0) and np.all(np.diff(mean2)>0)
err=max(abs(p3[:,-1]/expected-1))
assert err<.005
rows=[dict(radius_Mpc=float(b),mean_sphere_power_Lsun=float(s),mean_projected_power_Lsun=float(p),
    sphere_seed_min_Lsun=float(p3[:,i].min()),sphere_seed_max_Lsun=float(p3[:,i].max()),
    projected_seed_min_Lsun=float(p2[:,i].min()),projected_seed_max_Lsun=float(p2[:,i].max()),
    conditional_mean_radial_deflection_over_monopole_v2_c2=float(4*p/s),
    ratio_seed_min=float((4*p2[:,i]/p3[:,i]).min()),ratio_seed_max=float((4*p2[:,i]/p3[:,i]).max())) for i,(b,s,p) in enumerate(zip(apertures,mean3,mean2))]
result=dict(scope='814-source conditional power apertures. M87 extended; others point emitters. Randomized quadrature scatter is numerical, not observational. No fitted histories, lensing or motion data.',
    expected_total_power_Lsun=expected,max_fine_total_relative_error=float(err),
    coarse_mean_sphere_power_Lsun=np.mean([v['sphere_power_Lsun'] for v in coarse],axis=0).tolist(),
    coarse_mean_projected_power_Lsun=np.mean([v['projected_power_Lsun'] for v in coarse],axis=0).tolist(),
    catalog_sha256=hashlib.sha256((ROOT/'companion_causal_test/data/themis.dat').read_bytes()).hexdigest(),
    rows=rows,runs=runs)
(OUT/'results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2),flush=True)

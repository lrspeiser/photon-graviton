"""Check isotropic versus circular-orbit support in the fixed comparison well."""
from pathlib import Path
import ast,json
from functools import lru_cache
import numpy as np
from scipy.integrate import quad,solve_ivp
from scipy.special import roots_legendre
OUT=Path(__file__).resolve().parent
tree=ast.parse((OUT/'run.py').read_text());body=[]
for node in tree.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='rows' for t in node.targets):break
    body.append(node)
ns={'__file__':str(OUT/'run.py')};exec(compile(ast.Module(body=body,type_ignores=[]),'capture definitions','exec'),ns)
ns['roots_legendre']=lru_cache(None)(roots_legendre)
g=lambda r:r/(1+r*r)**1.5
rows=[]
for A in [.1,1.,10.,100.]:
    rho=lambda r:ns['profile'](r,A,256)
    checks=[]
    for r in [.1,.01,.001]:
        pressure=quad(lambda s:rho(s)*g(s),r,np.inf,epsabs=1e-12,epsrel=1e-9)[0]
        variance=pressure/rho(r)
        escape2=2/np.sqrt(1+r*r)
        checks.append(dict(radius_over_a=r,isotropic_radial_variance=variance,escape_speed_squared=escape2,variance_over_escape_squared=variance/escape2))
    assert checks[-1]['variance_over_escape_squared']>1
    assert checks[-1]['isotropic_radial_variance']>90*checks[-2]['isotropic_radial_variance']
    rows.append(dict(A=A,isotropic_checks=checks))

orbits=[]
for r in [.3,1.,3.,10.]:
    v=np.sqrt(r*g(r));period=2*np.pi*r/v
    def rhs(t,y):
        pos=y[:2];acc=-pos/(1+pos@pos)**1.5
        return np.r_[y[2:],acc]
    sol=solve_ivp(rhs,[0,10*period],[r,0,0,v],method='DOP853',rtol=1e-11,atol=1e-12,t_eval=np.linspace(0,10*period,1001))
    assert sol.success
    radius=np.linalg.norm(sol.y[:2],axis=0)
    energy=.5*np.sum(sol.y[2:]**2,axis=0)-1/np.sqrt(1+radius**2)
    err=max(abs(radius/r-1));drift=max(abs(energy/energy[0]-1))
    assert err<1e-8 and drift<1e-8
    epicycle=(r*r+4)/(1+r*r)**2.5
    assert epicycle>0 and v*v<2/np.sqrt(1+r*r)
    orbits.append(dict(radius_over_a=r,period=period,max_radius_relative_error=err,max_energy_relative_drift=drift,radial_epicyclic_frequency_squared=epicycle))
result=dict(scope='Fixed Plummer potential, nonrelativistic test-particle support. No self-gravity, capture-to-orbit mechanism, field binding or collective stability test.',isotropic=rows,circular_orbits=orbits)
(OUT/'support-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(dict(isotropic_center_ratios=[x['isotropic_checks'][-1]['variance_over_escape_squared'] for x in rows],circular_orbits=orbits),indent=2))

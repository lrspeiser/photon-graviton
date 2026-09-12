from pathlib import Path
import sys,json
import numpy as np
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'depleted-capture-feedback'))
from run import Model

def run(n,nq,dt):
    model=Model(100.,n,8)
    rho=np.load(HERE.parent/f'deposit-retention/profile-k0.1-n{n}.npz')['rho']
    dm,wd,_=model.gravity(rho);psi=model.wb+wd
    k=.1*(psi/(1+psi))**6
    captured=model.absorption(k)[0]/.1
    q=captured/model.V
    source_before=np.cumsum(captured)-captured
    before=np.cumsum(dm)-dm
    A=before-4*np.pi/3*rho*model.a**3
    shell=2*np.pi*rho*(model.b**2-model.a**2)
    outside=np.cumsum(shell[::-1])[::-1]-shell
    def field(x):
        r=np.abs(x);idx=np.clip(np.searchsorted(model.edges,r,side='right')-1,0,n-1)
        inside=r<=100.;safe=np.maximum(r,1e-100)
        enclosed=A[idx]+4*np.pi/3*rho[idx]*r**3
        gp=np.where(inside,enclosed/safe**2,dm.sum()/safe**2)
        depth=np.where(inside,enclosed/safe+2*np.pi*rho[idx]*(model.b[idx]**2-r**2)+outside[idx],dm.sum()/safe)
        acceleration=-np.sign(x)*(r/(1+r*r)**1.5+gp)
        return acceleration,-1/np.sqrt(1+r*r)-depth
    z,w=roots_legendre(nq)
    r=(model.a[:,None]**3+(model.b**3-model.a**3)[:,None]*(1+z)/2)**(1/3)
    cumulative=source_before[:,None]+4*np.pi/3*q[:,None]*(r**3-model.a[:,None]**3)
    beta=k[:,None]*cumulative/(4*np.pi*r*r*q[:,None])
    weights=(captured[:,None]*w[None,:]/2*np.sqrt(1-beta*beta)).ravel()
    source_rate=float(weights.sum());weights/=source_rate
    x=r.ravel();v=-1000*beta.ravel()
    acc,pot=field(x);energy0=.5*v*v+pot
    assert max(energy0)<0
    radii=np.array([.1,1.,3.,10.])
    def moments():
        fractions=(np.abs(x)[:,None]<radii[None,:]).T@weights
        ratio=np.minimum(radii[None,:]/np.maximum(abs(x)[:,None],1e-100),1.)
        projected=(1-np.sqrt(1-ratio*ratio)).T@weights
        energy=.5*v*v+field(x)[1]
        return np.r_[fractions,projected,float(weights@(.5*v*v)),float(weights@energy)]
    baseline=moments();previous=baseline.copy();integral=np.zeros(10)
    max_energy_error=0.;outputs=[]
    for step in range(1,round(30/dt)+1):
        v+=.5*dt*acc;x+=dt*v
        acc,pot=field(x);v+=.5*dt*acc
        now=moments();integral+=dt*(previous+now)/2;previous=now
        max_energy_error=max(max_energy_error,float(max(abs(.5*v*v+pot-energy0))))
        if step in [round(T/dt) for T in (1,10,30)]:
            T=step*dt
            outputs.append(dict(T=T,enclosed_fractions=(integral[:4]/T).tolist(),
                                projected_fractions=(integral[4:8]/T).tolist(),
                                mean_kinetic=float(integral[8]/T),mean_mechanical=float(integral[9]/T),
                                added_rest_mass_per_unit_source_amplitude=source_rate*T))
    mass_error=abs(float(weights.sum())-1)
    assert max_energy_error<1e-4 and mass_error<1e-10
    return dict(shells=n,launch_nodes=nq,dt=dt,source_rest_rate=source_rate,
                source_enclosed_fractions=baseline[:4].tolist(),source_projected_fractions=baseline[4:8].tolist(),source_mean_kinetic=float(baseline[8]),
                source_mean_mechanical=float(baseline[9]),max_orbit_energy_error=max_energy_error,
                mass_normalization_error=mass_error,outputs=outputs)

cases=[]
for args in ((512,2,.004),(512,2,.002),(1024,4,.002)):
    x=run(*args);cases.append(x);print(json.dumps(x),flush=True)
checks=[]
for label,a,b in (('time',cases[0],cases[1]),('combined',cases[1],cases[2])):
    errors=[]
    for x,y in zip(a['outputs'],b['outputs']):
        f=float(max(abs(np.array(x['enclosed_fractions'])-y['enclosed_fractions'])))
        projected=float(max(abs(np.array(x['projected_fractions'])-y['projected_fractions'])))
        kinetic=abs(x['mean_kinetic']/y['mean_kinetic']-1)
        errors.append(dict(T=x['T'],fraction_change=f,projected_fraction_change=projected,kinetic_relative_change=kinetic,passes=f<.005 and projected<.005 and kinetic<.01))
    checks.append(dict(kind=label,comparisons=errors))
(HERE/'results.json').write_text(json.dumps(dict(cases=cases,checks=checks,normalization='Linear test-population response: multiply source rate and mass coefficients by an infinitesimal amplitude eta. Coefficients are not a finite self-gravitating added mass.'),indent=2)+'\n',encoding='utf8',newline='\n')
assert all(c['passes'] for x in checks for c in x['comparisons']),'Retain failed numerical gate'

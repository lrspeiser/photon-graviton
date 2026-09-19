"""Exploratory sign, power, clock and constant changes; no retrospective selection."""
import itertools
from types import SimpleNamespace
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from common import save
from wild_model import WildModel
from two_dimensional import rays,bend_metrics

def homogeneous(config,em=1.):
    g=config.get('g',.2);b=config.get('b',2);p=config.get('power',1)
    mu2=config.get('mu2',.16);rest=config.get('rest',0.)
    clock=.003+4*rest;offset=mu2**2/.4 if mu2<0 else 0.
    def F(x):return g*x**p
    def Fp(x):return g*p*x**(p-1)
    def rhs(t,y):
        phi,pi,_,_=y;a=np.exp(-F(phi));q=np.exp(-b*F(phi))
        return [pi,-mu2*phi-.1*phi**3+Fp(phi)*(em*a+b*clock*q),a,q]
    sol=solve_ivp(rhs,[0,40],[0,0,0,0],method='DOP853',rtol=2e-11,atol=2e-13,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    time=np.linspace(0,40,2001);phi,pi,dist,tau=sol.sol(time)
    photon=em*np.exp(-F(phi));material=clock*np.exp(-b*F(phi))
    receiving=.5*pi*pi+.5*mu2*phi*phi+.025*phi**4+offset
    energy=photon+material+receiving;E0=em+clock+offset
    def arrival(te):
        target=sol.sol(te)[2]+5
        if target>sol.sol(40)[2]:return None
        return brentq(lambda t:sol.sol(t)[2]-target,te,40,xtol=1e-12)
    to=arrival(2.);factor=None;timing=None;error=None
    if to is not None and arrival(2.0003) is not None:
        delta=F(sol.sol(to)[0])-F(sol.sol(2.)[0])
        factor=np.exp((1-b)*delta)
        finite=[]
        for h in [3e-4,1.5e-4]:
            end=arrival(2.+h)
            finite.append((sol.sol(end)[3]-sol.sol(to)[3])/(sol.sol(2.+h)[3]-sol.sol(2.)[3]))
        timing=2*finite[1]-finite[0];error=abs(timing/factor-1)
    composition=(rest+.001)/(rest+.0005)-1
    nonzero_force=b!=0 and np.max(abs(Fp(phi)))>1e-12
    return dict(config=config,initial_em=em,initial_material=clock,initial_receiving=offset,
                energy_relative_error=np.max(abs(energy-E0))/E0,
                min_receiving_energy=np.min(receiving),max_phi=np.max(abs(phi)),
                final_em=photon[-1],final_material=material[-1],final_receiving=receiving[-1],
                arrival=to,spectral_stretch=factor,event_stretch=timing,timing_identity_error=error,
                fixed_ruler_speed_change=np.max(abs(np.exp((b-1)*F(phi))-1)),
                force_fractional_composition_difference=composition,nonzero_force=bool(nonzero_force),
                history=dict(time=time[::10],phi=phi[::10],em=photon[::10],
                             material=material[::10],receiving=receiving[::10]))

def derivative_check(config):
    rng=np.random.default_rng(221919)
    model=WildModel(n=24,L=12,dim=1,**config)
    for key,scale in [('A',.1),('D',.1),('phi',.3),('Pi',.03)]:
        setattr(model,key,rng.normal(0,scale,model.shape))
    model.X+=.117;model.P=rng.normal(0,.02,model.P.shape)
    state=model.state();rhs=model.rhs();errs=[]
    for _ in range(4):
        ds={k:rng.normal(size=v.shape) for k,v in state.items() if k!='theta'}
        expected=model.dv*np.sum(-rhs['D']*ds['A']+rhs['A']*ds['D']-rhs['Pi']*ds['phi']+rhs['phi']*ds['Pi'])
        expected+=np.sum(-rhs['P']*ds['X']+rhs['X']*ds['P'])
        eps=1e-6
        for k in ds:setattr(model,k,state[k]+eps*ds[k])
        plus=model.energy()
        for k in ds:setattr(model,k,state[k]-eps*ds[k])
        minus=model.energy();model.set_state(state)
        errs.append(abs((plus-minus)/(2*eps)-expected)/max(abs(expected),1e-3))
    return max(errs)

def spatial(config,out,label):
    model=WildModel(n=128,dim=2,L=32,**config)
    model.initialize_light(energy=1.,wavelength=3.)
    initial=model.sectors();E0=initial.sum();steps=int(np.ceil(18/(.08*model.dx/(np.sqrt(2)*max(1,model.v)))))
    dt=18/steps;history=[];error=0.
    for k in range(steps+1):
        if k%4==0 or k==steps:
            sec=model.sectors();error=max(error,abs(sec.sum()/E0-1))
            if k%max(4,steps//100//4*4)==0 or k==steps:history.append([k*dt,*sec])
        if k<steps:model.step(dt)
        if not np.all(np.isfinite(model.phi)):raise FloatingPointError(label)
    # For rays only, the generated optical phase is F(phi), not g*phi.
    rr=rays(SimpleNamespace(axis=model.axis,phi=model.F(model.phi),g=1.))
    result=dict(label=label,config=config,initial_sectors=initial,final_sectors=model.sectors(),
                energy_relative_error=error,history=history,rays=rr,bend=bend_metrics(rr),
                max_phi=np.max(abs(model.phi)),axis=model.axis[::2],phi=model.phi[::2,::2],
                counted_receiving_potential_offset=model.offset*model.L**2)
    save(out/(label+'.json'),result)
    return {k:v for k,v in result.items() if k not in ['history','axis','phi']}

def run(out):
    configs=[dict(g=g,power=p,b=b) for g,p,b in
             itertools.product([-.2,.2],[1,2,3,5],[-2,-1,0,.5,1,2,4])]
    configs += [dict(mu2=m) for m in [-.16,0,.16,1.]]
    configs += [dict(rest=r,b=b) for r,b in itertools.product([.01,1.,100.],[1,2])]
    rows=[]
    for i,config in enumerate(configs):
        result=homogeneous(config);off=homogeneous(config,em=0)
        science=dict(redshift=result['spectral_stretch'] is not None and result['spectral_stretch']>1.001,
                     timing=result['timing_identity_error'] is not None and result['timing_identity_error']<.01,
                     fixed_ruler_speed=result['fixed_ruler_speed_change']<.001,
                     universal_force=result['nonzero_force'] and result['force_fractional_composition_difference']<.001,
                     photon_exclusive=off['max_phi']<1e-10)
        save(out/('homogeneous-%03d.json'%i),dict(on=result,source_off=off,science=science))
        rows.append(dict(config=config,**{k:v for k,v in result.items() if k not in ['config','history']},
                         source_off_max_phi=off['max_phi'],source_off_energy_error=off['energy_relative_error'],
                         mechanism_gates=science,passed=all(science.values())))
    print('EXPLORATORY homogeneous',len(rows),'plus source-off controls complete',flush=True)
    tests=[('baseline',{}),('sign-reversed',dict(g=-.2)),('clock-reversed',dict(b=-1)),
           ('quadratic',dict(power=2)),('cubic',dict(power=3)),('speed-two',dict(v=2)),
           ('negative-mass-squared',dict(mu2=-.16)),('material-rest',dict(rest=1.))]
    derivatives={label:derivative_check(config) for label,config in tests}
    spatial_rows=[]
    for label,config in tests:
        row=spatial(config,out,label);spatial_rows.append(row)
        print('EXPLORATORY 2D',label,'energy',row['energy_relative_error'],'bend',row['bend']['max_bend'],flush=True)
    sign=max(abs(a['angle']-b['angle']) for a,b in zip(spatial_rows[0]['rays'],spatial_rows[1]['rays']))
    gates=dict(homogeneous_energy=all(max(r['energy_relative_error'],r['source_off_energy_error'])<1e-8 for r in rows),
               homogeneous_timing=all(r['timing_identity_error'] is None or r['timing_identity_error']<1e-6 for r in rows),
               hamiltonian_derivatives=max(derivatives.values())<1e-5,
               spatial_energy=all(r['energy_relative_error']<1e-3 for r in spatial_rows),
               sign_observable_symmetry=sign<1e-8)
    return dict(stage='explore',homogeneous=rows,spatial=spatial_rows,derivatives=derivatives,
                sign_max_absolute_ray_difference=sign,gates=gates,numerical_pass=all(gates.values()),
                mechanism_pass=any(r['passed'] for r in rows),
                note='Exploratory modified laws. Homogeneous screening cannot certify a spatial or observational theory; no historical novelty asserted.')

"""Analytic mode reductions and independent Hamiltonian controls; prior art credited."""
import itertools
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from model import Model

def homogeneous(g,b,m,em=1.,clock=.003):
    def rhs(t,y):
        ph,pi,distance,tau=y
        a=np.exp(-g*ph);q=np.exp(-b*g*ph)
        return [pi,g*em*a+b*g*clock*q-m*m*ph-.1*ph**3,a,q]
    sol=solve_ivp(rhs,[0,40],[0.,0.,0.,0.],dense_output=True,method='DOP853',rtol=2e-11,atol=2e-13)
    if not sol.success:raise RuntimeError(sol.message)
    times=np.linspace(0,40,1601);ph,pi,d,tau=sol.sol(times)
    energy=em*np.exp(-g*ph)+clock*np.exp(-b*g*ph)+.5*(pi*pi+m*m*ph*ph)+.025*ph**4
    te=2.
    def arrive(t):
        start=sol.sol(t)[2]
        return brentq(lambda tt:sol.sol(tt)[2]-start-5,t,40,xtol=1e-12)
    to=arrive(te)
    pe,po=sol.sol(te)[0],sol.sol(to)[0]
    carrier=np.exp(g*(po-pe))
    qratio=np.exp(-b*g*(po-pe))
    measured=carrier*qratio
    finite=[]
    for h in [3e-4,1.5e-4]:
        end=arrive(te+h)
        finite.append((sol.sol(end)[3]-sol.sol(to)[3])/(sol.sol(te+h)[3]-sol.sol(te)[3]))
    extrap=2*finite[1]-finite[0]
    speed=np.exp((b-1)*g*ph)
    return dict(g=g,b=b,m=m,initial_em=em,initial_clock=clock,
                relative_energy_error=np.max(abs(energy-em-clock))/max(em+clock,1e-30),
                arrival=to,coordinate_carrier_factor=carrier,clock_ratio=qratio,
                measured_spectral_factor=measured,measured_arrival_factor=extrap,
                timing_identity_error=abs(extrap/measured-1),local_speed_max_change=np.max(abs(speed-1)),
                source_off_max_phi=np.max(abs(ph)) if em==0 else None,
                final_em=em*np.exp(-g*ph[-1]),final_receiving=.5*(pi[-1]**2+m*m*ph[-1]**2)+.025*ph[-1]**4,
                final_clock=clock*np.exp(-b*g*ph[-1]),
                history=dict(time=times[::8],phi=ph[::8],em=(em*np.exp(-g*ph))[::8],
                             clock=(clock*np.exp(-b*g*ph))[::8],energy=energy[::8]),
                numerical_pass=np.max(abs(energy-em-clock))/max(em+clock,1e-30)<1e-9 and abs(extrap/measured-1)<1e-6)

def derivative_checks():
    rng=np.random.default_rng(190919)
    model=Model(n=32,L=12,dim=1)
    model.A=rng.normal(0,.1,model.shape);model.D=rng.normal(0,.1,model.shape)
    model.phi=rng.normal(0,.03,model.shape);model.Pi=rng.normal(0,.03,model.shape)
    model.X+=.113;model.P=rng.normal(0,.02,model.P.shape)
    initial=model.state();rhs=model.rhs();errors=[]
    for _ in range(8):
        direction={k:rng.normal(size=v.shape) for k,v in initial.items() if k!='theta'}
        pred=model.dv*np.sum(-rhs['D']*direction['A']+rhs['A']*direction['D']
                            -rhs['Pi']*direction['phi']+rhs['phi']*direction['Pi'])
        pred+=np.sum(-rhs['P']*direction['X']+rhs['X']*direction['P'])
        eps=1e-6
        for k,v in direction.items():setattr(model,k,initial[k]+eps*v)
        plus=model.energy()
        for k,v in direction.items():setattr(model,k,initial[k]-eps*v)
        minus=model.energy();model.set_state(initial)
        errors.append(abs((plus-minus)/(2*eps)-pred)/max(abs(pred),1e-3))
    energy=model.energy()
    for _ in range(40):model.step(.005)
    for _ in range(40):model.step(-.005)
    reverse=max(np.max(abs(model.state()[k]-v)) for k,v in initial.items())
    # Explicit negative control: a counted receiving wave changes EM energy,
    # but the reciprocal EM source is deliberately omitted.
    bad=Model(n=64,L=12,g=.5,b=0,clock_energy=0)
    bad.A=.1*np.cos(2*np.pi*bad.coords[0]/3);bad.D=.1*np.sin(2*np.pi*bad.coords[0]/3)
    bad.phi=.2*np.ones(bad.shape);bad.Pi=.1*np.ones(bad.shape)
    bad.omit_em_reciprocity=True;e0=bad.energy()
    for _ in range(400):bad.step(.005)
    defect=abs(bad.energy()/e0-1)
    return dict(max_gradient_relative_error=max(errors),time_reversal_error=reverse,
                missing_reciprocity_energy_defect=defect,
                passed=max(errors)<1e-5 and reverse<1e-10 and defect>1e-3)

def run():
    rows=[homogeneous(g,b,m) for g,b,m in itertools.product([.05,.2,.5],[0,1,2],[0,.4,1])]
    controls=[homogeneous(.2,b,.4,em=0) for b in [0,1,2]]
    controls+=[homogeneous(0,2,.4),homogeneous(.2,2,.4,clock=0)]
    checks=derivative_checks()
    # At the same phi gradient, acceleration is proportional to internal I/M.
    composition={str(b):dict(acceleration_ratio=2 if b else None,
                            universality_pass=False,explanation='Zero force' if b==0 else 'Force per mass changes by factor two')
                 for b in [0,1,2]}
    scientific=[]
    for r in rows:
        gates=dict(redshift=r['measured_spectral_factor']>1.001,
                   timing=r['timing_identity_error']<.01,
                   fixed_ruler_speed=r['local_speed_max_change']<.001,
                   universal_nonzero_matter_force=False)
        scientific.append(dict(g=r['g'],b=r['b'],m=r['m'],gates=gates,passed=all(gates.values())))
    return dict(stage='small',rows=rows,controls=controls,derivatives=checks,composition=composition,
                numerical_pass=all(r['numerical_pass'] for r in rows+controls) and checks['passed'],
                scientific=scientific,mechanism_pass=any(r['passed'] for r in scientific))

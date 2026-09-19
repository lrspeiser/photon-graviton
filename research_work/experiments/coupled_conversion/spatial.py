"""Run prepared finite-energy field fixtures; no cosmological or halo input."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.signal import hilbert
from model import Model

def pulse_metrics(series,reverse=False):
    t=np.array(series['time']);signal=np.array(series['signal']);q=np.array(series['clock_rate'])
    if np.max(abs(signal))<1e-12:return None
    analytic=hilbert(signal);power=abs(analytic)**2
    phase=np.unwrap(np.angle(analytic));frequency=np.gradient(phase,t)
    tau=cumulative_trapezoid(q,t,initial=0)
    rows=[]
    for lower,upper in [(12.,18.),(18.,25.5)]:
        mask=(t>=lower)&(t<=upper)
        mask&=power>.01*np.max(power[mask])
        if mask.sum()<8:return None
        w=power[mask];tm=t[mask];measure=q[mask]
        norm=np.trapezoid(w,tm)
        centroid=np.trapezoid(w*tau[mask],tm)/norm
        omega=np.trapezoid(w*frequency[mask]/measure,tm)/norm
        sigma=np.sqrt(np.trapezoid(w*(frequency[mask]/measure-omega)**2,tm)/norm)
        width=np.sqrt(np.trapezoid(w*(tau[mask]-centroid)**2,tm)/norm)
        rows.append(dict(clock_centroid=centroid,coordinate_centroid=np.trapezoid(w*tm,tm)/norm,
                         frequency=omega,frequency_std=sigma,pulse_width=width,signal_norm=norm))
    return dict(pulses=rows,clock_separation=rows[1]['clock_centroid']-rows[0]['clock_centroid'],
                frequency_mean=np.mean([r['frequency'] for r in rows]),
                phase_bandwidth_fraction=max(r['frequency_std']/abs(r['frequency']) for r in rows))

def simulate(n=1024,dim=1,g=.2,b=2,m=.4,v=1.,energy=1.,wavelength=2.,
             clock_energy=1.,reverse=False,step_factor=.08,free=False):
    L=48 if dim==1 else 32
    model=Model(n=n,dim=dim,L=L,g=g,b=b,m=m,v=v,clock_energy=clock_energy)
    model.initialize_light(energy=energy,wavelength=wavelength,reverse=reverse)
    if free:
        model.phi=.1*np.exp(-model.coords[0]**2/(2*1.5**2))*np.cos(2*np.pi*model.coords[0]/4)
        model.Pi=-v*model.grad(model.phi)[0]
        factor=1/np.sqrt(model.sectors()[1])
        model.phi*=factor;model.Pi*=factor
    initial=model.sectors();E0=initial.sum();mom0=model.momentum()
    T=26 if dim==1 else 18
    dt=step_factor*model.dx/np.sqrt(dim)
    steps=int(np.ceil(T/dt));dt=T/steps
    time=[];energies=[];momenta=[];retained=[];profile=[];directional=[]
    series=dict(time=[],signal=[],clock_rate=[])
    emax=0.;pmax=0.;min_a=1.;max_a=1.;max_clock_speed=0.
    observation=np.array([[(-6. if reverse else 6.)]]) if dim==1 else np.array([[8.,0.]])
    times_snap=[6.,10.,14.,18.] if dim==2 else [0.,6.,12.,18.,26.]
    snap_index=0
    for k in range(steps+1):
        t=k*dt
        if k%4==0 or k==steps:
            sec=model.sectors();mom=model.momentum()
            error=abs(sec.sum()-E0)/max(E0,1e-30)
            emax=max(emax,error);pmax=max(pmax,np.linalg.norm(mom-mom0)/max(initial[0],initial[1],initial[2],1e-30))
            a=np.exp(-g*model.phi);min_a=min(min_a,a.min());max_a=max(max_a,a.max())
            max_clock_speed=max(max_clock_speed,np.max(abs(np.exp((b-1)*g*model.phi)-1)))
            if k%max(4,steps//200//4*4)==0 or k==steps:
                time.append(t);energies.append(sec);momenta.append(mom)
                if dim==1:
                    grad_a=model.grad(model.A)[0]
                    directional.append([.25*np.sum(a*(model.D-sign*grad_a)**2)*model.dv for sign in [1,-1]])
                radius=np.sqrt(sum(c*c for c in model.coords))
                retained.append(np.sum(model.receiving_density()[radius<3])*model.dv)
            if dim==1:
                ph,_,_=model.interpolation(model.phi,observation)
                dd,_,_=model.interpolation(model.D,observation)
                grad,_,_=model.interpolation(model.grad(model.A)[0],observation)
                direction=-1 if reverse else 1
                signal=.5*np.exp(-g*ph[0])*(dd[0]-direction*grad[0])
                series['time'].append(t);series['signal'].append(signal);series['clock_rate'].append(np.exp(-b*g*ph[0]))
        if snap_index<len(times_snap) and t>=times_snap[snap_index]-dt/2:
            stride=max(1,n//256)
            if dim==1:
                profile.append(dict(time=t,x=model.axis[::stride],phi=model.phi[::stride],
                                    receiving_density=model.receiving_density()[::stride]))
            else:
                profile.append(dict(time=t,axis=model.axis[::max(1,n//64)],
                                    phi=model.phi[::max(1,n//64),::max(1,n//64)],
                                    receiving_density=model.receiving_density()[::max(1,n//64),::max(1,n//64)]))
            snap_index+=1
        if k<steps:model.step(dt)
        if not np.all(np.isfinite(model.phi)):raise FloatingPointError('Nonfinite receiving field')
    final=model.sectors()
    result=dict(config=dict(n=n,dim=dim,g=g,b=b,m=m,v=v,energy=energy,wavelength=wavelength,
                            clock_energy=clock_energy,reverse=reverse,step_factor=step_factor,free=free),
                initial_sectors=initial,final_sectors=final,sector_order=['electromagnetic','receiving','material_internal','material_motion'],
                energy_relative_error=emax,momentum_error_over_initial_em=pmax,
                min_optical_speed=min_a,max_optical_speed=max_a,local_clock_speed_max_change=max_clock_speed,
                em_energy_lost=initial[0]-final[0],material_internal_energy_lost=initial[2]-final[2],
                receiving_gain=final[1]-initial[1],
                material_positions=model.X,material_momenta=model.P,
                max_final_phi=np.max(abs(model.phi)),
                time=time,energy_history=energies,momentum_history=momenta,central_receiving_energy=retained,
                directional_energy_proxy=directional if dim==1 else None,
                snapshots=profile,receiver=series if dim==1 else None,
                pulse_metrics=pulse_metrics(series,reverse) if dim==1 and energy>0 else None)
    return result,model

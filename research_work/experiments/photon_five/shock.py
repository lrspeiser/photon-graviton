"""E4: conservative ordinary-gas collision, photon fuel, and local erasure."""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
GAMMA=5/3

def primitive(U):
    rho=U[0]; vel=U[1]/rho
    pressure=(GAMMA-1)*(U[2]-.5*rho*vel*vel)
    if np.min(rho)<=0 or np.min(pressure)<=0: raise ValueError('Gas positivity failed; no floors applied')
    return rho,vel,pressure

def fluxes(U):
    rho,vel,p=primitive(U)
    F=np.array([U[1],U[1]*vel+p,(U[2]+p)*vel])
    a=abs(vel)+np.sqrt(GAMMA*p/rho)
    pad=np.pad(U,((0,0),(1,1)),mode='edge')
    fp=np.pad(F,((0,0),(1,1)),mode='edge')
    speed=np.maximum(np.pad(a,(1,1),mode='edge')[:-1],np.pad(a,(1,1),mode='edge')[1:])
    flux=.5*(fp[:,:-1]+fp[:,1:])-.5*speed*(pad[:,1:]-pad[:,:-1])
    entropy=rho*np.log(p/rho**GAMMA)
    ent=np.pad(entropy,(1,1),mode='edge')
    entflux=np.pad(entropy*vel,(1,1),mode='edge')
    ef=.5*(entflux[:-1]+entflux[1:])-.5*speed*(ent[1:]-ent[:-1])
    return flux,ef,a

def update(U,dt,dx):
    rho,vel,p=primitive(U)
    entropy=rho*np.log(p/rho**GAMMA)
    F,ef,a=fluxes(U)
    new=U-dt/dx*np.diff(F,axis=1)
    adv=entropy-dt/dx*np.diff(ef)
    rn,vn,pn=primitive(new)
    residual=np.maximum(np.log(pn/rn**GAMMA)-adv/rn,0)/dt
    return new,residual,dt*(F[2,-1]-F[2,0]),vel

def evolve(n,eta,emitting=True):
    dx=12/n;x=-6+dx*(np.arange(n)+.5)
    left=np.exp(-.5*((x+2)/.3)**2);right=np.exp(-.5*((x-2)/.3)**2)
    rho=.01+left+right
    vel=.7*(left-right)/rho
    U=np.array([rho,rho*vel,.02/(GAMMA-1)+.5*rho*vel**2])
    u=np.zeros(n);fuel=.5;t=0.;boundary=0.;emitted=0.;returned=0.
    initial=U[2].sum()*dx+fuel
    maxledger=0.;maxfield=0.;minrho=rho.min();minpressure=.02
    entropy_total=0.;diffuse=0.;snapshots=[];nexttime=0.
    while t<6-1e-12:
        rho,vel,p=primitive(U)
        dt=min(.25*dx/np.max(abs(vel)+np.sqrt(GAMMA*p/rho)),6-t)
        U,rate,flux,oldvel=update(U,dt,dx)
        boundary+=flux
        entropy_total+=np.sum(rate*U[0])*dt*dx
        diffuse+=np.sum(rate*U[0]*(np.gradient(oldvel,dx)>=0))*dt*dx
        stars=np.array([-2+.7*(t+.5*dt),2-.7*(t+.5*dt)])
        profile=sum(np.exp(-.5*((x-s)/.2)**2) for s in stars)
        dose=min(.06*dt,fuel) if emitting else 0.
        q=profile/(profile.sum()*dx)*dose/dt
        decay=1/3+eta*rate
        fraction=-np.expm1(-decay*dt)
        unew=u*(1-fraction)+q*fraction/decay
        heat=u+q*dt-unew
        U[2]+=heat
        fuel-=dose;emitted+=dose;returned+=heat.sum()*dx;u=unew
        t+=dt
        rn,vn,pn=primitive(U)
        minrho=min(minrho,rn.min());minpressure=min(minpressure,pn.min())
        maxfield=max(maxfield,u.max())
        ledger=(U[2].sum()+u.sum())*dx+fuel+boundary
        maxledger=max(maxledger,abs(ledger-initial)/initial)
        if t>=nexttime or t>=6-1e-12:
            curvature=-gaussian_filter1d(u,.15/dx,order=2,mode='constant')/dx**2
            snapshots.append(dict(time=t,x=x[::max(1,n//256)],gas=rn[::max(1,n//256)],
                                  field=u[::max(1,n//256)],curvature=curvature[::max(1,n//256)]))
            nexttime+=2.
    curve=-gaussian_filter1d(u,.15/dx,order=2,mode='constant')/dx**2
    peaks,_=find_peaks(curve)
    chosen=[]
    for sign in [-1,1]:
        candidates=peaks[(x[peaks]*sign>0)&(curve[peaks]>0)]
        chosen.append(float(x[candidates[np.argmax(curve[candidates])]]) if len(candidates) else None)
    gas_peaks=[x[np.where(x*sign>0)[0][np.argmax(U[0,x*sign>0])]] for sign in [-1,1]]
    distance=max(abs(np.array(chosen)-np.array([-2.2,2.2]))) if None not in chosen else None
    return dict(grid=n,eta=eta,source_on=emitting,relative_energy_error=maxledger,
                min_density=minrho,min_pressure=minpressure,source_off_max_field=maxfield if not emitting else None,
                initial_total_energy=initial,final_gas_energy=U[2].sum()*dx,
                final_field_energy=u.sum()*dx,remaining_fuel=fuel,boundary_energy=boundary,
                donated=emitted,returned_heat=returned,entropy_residual_integral=entropy_total,
                noncompression_residual_fraction=diffuse/max(entropy_total,1e-30),
                optical_peaks=chosen,gas_peaks=gas_peaks,star_peak_distance=distance,snapshots=snapshots)

def run():
    rows=[]
    for eta in [0,1,5]:
        for n in [256,512,1024]:
            row=evolve(n,eta);rows.append(row)
            print('E4 collision',eta,n,row['relative_energy_error'],flush=True)
    off=evolve(256,5,False)
    U=np.array([np.ones(128),np.full(128,.2),np.full(128,.1/(GAMMA-1)+.02)])
    un,rate,_,_=update(U,.01,.1)
    uniform_error=np.max(abs(un-U))+np.max(rate)
    q=.07;t=2.;decay=.8;u=0.
    for _ in range(200):u=u*np.exp(-decay*.01)+q*(-np.expm1(-decay*.01))/decay
    analytic_error=abs(u-q*(-np.expm1(-decay*t))/decay)
    refinements=[]
    for eta in [0,1,5]:
        a,b=[r for r in rows if r['eta']==eta][-2:]
        refinements.append(max(abs(np.array(a['optical_peaks'])-b['optical_peaks'])))
    gates=dict(energy=max(r['relative_energy_error'] for r in rows+[off])<1e-8,
               positivity=all(r['min_density']>0 and r['min_pressure']>0 for r in rows+[off]),
               source_off=off['source_off_max_field']==0,uniform_entropy=uniform_error<1e-10,
               analytic_decay=analytic_error<1e-12,peak_refinement=max(refinements)<.15)
    baseline=rows[2]['star_peak_distance']
    mechanisms=[dict(eta=eta,pass_gate=all(r['star_peak_distance']<baseline and
                       max(abs(np.array(r['gas_peaks'])))<.75
                       for r in rows if r['eta']==eta and r['grid']==1024)) for eta in [1,5]]
    return dict(experiment='E4',gates=gates,numerical_pass=all(gates.values()),
                mechanism_pass=all(x['pass_gate'] for x in mechanisms),mechanisms=mechanisms,
                rows=rows,source_off=off,uniform_error=uniform_error,analytic_decay_error=analytic_error,
                peak_refinements=refinements,
                observational_status='Controlled 1D fixture only; numerical entropy contamination, full merger histories and coupled gravity unresolved.')

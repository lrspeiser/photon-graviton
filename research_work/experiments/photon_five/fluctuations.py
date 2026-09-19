"""E3: independently verified covariance and stochastic radial probe orbits."""
import numpy as np
from scipy.linalg import solve_continuous_lyapunov
from common import plain

def system(n):
    dx=4/(n+1); r=dx*np.arange(1,n+1)
    K=np.eye(n)*(1+2*(.3/dx)**2)
    K+=np.diag(np.full(n-1,-(.3/dx)**2),1)+np.diag(np.full(n-1,-(.3/dx)**2),-1)
    q=np.exp(-r*r/.5**2)
    D=np.sqrt(q[:,None]*q[None,:])*np.exp(-(r[:,None]-r[None,:])**2/(2*.15**2))
    eig,U=np.linalg.eigh(K)
    Dm=U.T@D@U
    Xmodal=Dm/(eig[:,None]+eig[None,:]+.5*(eig[:,None]-eig[None,:])**2)
    X=U@Xmodal@U.T
    V=(K@X+X@K)/2
    Y=(K@X-X@K)/2
    Sigma=np.block([[X,Y],[Y.T,V]])
    de,dv=np.linalg.eigh(D)
    B=dv*np.sqrt(np.maximum(de,0))[None,:]
    return r,dx,K,D,Sigma,B

def covariance_checks(n=64):
    r,dx,K,D,S,B=system(n)
    A=np.block([[np.zeros((n,n)),np.eye(n)],[-K,-np.eye(n)]])
    forcing=np.zeros_like(A); forcing[n:,n:]=D
    lyap=solve_continuous_lyapunov(A,-forcing)
    error=np.max(abs(lyap-S))/np.max(abs(S))
    injection=.5*np.trace(D)*dx
    heat=np.trace(S[n:,n:])*dx
    subdivision=np.max(abs(sum([D/37]*37)-D))/np.max(abs(D))
    luminosities=np.array([.01,.1,1,10,100])
    amplitudes=np.sqrt(luminosities*S[n//4,n//4])
    exponent=np.polyfit(np.log(luminosities),np.log(amplitudes),1)[0]
    return dict(relative_covariance_error=error,relative_energy_error=abs(injection-heat)/injection,
                injection=injection,damping=heat,subdivision_error=subdivision,
                scaling_exponent=exponent,zero_source_norm=np.max(abs(solve_continuous_lyapunov(A,np.zeros_like(A)))),
                wrong_amplitude_subdivision_ratio=np.sqrt(37),
                covariance_positive_min_eigenvalue=np.linalg.eigvalsh(S).min())

def orbit(n,steps_per_orbit,average,seed):
    r,dx,K,D,S,B=system(n)
    rms=np.sqrt(np.maximum(np.diag(S)[:n],0))
    derivative=np.gradient(rms,r)
    ordinary=1/(1+.2**2)**1.5
    gain=-.1*ordinary/np.interp(1.,r,derivative)
    if gain<=0: raise ValueError('No inward mean response at launch')
    period=2*np.pi/np.sqrt(1.1*ordinary)
    dt=period/steps_per_orbit
    Tav=average*period
    burn=int(np.ceil(max(10.,5*Tav)/dt))
    total=burn+20*steps_per_orbit
    rng=np.random.default_rng(seed)
    a=np.zeros(n); v=np.zeros(n); power=np.zeros(n)
    x=np.array([1.,0.]); vel=np.array([0.,np.sqrt(1.1*ordinary)])
    initialL=x[0]*vel[1]
    E0=0.; supplied=0.; heat=0.; max_budget=0.
    radii=[]; angular=[]; samples=[]; profile_power=np.zeros(n)
    energy=lambda aa,vv: .5*dx*(vv@vv+aa@K@aa)
    def force(xx,aa):
        R=np.linalg.norm(xx)
        mem=gain*np.interp(R,r,np.gradient(np.sqrt(np.maximum(aa,0)),r))
        return (-R/(R*R+.04)**1.5+mem)*xx/R
    decay=np.exp(-dt)
    smooth=-np.expm1(-dt/Tav)
    for k in range(total):
        if k>=burn: vel+=.5*dt*force(x,power); x+=dt*vel
        v-=.5*dt*(K@a)
        a+=.5*dt*v
        before=energy(a,v)
        v*=decay
        heat+=before-energy(a,v)
        before=energy(a,v)
        v+=np.sqrt((1-decay*decay)/2)*(B@rng.standard_normal(n))
        supplied+=energy(a,v)-before
        a+=.5*dt*v
        v-=.5*dt*(K@a)
        power+=smooth*(a*a-power)
        if k>=burn:
            vel+=.5*dt*force(x,power)
            if k%8==0:
                radii.append(np.linalg.norm(x)); angular.append(x[0]*vel[1]-x[1]*vel[0])
                profile_power+=a*a; samples.append(a[n//4])
        if k%128==0: max_budget=max(max_budget,abs(energy(a,v)+heat-supplied))
    rr=np.array(radii)
    spread=np.std(rr)
    drift=abs(np.median(rr)-1)
    Lerr=np.max(abs(np.array(angular)/initialL-1))
    sampled=profile_power/len(samples)
    cov_error=np.linalg.norm(sampled-np.diag(S)[:n])/np.linalg.norm(np.diag(S)[:n])
    return dict(grid=n,steps_per_orbit=steps_per_orbit,average_periods=average,seed=seed,
                gain=gain,burn_periods=burn/steps_per_orbit,radial_std=spread,
                radial_5_50_95=np.quantile(rr,[.05,.5,.95]),median_radius_drift=drift,
                angular_momentum_relative_error=Lerr,
                quiet=spread<.05 and drift<.05 and Lerr<1e-5,
                stationary_sample_diagonal_relative_error=cov_error,
                field_split_ledger_absolute_error=max_budget,
                net_stochastic_input=supplied,damping_heat=heat,final_field_energy=energy(a,v),
                field_ledger_relative_error=max_budget/max(abs(supplied),1e-30),
                sample_radii=rr[::max(1,len(rr)//240)])

def run():
    checks=covariance_checks()
    orbits=[]
    for avg in [.1,1,10]:
        for seed in [1909,1910,1911]:
            row=orbit(64,1024,avg,seed);orbits.append(row)
            print('E3 orbit',avg,seed,row['quiet'],flush=True)
    refined=[orbit(128,2048,avg,1909) for avg in [.1,1,10]]
    gates=dict(covariance=checks['relative_covariance_error']<1e-9,
               energy=checks['relative_energy_error']<1e-9,
               subdivision=checks['subdivision_error']<1e-9,
               luminosity=abs(checks['scaling_exponent']-.5)<.01,
               zero_source=checks['zero_source_norm']==0,
               angular_momentum=all(x['angular_momentum_relative_error']<1e-5 for x in orbits+refined))
    return dict(experiment='E3',checks=checks,gates=gates,numerical_pass=all(gates.values()),
                mechanism_pass=all(x['quiet'] for x in orbits+refined),orbits=orbits,refined=refined,
                observational_status='Unresolved finite photon reservoir, temporal bandwidth, matter backreaction and real-galaxy profiles.')

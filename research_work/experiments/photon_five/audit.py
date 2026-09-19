"""Independent transport moments and isolated optics verification."""
import numpy as np
from scipy.integrate import quad
from common import PM, track
from transport import packets, response
from cluster_lensing_checks import checks
from strain import ray

def run():
    track(PM/'cluster_lensing_checks.py')
    rng=np.random.default_rng(1921)
    count=262144;tau=3.;T=15.
    age=rng.random(count)*T
    alive=rng.random(count)<np.exp(-age/tau)
    packet=T/count
    stored=alive.sum()*packet;decayed=(~alive).sum()*packet
    probability=tau/T*(-np.expm1(-T/tau))
    sigma=T*np.sqrt(probability*(1-probability)/count)
    survival_z=abs(stored-T*probability)/sigma
    moments=[]
    for lam in [0.,3.]:
        pos,diag=packets(32768,lam,0,3,1921)
        square=np.sum(pos*pos,axis=1)
        def msd(a):
            return a*a if lam==0 else 2*(a/lam-(-np.expm1(-lam*a))/lam**2)
        expected=quad(lambda a:msd(a)*np.exp(-a/tau)/(tau*(-np.expm1(-5))),0,T,
                      epsabs=1e-11,epsrel=1e-11)[0]
        se=square.std(ddof=1)/np.sqrt(len(square))
        moments.append(dict(scattering_rate=lam,predicted_msd=expected,packet_msd=square.mean(),
                            standard_error=se,z=abs(square.mean()-expected)/se))
    radius=np.linspace(.1,1.,20);h=.25
    numeric=response(np.zeros((32,3)),radius)
    exact=radius**2/h**4*np.exp(-radius**2/(2*h*h))/(2*np.pi*h*h)
    gaussian_error=np.max(abs(numeric-exact))/np.max(abs(exact))
    optics=checks()
    timestep_rays=[ray(1e-6,1e-8,.5),ray(1e-6,1e-11,.125)]
    ray_change=abs(timestep_rays[0]['deflection']/timestep_rays[1]['deflection']-1)
    gates=dict(finite_packet_budget=abs(T-stored-decayed)<1e-10,
               survival=survival_z<5,flight_moments=all(m['z']<5 for m in moments),
               gaussian_curvature=gaussian_error<1e-4,
               independent_cluster_optics=optics['numerical_verification_passed'],
               ray_timestep_refinement=ray_change<1e-4)
    return dict(experiment='independent audit',gates=gates,numerical_pass=all(gates.values()),
                ray_timestep_relative_change=ray_change,timestep_rays=timestep_rays,
                survival=dict(supplied=T,stored=stored,decayed=decayed,expected_stored=T*probability,z=survival_z),
                moments=moments,gaussian_curvature_relative_error=gaussian_error,optics=optics)

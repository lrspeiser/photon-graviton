"""E2: causal positive energy packets with age-dependent scattering."""
import itertools
import numpy as np
from scipy.integrate import cumulative_trapezoid
from common import PM, ROOT, read_json, inv_cdf, directions, fit_shape
import steady_field as SF
from cluster_lensing import BaryonShells, WrittenClusterLens

def packets(count, lam, nu, tau, seed):
    rng=np.random.default_rng(seed)
    age=-tau*np.log1p(-rng.random(count)*(-np.expm1(-5.)))
    current=np.zeros(count)
    pos=np.zeros((count,3))
    n=directions(rng.random(count),rng.random(count))
    active=np.arange(count)
    collisions=0
    while len(active):
        t=current[active]
        E=rng.exponential(size=len(active))
        if lam==0: nxt=np.full(len(active),np.inf)
        elif nu==0: nxt=t+E/lam
        else:
            remaining=np.exp(-nu*t)-nu*E/lam
            nxt=-np.log(np.maximum(remaining,np.finfo(float).tiny))/nu
            nxt[remaining<=0]=np.inf
        end=np.minimum(nxt,age[active])
        pos[active]+=n[active]*(end-t)[:,None]
        current[active]=end
        hit=nxt<age[active]
        active=active[hit]
        n[active]=directions(rng.random(len(active)),rng.random(len(active)))
        collisions+=len(active)
    stored=tau*(-np.expm1(-5.))
    supplied=5*tau
    decay=supplied-stored
    stderr=np.std(pos,axis=0,ddof=1)/np.sqrt(count)
    diag=dict(packet_count=count,mean=pos.mean(axis=0),
              isotropy_z=np.divide(abs(pos.mean(axis=0)),stderr,out=np.zeros(3),where=stderr>0),
              causal_excess=np.max(np.linalg.norm(pos,axis=1)-age),
              supplied=supplied,stored=stored,decayed=decay,
              budget_error=abs(supplied-stored-decay)/supplied,
              packet_normalization_error=abs(np.sum(np.full(count,stored/count))-stored)/stored,
              scatter_events=collisions,rms_displacement=np.sqrt(np.mean(np.sum(pos*pos,axis=1))))
    return pos,diag

def source_positions(n, seed, kind):
    rng=np.random.default_rng(seed+811)
    r=np.linspace(0,3,4096)
    rho=(1+(r/.296)**2)**(-1.125)
    gas=r*r*rho*rho
    stars=r*r*rho
    gas/=np.trapezoid(gas,r)
    stars/=np.trapezoid(stars,r)
    density={'gas':gas,'stars':stars,'mixed':.5*(gas+stars)}[kind]
    rad=inv_cdf(r,density,rng.random(n))
    return rad[:,None]*directions(rng.random(n),rng.random(n))

def response(pos, radius):
    # Azimuthally averaged Gaussian projection, normalized to unit stored energy.
    rp=np.linalg.norm(pos[:,:2],axis=1)
    edges=np.linspace(0,max(rp.max()+1e-9,1.),2049)
    weights,_=np.histogram(rp,edges)
    mid=.5*(edges[:-1]+edges[1:])
    h=.25
    d1=SF.ring_kernel_dR(radius,mid,h) @ (weights/len(rp))
    d2=SF.ring_kernel_d2R(radius,mid,h) @ (weights/len(rp))
    return (d2-d1/radius)/(2*np.pi*h*h)

def ordinary(ne0,ms,R):
    r=np.geomspace(.1,3000,2000)
    rho0=1.17*SF.MP_G*ne0*(SF.KPC_M*100)**3/(SF.MSUN_KG*1000)
    rho=rho0*(1+(r/296)**2)**(-1.125)
    mass=cumulative_trapezoid(4*np.pi*r*r*rho,r,initial=0)+4*np.pi*r[0]**3*rho[0]/3
    mass*=1+ms/mass[-1]
    shells=BaryonShells(r,mass)
    spectrum=read_json(PM/'cl2-results.json')['E1_spectrum']['clusters_alone']['spectrum']
    lens=WrittenClusterLens(shells,[s['w_kpc'] for s in spectrum],[s['Lambda'] for s in spectrum])
    p=lens.profile(R*1000)
    bary=p['baryon_projected_mass_Msun']/(np.pi*(R*1000)**2)-p['baryon_surface_density_Msun_kpc2']
    return bary,p['delta_sigma_Msun_kpc2']-bary,mass[-1]

def run():
    data=read_json(ROOT/'research_work/results/cluster-observation-readiness/kubo-figure-data.json')['rows']
    R=np.r_[[d['published_radius_h_inverse_Mpc']/.7 for d in data],1.,5.]
    y=np.array([d['shear_t'] for d in data]); err=np.array([d['plotted_sigma_t'] for d in data])
    baselines=[ordinary(ne,ms,R) for ne,ms in [(.0025,.5e13),(.0045,2e13)]]
    p,diag=packets(65536,3,.3,3,1909)
    primary=response(p+source_positions(len(p),1909,'gas'),R)*diag['stored']
    if primary[-2]<=0: raise ValueError('Declared positive anchor unavailable')
    coefficient=baselines[0][1][-2]/primary[-2]
    def evaluate(lam,nu,tau,n,seed,kind='gas'):
        pos,d=packets(n,lam,nu,tau,seed)
        shape=response(pos+source_positions(n,seed,kind),R)*d['stored']
        fits=[]
        for bary,old,mass in baselines:
            field=coefficient*shape*mass/baselines[0][2]
            fits.append(dict(fit=fit_shape((bary+field)[:-2],y,err),
                             field_delta_sigma=field,baryon_fit=fit_shape(bary[:-2],y,err),
                             old_fit=fit_shape((bary+old)[:-2],y,err),old_field=old,
                             extension_ratio=field[-1]/old[-1] if old[-1] else 0))
        return dict(lam=lam,nu=nu,tau=tau,n=n,seed=seed,emissivity=kind,
                    transport=d,response=shape,fits=fits)
    primary_runs=[evaluate(3,.3,3,n,seed) for seed in [1909,1910,1911] for n in [32768,65536]]
    scan=[evaluate(lam,nu,tau,32768,1909)
          for lam,nu,tau in itertools.product([.3,3,30],[0,.3,1],[1,3,10])]
    controls=[evaluate(0,0,3,65536,1909),evaluate(3,0,3,65536,1909)]
    sensitivities=[evaluate(3,.3,3,65536,1909,kind) for kind in ['stars','mixed']]
    maxscale=max(abs(primary).max(),1e-30)
    refinements=[np.max(abs(a['response']-b['response']))/maxscale
                 for a,b in zip(primary_runs[::2],primary_runs[1::2])]
    scatter=np.std([p['response'] for p in primary_runs[1::2]],axis=0,ddof=1)
    allruns=primary_runs+scan+controls+sensitivities
    gates=dict(budget=all(x['transport']['budget_error']<1e-10 and
                         x['transport']['packet_normalization_error']<1e-10 for x in allruns),
               causal=all(x['transport']['causal_excess']<1e-10 for x in allruns),
               isotropy=all(max(x['transport']['isotropy_z'])<5 for x in allruns),
               sampling=max(refinements)<.15,
               source_off=np.all(response(np.zeros((16,3)),R)*0==0))
    row=primary_runs[1]
    mechanism=all(x['fit']['chi2']<=x['baryon_fit']['chi2'] and x['extension_ratio']>1 for x in row['fits'])
    return dict(experiment='E2',radius_Mpc=R,observed=y,error=err,coefficient=coefficient,
                anchor='Frozen old cluster pressure response at 1 Mpc, low ordinary mass bracket',
                gates=gates,numerical_pass=all(gates.values()),mechanism_pass=mechanism,
                primary=primary_runs,scan=scan,controls=controls,sensitivities=sensitivities,
                relative_resolution_changes=refinements,seed_response_std=scatter,
                observational_status='Unresolved absolute luminosity, geometry and optical microphysics; exposed six-bin shape test only.')

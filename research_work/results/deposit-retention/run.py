from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'depleted-capture-feedback'))
from run import Model

def profile(n,k0):
    m=Model(100.,n,8);C=1.;p=6;B=C/k0
    def rhs(s,rho):
        _,wd,_=m.gravity(rho);W=m.wb+wd
        return B*m.absorption(k0*(W/(1+W))**p)[0]/m.V
    cache=HERE/f'profile-k{k0:g}-n{n}.npz'
    if cache.exists() and '--recompute' not in sys.argv:
        rho=np.load(cache)['rho']
    else:
        sol=solve_ivp(rhs,(0,1),np.zeros(n),rtol=1e-8,atol=1e-12)
        assert sol.success and sol.y.min()>=0
        rho=sol.y[:,-1]
    dm,wd,_=m.gravity(rho)
    before=np.cumsum(dm)-dm
    A=before-4*np.pi/3*rho*m.a**3
    def integral_g(lo,hi):
        lo=np.asarray(lo);hi=np.asarray(hi)
        inverse=np.divide(1.,lo,out=np.zeros_like(lo),where=lo>0)-1/hi
        return 1/np.sqrt(1+lo*lo)-1/np.sqrt(1+hi*hi)+A*inverse+2*np.pi/3*rho*(hi*hi-lo*lo)
    shellP=rho*integral_g(m.a,m.b)
    outerP=np.cumsum(shellP[::-1])[::-1]-shellP
    # Integrate independently over each uniform-density shell.
    z,w=roots_legendre(12)
    Kc=0.;Kj=0.;Uself=0.;Uexternal=0.;minbind=float('inf');minep=float('inf')
    for j in range(len(z)):
        r=m.a+(m.b-m.a)*(1+z[j])/2
        weight=(m.b-m.a)*w[j]/2
        enclosed=A+4*np.pi/3*rho*r**3
        mb=r**3/(1+r*r)**1.5
        g=(mb+enclosed)/r**2
        P=outerP+rho*integral_g(r,m.b)
        dv=4*np.pi*r*r*weight
        Kc+=float(np.sum(.5*rho*r*g*dv))
        Kj+=float(np.sum(1.5*P*dv))
        Uself-=float(np.sum(4*np.pi*rho*enclosed*r*weight))
        Uexternal-=float(np.sum(rho/np.sqrt(1+r*r)*dv))
        shell=2*np.pi*rho*(m.b*m.b-m.a*m.a)
        outside=np.cumsum(shell[::-1])[::-1]-shell
        psi=1/np.sqrt(1+r*r)+enclosed/r+2*np.pi*rho*(m.b*m.b-r*r)+outside
        minbind=min(minbind,float(np.min(psi-.5*r*g)))
        # omega_rad^2=(M_total+r M_total')/r^3 for circular test orbits.
        massprime=3*r*r/(1+r*r)**2.5+4*np.pi*rho*r*r
        minep=min(minep,float(np.min((mb+enclosed+r*massprime)/r**3)))
    sample=np.array([.1,.3,1.,3.,10.])
    sampled=np.interp(sample,m.r,rho)
    pressure=np.interp(sample,m.r,outerP+rho*integral_g(m.r,m.b))
    psi_samples=np.interp(sample,m.r,m.wb+wd)
    escape_ratios=3*pressure/(2*sampled*psi_samples)
    deposited_enclosed=np.interp(sample**3,m.edges**3,np.r_[0.,np.cumsum(dm)])
    total_enclosed=deposited_enclosed+sample**3/(1+sample*sample)**1.5
    circular_speed_squared=total_enclosed/sample
    angular_momentum=np.sqrt(sample*total_enclosed)
    selected=m.r>.05;peak=np.where(selected)[0][np.argmax(rho[selected])]
    # A robust outward increase between widely separated registered samples.
    increases=[dict(inner=float(sample[i]),outer=float(sample[j]),density_ratio=float(sampled[j]/sampled[i]))
               for i in range(5) for j in range(i+1,5) if sampled[j]>1.02*sampled[i]]
    result=dict(kappa0=k0,nodes=n,mass=float(dm.sum()),K_circular=Kc,K_Jeans=Kj,
                virial_relative_difference=abs(Kj/Kc-1),U_self=Uself,U_external=Uexternal,
                min_circular_binding=minbind,min_radial_epicyclic_squared=minep,
                density_samples=sampled.tolist(),pressure_samples=pressure.tolist(),
                isotropic_mean_speed_squared_over_escape_squared=escape_ratios.tolist(),
                circular_speed_squared_samples=circular_speed_squared.tolist(),
                circular_specific_angular_momentum_samples=angular_momentum.tolist(),
                outward_density_increases=increases,peak_radius=float(m.r[peak]),
                density_at_point1_over_peak=float(sampled[0]/rho[peak]),
                diagnostics=[float(dm.sum()),Kc,Kj]+sampled.tolist()+pressure.tolist())
    assert result['virial_relative_difference']<1e-6 and minbind>0 and minep>0
    np.savez_compressed(HERE/f'profile-k{k0:g}-n{n}.npz',r=m.r,rho=rho,edges=m.edges,
                        pressure=outerP+rho*integral_g(m.r,m.b))
    return result

results=[]
for k in (.1,10.,1000.):
    a=profile(256,k);b=profile(512,k)
    error=float(max(abs(np.array(a['diagnostics'])/b['diagnostics']-1)))
    entry=dict(coarse=a,refined=b,initial_relative_change=error)
    if error>=.01:
        c=profile(1024,k)
        error=float(max(abs(np.array(b['diagnostics'])/c['diagnostics']-1)))
        entry['refined_again']=c
    entry.update(max_relative_change=error,passes=error<.01)
    results.append(entry)
    print(k,entry.get('refined_again',b),error,flush=True)
(HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf8',newline='\n')
assert all(x['passes'] for x in results),'Retained failure; refine affected diagnostics'

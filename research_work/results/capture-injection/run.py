from pathlib import Path
import sys,json
import numpy as np
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'depleted-capture-feedback'))
from run import Model

def angular(m,k,r,nmu):
    mu,w=roots_legendre(nmu)
    impact2=r*r*(1-mu*mu)
    upper=np.sqrt(np.maximum(m.b[None,:]**2-impact2[:,None],0))
    lower=np.sqrt(np.maximum(m.a[None,:]**2-impact2[:,None],0))
    half=(upper-lower)@k
    along=r*mu
    partial=(np.minimum(upper,abs(along)[:,None])-np.minimum(lower,abs(along)[:,None]))@k
    tau=half+np.sign(along)*partial
    assert tau.min()>-1e-12
    I=np.exp(-tau)
    J=float(w@I/2);F=float(w@(mu*I)/2)
    return dict(J=J,F=F,beta=abs(F)/J)

results=[];checks=[];samples=np.array([.1,.3,1.,3.,10.])
for k0 in (.1,10.,1000.):
    pair=[]
    for n in (512,1024):
        m=Model(100.,n,8)
        rho=np.load(HERE.parent/f'deposit-retention/profile-k{k0:g}-n{n}.npz')['rho']
        _,wd,_=m.gravity(rho);psi=m.wb+wd
        k=k0*(psi/(1+psi))**6
        # B=1/k0, q in mass-equivalent source per normalized exposure.
        absorbed,_,_=m.absorption(k);source=absorbed/k0
        q=source/m.V
        cumulative=np.cumsum(source)-source+source*(m.r**3-m.a**3)/(m.b**3-m.a**3)
        F=-cumulative/(4*np.pi*m.r*m.r)
        J=q/k
        beta=abs(F)/J
        assert np.all(beta>=0) and np.all(beta<1)
        rest=np.sqrt(1-beta*beta);kinetic=1-rest
        # Gamma M/E =1; total radial momentum per E/c is -beta.
        gamma=1/rest
        ledger=max(float(max(abs(rest+kinetic-1))),float(max(abs(gamma*rest*beta-beta))))
        assert ledger<1e-12
        comparisons=[]
        for r in samples:
            a=angular(m,k,r,256);b=angular(m,k,r,512)
            flux_beta=float(np.interp(r,m.r,beta))
            comparisons.append(dict(r=float(r),angular=b,angle_difference=abs(a['beta']-b['beta']),
                                    continuity_beta=flux_beta,continuity_difference=abs(flux_beta-b['beta'])))
        fractions=[]
        for cbar in (300.,1000.,3000.):
            bound=beta*beta<2*psi/cbar**2
            fractions.append(dict(c_over_v0=cbar,bound_energy_fraction=float(source[bound].sum()/source.sum()),
                                  bound_rest_fraction=float((source*rest)[bound].sum()/(source*rest).sum()),
                                  max_weak_field_depth=float(max(psi)/cbar**2)))
        pair.append(dict(kappa0=k0,nodes=n,source_rate=float(source.sum()),
                         rest_energy_fraction=float(source@rest/source.sum()),kinetic_energy_fraction=float(source@kinetic/source.sum()),
                         unresolved_central_source_fraction=float(source[0]/source.sum()),first_shell_outer_radius=float(m.b[0]),
                         max_beta=float(max(beta)),beta_samples=np.interp(samples,m.r,beta).tolist(),
                         samples=comparisons,binding=fractions,local_ledger_error=ledger))
        print(k0,n,pair[-1]['rest_energy_fraction'],fractions,flush=True)
    a,b=pair
    changes=[abs(x['bound_energy_fraction']-y['bound_energy_fraction']) for x,y in zip(a['binding'],b['binding'])]
    beta_change=float(max(abs(np.array(a['beta_samples'])-b['beta_samples'])))
    angular_error=max(x['angle_difference'] for z in pair for x in z['samples'])
    continuity_error=max(x['continuity_difference'] for z in pair for x in z['samples'])
    check=dict(kappa0=k0,bound_fraction_changes=changes,beta_change=beta_change,angular_error=angular_error,
               continuity_error=continuity_error,passes=max(changes)<.01 and beta_change<.005 and angular_error<.005 and continuity_error<.005)
    checks.append(check);results.append(dict(coarse=a,refined=b));print(check,flush=True)
(HERE/'results.json').write_text(json.dumps(dict(results=results,checks=checks),indent=2)+'\n',encoding='utf8',newline='\n')
assert all(x['passes'] for x in checks),'Retain and inspect failed gate'

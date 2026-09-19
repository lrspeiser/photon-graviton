"""E5: reciprocal STF mode, Hamiltonian rays, and static six-lens screen."""
import importlib.util
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_banded
from scipy.optimize import minimize, minimize_scalar
from common import ROOT, PM, read_json, track

D=np.array([2.,-1.,-1.])/np.sqrt(6)
NAMES=['J0037-0942','J1112+0826','J1204+0358','J1402+6321','J1621+3931','J1630+4520']

def homogeneous(weights,coupling=.2):
    weights=np.array(weights)
    def terms(q):
        h=weights*np.exp(-coupling*q*D)
        return h,-coupling*D*h
    def rhs(t,y):
        q,v,heat=y
        h,dh=terms(q)
        return [v,-q-q**3-.4*v-dh.sum(),.4*v*v]
    sol=solve_ivp(rhs,[0,30],[0.,0.,0.],method='DOP853',rtol=1e-11,atol=1e-13,
                  t_eval=np.linspace(0,30,1201))
    if not sol.success:raise RuntimeError(sol.message)
    q,v,heat=sol.y
    photons=(weights[:,None]*np.exp(-coupling*D[:,None]*q)).sum(axis=0)
    total=.5*(q*q+v*v)+q**4/4+photons+heat
    return dict(weights=weights,coupling=coupling,
                relative_energy_error=np.max(abs(total-weights.sum()))/max(weights.sum(),1),
                max_q=np.max(abs(q)),minimum_metric_eigenvalue=np.min(np.exp(-2*coupling*D[:,None]*q)),
                final_photon_energy=photons[-1],final_field_energy=.5*(q[-1]**2+v[-1]**2)+q[-1]**4/4,
                final_heat=heat[-1])

def ray(amplitude,tolerance):
    def rhs(t,y):
        x,z,px,pz=y
        q=amplitude*np.exp(-.5*(x*x+z*z))
        ex,ez=np.exp(-2*q*D[[0,2]])
        H=np.sqrt(ex*px*px+ez*pz*pz)
        factor=(D[0]*ex*px*px+D[2]*ez*pz*pz)/H
        return [ex*px/H,ez*pz/H,-factor*x*q,-factor*z*q]
    sol=solve_ivp(rhs,[0,16],[1.1,-8,0.,1.],method='DOP853',rtol=tolerance,atol=tolerance*1e-3,max_step=.25)
    if not sol.success:raise RuntimeError(sol.message)
    x,z,px,pz=sol.y
    q=amplitude*np.exp(-.5*(x*x+z*z))
    H=np.sqrt(np.exp(-2*D[0]*q)*px*px+np.exp(-2*D[2]*q)*pz*pz)
    return dict(deflection=-np.arctan2(px[-1],pz[-1]),hamiltonian_error=np.max(abs(H-H[0])),
                min_metric=np.min(np.exp(-2*D[:,None]*q)))

def load_lenses():
    folder=ROOT/'research_work/results'
    modelpath=track(folder/'slacs-component-refit/model.py')
    track(folder/'slacs-resolved-fit/model.py')
    spec=importlib.util.spec_from_file_location('pf5_stellar_components',modelpath)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    observed={x['Name']:x for x in read_json(folder/'slacs-resolved-input-audit/results.json')['systems']}
    lights={x['Name']:x for x in read_json(folder/'slacs-light-profile-audit/results.json')['rows']}
    pilot={x['Name']:x for x in read_json(folder/'slacs-motion-lensing-pilot/results.json')['rows']
           if x['model']=='baryons' and x['cutoff_in_a']==20}
    masses=read_json(folder/'lens-photometric-audit/normalization-sensitivity.json')
    geometries=read_json(PM/'cl2-geometry.json')
    lenses=[]
    for name in NAMES:
        obs=observed[name]; light=lights[name]
        geom=geometries[name]['PF1_static_euclidean']
        dl=geom['Dl_kpc']; perarc=dl/module.ARCSEC
        edges=np.r_[obs['inner_arcsec'][0],obs['outer_arcsec']]*perarc
        if not np.allclose(obs['inner_arcsec'][1:],obs['outer_arcsec'][:-1]):raise ValueError('Noncontiguous apertures')
        components=[dict(R=x['R_arcsec']*perarc,n=x['n'],amp=x['amp_at_R'],bn=x['bn']) for x in light['components']]
        a=light['computed_equal_area_half_light_arcsec']*perarc/1.8153
        model=module.ComponentModel(a,edges,obs['psf_fwhm_arcsec']*perarc/2.354820045,
                                    0.,.5,1.,20,components)
        model.forces=model.forces[:1] # ordinary stellar gravity only
        mass={}
        for imf in ['Chabrier','Salpeter']:
            row=next(x for x in masses if x['Name']==name and x['imf']==imf)
            mass[imf]=10**row['published_log10_stellar_mass']*geom['mass_factor']
        y=np.array(obs['vrms_kms']); cov=np.array(obs['covariance_kms_squared'])
        chol=np.linalg.cholesky(cov)
        def objective(par):
            logm,beta=par
            pred=np.sqrt(np.maximum(model.coefficients(beta)[0]*10**logm/1e11,0))
            whiten=np.linalg.solve(chol,pred-y)
            return float(whiten@whiten)
        bounds=[(np.log10(mass['Chabrier']),np.log10(mass['Salpeter'])),(-2,.45)]
        trials=[minimize(objective,[bounds[0][0],b],bounds=bounds,method='L-BFGS-B')
                for b in [-1,0,.4]]
        best=min(trials,key=lambda v:v.fun)
        chab=minimize_scalar(lambda b:objective([bounds[0][0],b]),bounds=(-2,.45),method='bounded')
        mfit=10**best.x[0]
        b=pilot[name]['catalog_SIE_arcsec']*perarc
        theta=pilot[name]['catalog_SIE_arcsec']/module.ARCSEC
        need=theta/geom['Dls_over_Ds']
        def baryon_integrand(t):
            rad=b/np.cos(t)
            return module.G*mfit*model.mass_fraction(rad)/rad
        bary=4/module.C**2*quad(baryon_integrand,0,np.pi/2,epsabs=1e-4,epsrel=1e-8)[0]
        lenses.append(dict(name=name,model=model,a=a,b=b,need=need,baryon=bary,luminosity_proxy=mass['Chabrier']/1e11,
                           geometry=geom,motion=dict(chi2=best.fun,chi2_per_bin=best.fun/len(y),mass=mfit,
                                                    beta=best.x[1],optimizer_success=best.success,
                                                    predicted=np.sqrt(model.coefficients(best.x[1])[0]*mfit/1e11),
                                                    observed=y,chabrier_chi2=chab.fun,population_bounds=mass)))
    return lenses

def tensor_source(lens,nrad=601,nshell=700,angular=64):
    model=lens['model'];a=lens['a']
    r=np.geomspace(max(1e-4,a*1e-4),max(1000,100*a),nrad)
    edges=np.geomspace(model.r[0],r[-1],nshell+1)
    mid=np.sqrt(edges[:-1]*edges[1:])
    enclosed=np.array([model.mass_fraction(v) for v in edges])
    dm=np.diff(np.r_[0,enclosed])
    # Central unresolved stellar cell has its actual integrated luminosity.
    mid=np.r_[edges[0]/2,mid]
    dm*=lens['luminosity_proxy']
    mu,w=np.polynomial.legendre.leggauss(angular)
    source=[]
    for rr in r:
        distance2=rr*rr+mid[:,None]**2-2*rr*mid[:,None]*mu
        cosine2=(rr-mid[:,None]*mu)**2/np.maximum(distance2,1e-30)
        integrand=(1.5*cosine2-.5)/(distance2+.01**2)
        source.append(np.sum(dm*(integrand@(w/2)))/(4*np.pi))
    return r,np.array(source)

def tensor_response(r,source,ell,b,los=1024):
    ds=np.log(r[1]/r[0]);s=ell**2/r**2
    ab=np.zeros((3,len(r)))
    ab[1]=1+s*(2/ds**2+6)
    ab[0,1:]=-s[:-1]*(1/ds**2+1/(2*ds))
    ab[2,:-1]=-s[1:]*(1/ds**2-1/(2*ds))
    ab[1,0]=ab[1,-1]=1;ab[0,1]=0;ab[2,-2]=0
    src=source.copy();src[0]=src[-1]=0
    q=solve_banded((1,1),ab,src)
    # Residual verifies assembled spatial operator.
    product=ab[1]*q
    product[:-1]+=ab[0,1:]*q[1:];product[1:]+=ab[2,:-1]*q[:-1]
    residual=np.max(abs(product-src))/max(abs(src).max(),1e-30)
    interpolator=PchipInterpolator(np.log(r),q,extrapolate=False)
    nodes,weights=np.polynomial.legendre.leggauss(los)
    t=(nodes+1)*np.pi/4
    z=b*np.tan(t);rad=b/np.cos(t)
    active=rad<=r[-1]
    val=np.zeros_like(rad);der=np.zeros_like(rad)
    val[active]=interpolator(np.log(rad[active]))
    der[active]=interpolator.derivative()(np.log(rad[active]))/rad[active]
    dQ=der*b/rad*(z*z/rad**2-1/3)-2*val*z*z*b/rad**4
    bend=-2*np.sum(dQ*b/np.cos(t)**2*weights*np.pi/4)
    return dict(bend=bend,max_abs_q=np.max(abs(q)),operator_residual=residual)

def run():
    homogeneous_runs=[homogeneous([1,0,0]),homogeneous([1/3]*3),homogeneous([0,0,0]),homogeneous([1,0,0],0)]
    rays=[ray(1e-6,1e-8),ray(1e-6,1e-11),ray(0,1e-11)]
    expected=D[2]*1.1*1e-6*np.sqrt(2*np.pi)*np.exp(-1.1**2/2)
    ray_error=abs(rays[0]['deflection']-rays[1]['deflection'])/abs(rays[1]['deflection'])
    analytic_error=abs(rays[1]['deflection']/expected-1)
    lenses=load_lenses()
    sources=[tensor_source(l) for l in lenses]
    need=np.array([l['need'] for l in lenses]);bary=np.array([l['baryon'] for l in lenses])
    scan=[]
    for ell in [.3,1,3,10,30]:
        fields=[tensor_response(r,s,ell,l['b']) for (r,s),l in zip(sources,lenses)]
        unit=np.array([f['bend'] for f in fields])
        w=1/(.03*need[:5])**2
        eta=max(np.sum(unit[:5]*(need-bary)[:5]*w)/np.sum(unit[:5]**2*w),0)
        fractional=(bary+eta*unit)/need-1
        scan.append(dict(length_kpc=ell,coupling=eta,training_score=np.sum((fractional[:5]/.03)**2),
                         fractional_bend_error=fractional,fields=fields,
                         max_abs_tensor=eta*max(f['max_abs_q'] for f in fields)*2/3))
    selected=min(scan,key=lambda x:x['training_score'])
    refined=[]
    for lens,coarse in zip(lenses,selected['fields']):
        r,s=tensor_source(lens,1201,1400,128)
        fine=tensor_response(r,s,selected['length_kpc'],lens['b'],2048)
        fine['relative_bend_change']=abs(fine['bend']-coarse['bend'])/max(abs(fine['bend']),1e-30)
        fine['name']=lens['name'];refined.append(fine)
        print('E5 lens',lens['name'],lens['motion']['chi2_per_bin'],flush=True)
    gates=dict(energy=max(x['relative_energy_error'] for x in homogeneous_runs)<1e-8,
               isotropic=homogeneous_runs[1]['max_q']<1e-10,
               source_off=homogeneous_runs[2]['max_q']<1e-10,
               zero_coupling=homogeneous_runs[3]['max_q']<1e-10 and abs(rays[2]['deflection'])<1e-10,
               metric=all(x['minimum_metric_eigenvalue']>0 for x in homogeneous_runs),
               ray_refinement=ray_error<1e-4,analytic_ray=analytic_error<1e-4,
               spatial_refinement=max(x['relative_bend_change'] for x in refined)<.01,
               motion_optimizer=all(l['motion']['optimizer_success'] for l in lenses))
    screening=dict(bends=np.max(abs(selected['fractional_bend_error']))<=.03,
                   motion=all(l['motion']['chi2_per_bin']<=3 for l in lenses),
                   weak_tensor=selected['max_abs_tensor']<1e-3)
    return dict(experiment='E5',homogeneous=homogeneous_runs,rays=rays,ray_relative_change=ray_error,
                analytic_ray_error=analytic_error,gates=gates,numerical_pass=all(gates.values()),
                screening=screening,mechanism_pass=all(screening.values()),scan=scan,selected=selected,
                refined=refined,lenses=[{k:v for k,v in l.items() if k!='model'} for l in lenses],
                observational_status='Exposed conditional six-lens screen; spatial photon funding and reciprocal matter backreaction unresolved.')

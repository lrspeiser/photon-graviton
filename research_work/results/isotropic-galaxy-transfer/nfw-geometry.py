"""Same-data NFW comparison in alternative and standard benchmark geometries."""
import ast,sys
from pathlib import Path
source=Path(__file__).with_name('lensing.py');tree=ast.parse(source.read_text())
loop=next(n for n in tree.body if isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id=='item')
branchloop=next(n for n in loop.body if isinstance(n,ast.For) and isinstance(n.target,ast.Tuple))
sys.argv.append('--third-retention-optics')
exec(compile(ast.Module(body=tree.body[:tree.body.index(loop)],type_ignores=[]),str(source),'exec'))
sys.argv.remove('--third-retention-optics')
setup=[n for n in loop.body[:loop.body.index(branchloop)] if not (isinstance(n,ast.If) and any(isinstance(v,ast.Continue) for v in n.body))]
ire=next(i for i,n in enumerate(setup) if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Name) and n.targets[0].id=='Re')
def chi(z):return C/70*quad(lambda t:1/np.sqrt(.3*(1+t)**3+.7),0,z,epsabs=1e-9)[0]
def pmass(x):
    if x<1:return np.log(x/2)+np.arccosh(1/x)/np.sqrt(1-x*x)
    if x>1:return np.log(x/2)+np.arccos(1/x)/np.sqrt(x*x-1)
    return 1-np.log(2)
results=[]
for item in data['systems']:
    if item['Name'] not in allowed:continue
    for geometry in ['companion_regular','standard_flat_FLRW']:
        exec(compile(ast.Module(body=setup,type_ignores=[]),str(source),'exec'))
        if geometry=='standard_flat_FLRW':
            dl=chi(zl)/(1+zl)*1000;ratio=1-chi(zl)/chi(zs)
            a=pilot[name]['scale_a_kpc']*dl/(geo[name]['conditional_Dl_Mpc']*1000)
            exec(compile(ast.Module(body=setup[ire:],type_ignores=[]),str(source),'exec'))
        facall=cho_factor(cov);theta=pilot[name]['catalog_SIE_arcsec']/ARCSEC;bE=theta*dl;need=theta/ratio
        starforce=model.forces[0].copy();starbend=4*G*1e11/C**2*quad(lambda t:model.mass_fraction(bE/np.cos(t))/(bE/np.cos(t)),0,np.pi/2,epsabs=1e-9,epsrel=1e-8)[0]
        mb=np.array(cfg['mass_Msun_bounds'])/1e11
        def evaluate(p):
            f,logscale,beta=p;rs=Re*10**logscale;x=r/rs
            halo=1e11*(np.log1p(x)-x/(1+x))
            bend=4*G*1e11/C**2/bE*pmass(bE/rs)
            mstar=(1-f)*need/starbend;amplitude=f*need/bend
            model.forces=np.array([starforce,G*halo/r**2]);cb,ch=model.coefficients(beta)
            pred=np.sqrt(np.maximum(cb*mstar+ch*amplitude,1e-100));e=pred-y
            return float(e@cho_solve(facall,e)),mstar,amplitude,pred,rs,bend
        flo=max(0.,1-mb[1]*starbend/need);fhi=min(1.,1-mb[0]*starbend/need)
        assert flo<fhi
        for kind in ['stars_only_exact_lens','NFW_exact_lens']:
            if kind=='stars_only_exact_lens':
                if not mb[0]<=need/starbend<=mb[1]:continue
                opts=[minimize(lambda p:evaluate([0.,0.,p[0]])[0],[b],bounds=[cfg['constant_beta_bounds']],method='L-BFGS-B') for b in cfg['starts_beta']]
                good=[o for o in opts if o.success];assert good
                best=min(good,key=lambda o:o.fun);p=np.array([0.,0.,best.x[0]])
            else:
                starts=[[flo+(fhi-flo)*f,s,b] for f in [.1,.7] for s in [-1.,0.,1.] for b in [-.5,.3]]
                opts=[minimize(lambda p:evaluate(p)[0],p,bounds=[(flo,fhi),(-2,2),cfg['constant_beta_bounds']],method='L-BFGS-B',options={'ftol':1e-11,'maxiter':400}) for p in starts]
                good=[o for o in opts if o.success and np.isfinite(o.fun)];assert good
                best=min(good,key=lambda o:o.fun);p=best.x
            score,ms,amp,pred,rs,bend=evaluate(p)
            # Analytic NFW projected mass compared with independent force integral.
            check=4*G*1e11/C**2*quad(lambda t:(np.log1p(bE/np.cos(t)/rs)-(bE/np.cos(t)/rs)/(1+bE/np.cos(t)/rs))/(bE/np.cos(t)),0,np.pi/2,epsabs=1e-10,epsrel=1e-8)[0]
            assert abs(check/bend-1)<1e-5
            closure=ratio*(ms*starbend+amp*bend)/theta-1;assert abs(closure)<1e-10
            row=dict(Name=name,geometry=geometry,model=kind,stellar_chi2=score,stellar_mass_Msun=float(ms*1e11),halo_A_Msun=float(amp*1e11),rs_kpc=float(rs),rs_over_Re=float(10**p[1]),halo_deflection_fraction=float(p[0]),beta=float(p[2]),beta_boundary=bool(min(abs(p[2]-np.array(cfg['constant_beta_bounds'])))<1e-4),scale_boundary=bool(abs(abs(p[1])-2)<1e-4),fraction_boundary=bool(min(abs(p[0]-np.array([flo,fhi])))<1e-5),predicted_vrms=pred.tolist(),observed_vrms=y.tolist(),max_marginal_sigma=float(max(abs(pred-y)/np.sqrt(np.diag(cov)))),Dl_Mpc=float(dl/1000),Dls_over_Ds=float(ratio),lens_constraint_residual=float(closure),nfw_projection_relative_error=float(abs(check/bend-1)),optimizer_success_count=len(good))
            results.append(row);print(name,geometry,kind,round(score,3),'f',round(p[0],3),'scale',round(10**p[1],3),flush=True)
(HERE/'nfw-geometry-results.json').write_text(json.dumps(dict(rows=results,standard_geometry=dict(H0=70,Omega_m=.3,Omega_Lambda=.7,role='Alternative benchmark only'),scope='Target-fitted exact lens constraint; no independent prediction'),indent=2,allow_nan=False)+'\n')

"""Inverse compatibility at fixed one-third inventory, not an out-of-sample test."""
import ast
import sys
from pathlib import Path
source=Path(__file__).with_name('lensing.py')
tree=ast.parse(source.read_text())
system_loop=next(n for n in tree.body if isinstance(n,ast.For) and isinstance(n.target,ast.Name) and n.target.id=='item')
branch_loop=next(n for n in system_loop.body if isinstance(n,ast.For) and isinstance(n.target,ast.Tuple))
# Reuse precisely the inspected lens input/geometry setup, without its fit/output driver.
sys.argv.append('--third-retention-optics')
exec(compile(ast.Module(body=tree.body[:tree.body.index(system_loop)],type_ignores=[]),str(source),'exec'))
sys.argv.remove('--third-retention-optics')
scales=np.array([.1,.25,.5,1.,2.,4.,10.]);result=[];verification=[]
prior_run=json.loads((HERE/'third-retention-optics-results.json').read_text())
setup=system_loop.body[:system_loop.body.index(branch_loop)]
# The standalone loop guard includes continue; selection is performed explicitly below.
setup=[n for n in setup if not (isinstance(n,ast.If) and any(isinstance(v,ast.Continue) for v in n.body))]
cut=next(i for i,n in enumerate(branch_loop.body) if isinstance(n,ast.If) and isinstance(n.test,ast.Name) and n.test.id=='REDISTRIBUTION')
profile_setup=branch_loop.body[:cut]
for item in data['systems']:
    if item['Name'] not in allowed:continue
    exec(compile(ast.Module(body=setup,type_ignores=[]),str(source),'exec'))
    for branch,cp in capture.items():
        exec(compile(ast.Module(body=profile_setup,type_ignores=[]),str(source),'exec'))
        # Extend the baseline cumulative mass table to cover every r/s.
        rg=np.geomspace(r[0]/100,r[-1]*100,12001);xx=rg/ac;tt=xx[:,None]*mu;bb2=1+xx[:,None]**2*(1-mu**2);bb=np.sqrt(bb2)
        tau0=cp['k0_per_kpc']*ac*(tt/(2*bb2*(bb2+tt*tt))+(np.arctan(tt/bb)+np.pi/2)/(2*bb**3))
        dens=2*eta*cp['C_Msun_kpc3']*.5*np.sum(np.exp(-np.maximum(tau0,0))*w,axis=1)/(1+xx*xx)**2
        mass0=4*np.pi*(dens[0]*rg[0]**3/3+cumulative_trapezoid(dens*rg*rg,rg,initial=0));M=PchipInterpolator(np.log(rg),mass0,extrapolate=False)
        masses=np.array([M(np.log(r/s)) for s in scales]);assert np.isfinite(masses).all()
        model.forces=np.vstack([model.forces[0],G*masses/r**2])
        old=next(v for v in prior_run['rows'] if v['Name']==name and v['model']==branch)
        oldcoeff=model.coefficients(old['beta']);reproduced=np.sqrt(oldcoeff[0]*old['mass_Msun']/1e11+oldcoeff[4])
        drift=float(max(abs(reproduced-np.array(old['predicted_stellar_vrms']))));assert drift<.05
        verification.append(dict(Name=name,population=cp['population'],baseline_max_velocity_drift_kms=drift))
        beta_bounds=cfg['constant_beta_bounds'];betas=np.linspace(*beta_bounds,101)
        coef=PchipInterpolator(betas,np.array([model.coefficients(b) for b in betas]),axis=0)
        fullfac=cho_factor(cov);angle=pilot[name]['catalog_SIE_arcsec'];thetaE=angle/ARCSEC;bE=thetaE*dl
        starbend=4*G*1e11/C**2*quad(lambda t:model.mass_fraction(bE/np.cos(t))/(bE/np.cos(t)),0,np.pi/2,epsabs=1e-8,epsrel=1e-8)[0]
        def massat(rad,s):
            if rad/s>rg[-1]:return mass0[-1]
            return float(M(np.log(rad/s)))
        bends=np.array([4*G/C**2*quad(lambda t:massat(bE/np.cos(t),s)/(bE/np.cos(t)),0,np.pi/2,epsabs=1e-8,epsrel=1e-8)[0] for s in scales])
        need=thetaE/ratio;massbounds=np.array(cfg['mass_Msun_bounds'])/1e11
        unit=np.eye(7)[3]
        def residual_score(pred):
            e=y-pred;return float(e@cho_solve(fullfac,e))
        def free_obj(p):
            c=coef(p[1]);return residual_score(np.sqrt(np.maximum(c[0]*np.exp(p[0])+unit@c[1:],1e-100)))
        freefits=[minimize(free_obj,[0,b],bounds=[np.log(massbounds),beta_bounds],method='L-BFGS-B') for b in cfg['starts_beta']]
        free=min([o for o in freefits if o.success],key=lambda o:o.fun)
        def fitted_mass(weights):return (need-weights@bends)/starbend
        for kind in ['original_free','original_exact_lens','mixture_exact_lens']:
            exact=False
            def evaluate(p,exact=False):
                weights=unit if kind=='original_exact_lens' else p[:-1]
                mass=fitted_mass(weights);c=model.coefficients(p[-1]) if exact else coef(p[-1]);pred=np.sqrt(np.maximum(c[0]*mass+weights@c[1:],1e-100))
                return residual_score(pred),mass,pred,weights
            if kind=='original_free':
                mass=np.exp(free.x[0]);beta=free.x[1];weights=unit;c=model.coefficients(beta);pred=np.sqrt(c[0]*mass+weights@c[1:]);score=residual_score(pred);success=True
            else:
                if kind=='original_exact_lens':
                    bounds=[beta_bounds];starts=[[b] for b in cfg['starts_beta']];cons=[]
                else:
                    bounds=[(0,1)]*7+[beta_bounds];starts=[np.r_[v,b] for v in [unit,np.ones(7)/7,*np.eye(7)] for b in [-.5,.3]]
                    cons=[{'type':'eq','fun':lambda p:sum(p[:-1])-1}, {'type':'ineq','fun':lambda p:fitted_mass(p[:-1])-massbounds[0]}, {'type':'ineq','fun':lambda p:massbounds[1]-fitted_mass(p[:-1])}]
                opts=[minimize(lambda p:evaluate(p)[0],p,method='SLSQP',bounds=bounds,constraints=cons,options={'ftol':1e-9,'maxiter':400}) for p in starts]
                good=[o for o in opts if o.success and massbounds[0]<=evaluate(o.x)[1]<=massbounds[1]]
                if not good:
                    result.append(dict(Name=name,population=cp['population'],fit=kind,feasible=False));continue
                opt=min(good,key=lambda o:o.fun)
                refined=minimize(lambda p:evaluate(p,True)[0],opt.x,method='SLSQP',bounds=bounds,constraints=cons,options={'ftol':1e-9,'maxiter':200})
                chosen=refined if refined.success and evaluate(refined.x,True)[0]<=evaluate(opt.x,True)[0] else opt
                score,mass,pred,weights=evaluate(chosen.x,True);beta=chosen.x[-1];success=bool(chosen.success)
            closure=ratio*(starbend*mass+weights@bends)/thetaE-1
            if kind!='original_free':assert abs(closure)<1e-9
            assert abs(sum(weights)-1)<1e-6 and min(weights)>-1e-7
            changes=[dict(radius_Re=z,baseline_mass=float(M(np.log(z*Re))),new_mass=float(sum(wt*M(np.log(z*Re/s)) for wt,s in zip(weights,scales)))) for z in [.5,1,2,5]]
            row=dict(Name=name,population=cp['population'],fit=kind,feasible=True,stellar_chi2=score,stellar_mass_Msun=float(mass*1e11),beta=float(beta),beta_boundary=bool(min(abs(beta-np.array(beta_bounds)))<1e-4),mass_boundary=bool(min(abs(mass-massbounds))<1e-6),weights=weights.tolist(),predicted_vrms=pred.tolist(),observed_vrms=y.tolist(),max_marginal_sigma=float(max(abs(pred-y)/np.sqrt(np.diag(cov)))),lens_equation_fractional_residual=float(closure),catalog_angle_arcsec=angle,cumulative_mass_changes=changes,eta=float(eta),optimizer_success=success)
            result.append(row);print(name,cp['population'],kind,round(score,3),round(row['max_marginal_sigma'],3),flush=True)
out=dict(scope='Target-fitted profile compatibility; no predictive success claim',scales=scales.tolist(),rows=result,verification=verification,capture_input_sha256=hashlib.sha256((HERE/'third-radiation-retention-results.json').read_bytes()).hexdigest(),optical_input_sha256=prior_run['optical_input_sha256'],photometric_input_sha256=prior_run['photometric_input_sha256'])
(HERE/'lens-profile-compatibility-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')

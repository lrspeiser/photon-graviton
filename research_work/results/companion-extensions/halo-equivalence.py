"""Frozen halo-to-companion substitution; execute setup, never the optimizer."""
from pathlib import Path
import json,hashlib,ast
import numpy as np
from scipy.integrate import cumulative_trapezoid,quad
from scipy.interpolate import PchipInterpolator
from scipy.special import eval_legendre
P=Path(__file__).resolve().parent
source=P/'free-nfw.py'
prefix=source.read_text().split('    seeds=[];reproduction=[]')[0]
assert 'trials=[minimize' not in prefix
test=r'''
    frozen=next(v for v in json.loads((P/'free-nfw-results.json').read_text())['rows'] if v['Name']==name)
    rs=frozen['nfw_rs_kpc'];amp=frozen['nfw_mass_amplitude_Msun'];mass=frozen['stellar_mass_Msun']
    h=frozen['h'];b0=frozen['beta0'];bi=frozen['beta_infinity'];ra=Re*frozen['orbit_transition_radius_Re']
    weight=h*meanH/(1+h*meanH)
    beta=b0+(bi-b0)*r*r/(r*r+ra*ra)
    factor=np.exp(2*b0*np.log(r/model.a)+(bi-b0)*np.log((r*r+ra*ra)/(model.a*model.a+ra*ra)))
    starforce=mass/1e11*((1-weight)*model.forces[0]+weight*model.forces[1])
    records=[]
    for count in [8193,16385]:
        # Integrate density in log radius; analytic inner cusp correction only.
        xx=np.geomspace(1e-12,max(1e8,float(max(r/rs))*2),count)
        integrand=xx**2/(1+xx)**2
        enclosed=.5*xx[0]**2+cumulative_trapezoid(integrand,np.log(xx),initial=0)
        lut=PchipInterpolator(np.log(xx),enclosed)
        num=lut(np.log(r/rs));ana=nfw_fraction(r/rs)
        force=starforce+G*amp*num/r**2
        pressure=-cumulative_trapezoid((nu*force*factor)[::-1],r[::-1],initial=0)[::-1]/factor
        pred=np.sqrt(np.trapezoid(r*r*pressure*(model.W-beta*model.T),r,axis=1)/model.den)
        err=y-pred;score=float(err@cho_solve(full_fac,err))
        records.append(dict(grid=count,max_halo_force_relative_error=float(max(abs(num/ana-1))),max_stellar_speed_difference_kms=float(max(abs(pred-frozen['predicted_stellar_vrms']))),chi2=score,predicted_stellar_vrms=pred.tolist()))
    # Independent density-to-cylinder projection, without using enclosed-mass bending.
    b=impact/rs
    inner=quad(lambda t:t/(1+t)**2,0,b,epsabs=1e-12,epsrel=1e-10)[0]
    def shell(t):
        v=(b/t)**2
        fraction=v/(1+np.sqrt(max(0,1-v)))
        return t/(1+t)**2*fraction
    outer=quad(shell,b,np.inf,epsabs=1e-12,epsrel=1e-10,limit=300)[0]
    density_bending=4*G*amp*(inner+outer)/(C**2*impact)
    analytic_bending=amp*unit_deflection(np.log(rs/Re))
    row=dict(Name=name,nfw_rs_kpc=rs,amplitude_Msun=amp,baseline_chi2=frozen['all_motion_chi2'],refinement=records,lensing_relative_difference=float(abs(density_bending/analytic_bending-1)) if analytic_bending else float(abs(density_bending)))
    assert records[-1]['max_stellar_speed_difference_kms']<.01,row
    assert row['lensing_relative_difference']<1e-6,row
    out['rows'].append(row)
    print(name,records[-1]['max_stellar_speed_difference_kms'],row['lensing_relative_difference'],flush=True)
'''
env={'__file__':str(source),'__name__':'halo_equivalence_replay'}
exec(compile(prefix+test,str(source)+' [frozen diagnostic]','exec'),env)
out=env['out'];out['scope']='Frozen six spherical NFW halo fits replaced by identical companion density. No capture law, parameter fit, new data or physical origin inference.'
# Exercise the existing nonspherical companion field solver with the same NFW density.
tree=ast.parse((P/'directional-capture.py').read_text())
node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='field')
checks=[]
for nr in [768,1536]:
 x=np.geomspace(1e-8,1e5,nr);u,uw=np.polynomial.legendre.leggauss(48);ls=np.arange(0,25,2)
 basis=np.array([eval_legendre(l,u)*(2*l+1)*uw/2 for l in ls]);p0=np.array([eval_legendre(l,0) for l in ls])
 ns=dict(np=np,nr=nr,x=x,ls=ls,basis=basis,p0=p0,PchipInterpolator=PchipInterpolator)
 exec(compile(ast.Module(body=[node],type_ignores=[]),'production-field','exec'),ns)
 rho=np.broadcast_to((1/(x*(1+x)**2))[:,None],(nr,len(u)))
 fn,_=ns['field'](rho);rr=np.geomspace(.001,1000,150)
 exact=(np.log1p(rr)-rr/(1+rr))/rr
 error=float(max(abs(fn(np.log(rr))/exact-1)));checks.append(dict(radial_grid=nr,max_relative_circular_speed_squared_error=error))
assert checks[-1]['max_relative_circular_speed_squared_error']<.001
out['production_axisymmetric_field_checks']=checks
out['summary']=dict(baseline_chi2=sum(r['baseline_chi2'] for r in out['rows']),substituted_chi2=sum(r['refinement'][-1]['chi2'] for r in out['rows']),max_speed_difference_kms=max(r['refinement'][-1]['max_stellar_speed_difference_kms'] for r in out['rows']),max_lensing_relative_difference=max(r['lensing_relative_difference'] for r in out['rows']))
paths=[source,P/'free-nfw-results.json',P/'directional-capture.py',Path(__file__)]
out['source_sha256']={str(p.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
(P/'halo-equivalence-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(out['summary']);print(checks)

"""Fixed common-parameter, fixed-photometric-mass training comparison."""
from pathlib import Path
import importlib.util
import hashlib,json,argparse
import numpy as np
from scipy.integrate import cumulative_trapezoid,quad
from scipy.optimize import brentq
from scipy.stats import ncx2
from numpy.polynomial.legendre import leggauss

H=Path(__file__).resolve().parent;R=H.parent
parser=argparse.ArgumentParser();parser.add_argument('--verify-resolution',action='store_true');parser.add_argument('--grid-index',type=int,choices=range(9));parser.add_argument('--role',choices=['training','validation','test'],default='training');parser.add_argument('--orbit-audit',action='store_true');args=parser.parse_args()
if args.role!='training' and args.grid_index is not None:raise ValueError('Parameter grid is training-only.')
if args.role=='test' and args.orbit_audit:raise ValueError('No exploratory orbital alternatives on the final test sample.')
spec=importlib.util.spec_from_file_location('wave_equilibrium',R/'self-consistent-wave/run.py')
wave=importlib.util.module_from_spec(spec);spec.loader.exec_module(wave)
files=[H/'protocol.json',R/'lens-photometric-audit/normalization-sensitivity-updated-profile.json',
       R/'lensing-data-readiness/lens-observations-and-image-models.json',R/'lensing-data-readiness/conditional-geometry.json']
protocol,masses,observed,geometry=[json.loads(f.read_text()) for f in files]
OUT=H
if args.grid_index is not None:
    gpfile=H/'grid-protocol.json';gp=json.loads(gpfile.read_text());files.append(gpfile)
    protocol['field_mass_eV_c2']=gp['field_mass_values_eV_c2'][args.grid_index//3]
    protocol['source_to_stellar_mass_ratio']=gp['source_to_stellar_mass_values'][args.grid_index%3]
    protocol['choice']='Declared shared-parameter training grid; not tuned per galaxy'
    protocol['grid_protocol']=gp
    OUT=H/'parameter-grid'/f'cell-{args.grid_index}'
    OUT.mkdir(parents=True,exist_ok=True)
if args.role in ['validation','test']:
    OUT=R/f'wave-lens-{args.role}'
    vpfile=OUT/'protocol.json';vp=json.loads(vpfile.read_text())
    selectedfile=H/'parameter-grid/selected-calibration.json';selected=json.loads(selectedfile.read_text())
    for key in ['field_mass_eV_c2','source_to_stellar_mass_ratio']:protocol[key]=selected[key]
    protocol.update({'role':args.role,'choice':'Frozen training-selected common parameters','evaluation_protocol':vp})
    files[1]=OUT/'mass-inputs.json';masses=json.loads(files[1].read_text())
    files.extend([vpfile,selectedfile])
    if args.role=='test':
        for file,digest in vp['frozen_inputs'].items():
            assert hashlib.sha256((R/file).read_bytes()).hexdigest()==digest, f'Frozen input changed: {file}'
        masscalfile=R/'stellar-mass-degeneracy/frozen-calibration.json';files.append(masscalfile)
        shared_mass_factor=json.loads(masscalfile.read_text())['lambda_mass']
if args.orbit_audit:
    opfile=R/'stellar-orbit-degeneracy/protocol.json'
    protocol['orbital_audit']=json.loads(opfile.read_text());files.append(opfile)
masses=[r for r in masses if r['model']=='baryons' and r['imf']=='Salpeter' and r['propagation_branch']=='energy_loss_and_event_stretch']
expected=32 if args.role=='training' else (7 if args.role=='validation' else json.loads((OUT/'preparation.json').read_text())['eligible'])
assert len(masses)==expected and all(r['role']==args.role for r in masses)
if args.verify_resolution:
    old=json.loads((OUT/'results.json').read_text())['equilibrium_checks']
    names={old[0]['Name'],min(old,key=lambda r:r['eta'])['Name'],max(old,key=lambda r:r['eta'])['Name']}
    masses=[r for r in masses if r['Name'] in names]
observed={r['Name']:r for r in observed if r['role']==args.role}
geometry={r['Name']:r for r in geometry if r['role']==args.role}
G=4.30091727003628e-6;KPC=3.085677581491367e19;RAD=np.pi/(180*3600);C=299792.458
particle=protocol['field_mass_eV_c2']*1.7826619216279e-36
source_fraction=protocol['source_to_stellar_mass_ratio']
t,w=leggauss(256 if args.verify_resolution else 128)
x=np.geomspace(1e-6,1e5,24000 if args.verify_resolution else 12000);j=1/(x*(1+x)**3)
rows=[];checks=[];orbitrows=[]
for mrow in masses:
    name=mrow['Name'];obs=observed[name];geo=geometry[name]
    dl=geo['conditional_Dl_Mpc']*1000;ratio=geo['conditional_Dls_over_Ds']
    re=mrow['new_I_Re_arcsec'];a=dl*re*RAD/1.8153;M=10**mrow['conditional_log10_stellar_mass']
    eta=1.054571817e-34**2/(particle**2*(G*M*1e6*KPC)*(a*KPC))
    prev=None
    for f in np.geomspace(.03,source_fraction,50):
        prev=wave.solve(eta,float(f),previous=prev)
    sol=wave.solve(eta,source_fraction,L=90,tol=2e-9,previous=prev,eps=5e-4)
    check=wave.diagnostics(sol,eta,source_fraction,90)
    assert check['virial_relative_residual']<2e-5 and check['normalization_relative_error']<1e-6
    check['Name']=name
    check['max_dimensionless_potential_depth_over_c2']=float((G*M/a)/C**2*(1-sol.sol(sol.x[0])[2]))
    checks.append(check)
    def enclosed(z):
        arr=np.asarray(z);clipped=np.clip(arr,sol.x[0],90)
        value=sol.sol(clipped)[3]
        return np.where(arr<sol.x[0],sol.sol(sol.x[0])[3]*(arr/sol.x[0])**3,
                        np.where(arr>90,source_fraction,value))
    ap=1.5/re*1.8153;s=protocol['seeing_fwhm_arcsec']/2.354820045/re*1.8153
    top=np.arcsin(np.minimum(1,(ap+8*s)/x))
    theta=(t[None,:]+1)*top[:,None]/2
    prob=ncx2.cdf((ap/s)**2,2,(x[:,None]*np.sin(theta)/s)**2)
    W=np.sum(prob*np.sin(theta)*w[None,:],axis=1)*top/2
    kernel=cumulative_trapezoid(x*x*W,x,initial=0)
    kernels={0.:kernel};angular_weights={0.:W}
    if args.orbit_audit:
        for beta in [-.3,.3]:
            Wbeta=np.sum(prob*np.sin(theta)*(1-beta*np.sin(theta)**2)*w[None,:],axis=1)*top/2
            angular_weights[beta]=Wbeta
            kernels[beta]=x**(2*beta)*cumulative_trapezoid(x**(2-2*beta)*Wbeta,x,initial=0)
    denom=np.trapezoid(j*x*x*W,x)
    for model in ['baryons','stationary_wave']+(['shared_stellar_mass'] if args.role=='test' else []):
        multiplier=shared_mass_factor if model=='shared_stellar_mass' else 1.
        gc=enclosed(x)/x**2 if model=='stationary_wave' else np.zeros_like(x)
        sigma=np.sqrt(multiplier*G*M/a*np.trapezoid(j*(1/(1+x)**2+gc)*kernel,x)/denom)
        def bend(angle):
            b=dl*angle*RAD/a
            def integrand(t):
                rr=b/np.cos(t)
                mm=rr**2/(1+rr)**2
                if model=='stationary_wave': mm+=float(enclosed(rr))
                return mm*np.cos(t)/b
            return multiplier*4*G*M/(a*C*C)*quad(integrand,0,np.pi/2,epsabs=1e-9,epsrel=2e-8)[0]*ratio/RAD
        angle=np.exp(brentq(lambda log:bend(np.exp(log))/np.exp(log)-1,-24,4))
        b=dl*angle*RAD/a
        def cap(r):return (b/r)**2/(1+np.sqrt(1-(b/r)**2))
        projected=b*b/(1+b)**2+quad(lambda r:2*r/(1+r)**3*cap(r),b,np.inf,epsabs=1e-10)[0]
        if model=='stationary_wave':
            projected+=float(enclosed(b))
            if b<90:projected+=quad(lambda r:sol.sol(r)[0]**2*cap(r),b,90,epsabs=1e-10)[0]
        projected_bend=multiplier*4*G*M/(a*C*C)*projected/b*ratio/RAD
        lens_error=abs(projected_bend/angle-1)
        assert lens_error<1e-6
        rows.append({'Name':name,'role':args.role,'model':model,'eta':eta,'stellar_mass_Msun':M*multiplier,
                     'sigma_pred_km_s':float(sigma),'sigma_observed_km_s':obs['sigma'],
                     'sigma_error_km_s':obs['e_sigma'],'theta_pred_arcsec':float(angle),
                     'theta_SIE_arcsec':obs['bSIE'],'projected_mass_lensing_identity_relative_error':lens_error})
        if args.orbit_audit:
            for beta,kbeta in sorted(kernels.items()):
                entry=dict(rows[-1]);entry['beta']=beta
                entry['sigma_pred_km_s']=float(np.sqrt(G*M/a*np.trapezoid(j*(1/(1+x)**2+gc)*kbeta,x)/denom))
                pressure_integrand=j*(1/(1+x)**2+gc)*x**(2*beta)
                pressure=-cumulative_trapezoid(pressure_integrand[::-1],x[::-1],initial=0)[::-1]*x**(-2*beta)
                variance_direct=G*M/a*np.trapezoid(x*x*angular_weights[beta]*pressure,x)/denom
                entry['independent_pressure_variance_relative_error']=abs(variance_direct/entry['sigma_pred_km_s']**2-1)
                assert entry['independent_pressure_variance_relative_error']<1e-4
                orbitrows.append(entry)
    print(name,flush=True)
summary={'classification':f'Fixed shared-parameter conditional {args.role} comparison; not capture theory',
         'systems':len(masses),'protocol':protocol,'scores':{},'equilibrium_checks':checks,
         'input_sha256':{str(f.relative_to(R)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}}
for model in ['baryons','stationary_wave']+(['shared_stellar_mass'] if args.role=='test' else []):
    ss=[r for r in rows if r['model']==model]
    ds=np.array([r['sigma_pred_km_s']-r['sigma_observed_km_s'] for r in ss])
    da=np.array([r['theta_pred_arcsec']-r['theta_SIE_arcsec'] for r in ss])
    summary['scores'][model]={'sigma_rmse_km_s':float(np.sqrt(np.mean(ds**2))),
        'sigma_mean_residual_km_s':float(ds.mean()),'theta_rmse_arcsec':float(np.sqrt(np.mean(da**2))),
        'theta_mean_residual_arcsec':float(da.mean()),
        'median_sigma_pred_over_observed':float(np.median([r['sigma_pred_km_s']/r['sigma_observed_km_s'] for r in ss])),
        'median_theta_pred_over_SIE':float(np.median([r['theta_pred_arcsec']/r['theta_SIE_arcsec'] for r in ss]))}
for filename,data in [('predictions.json',rows),('results.json',summary)]:
    if args.orbit_audit:
        if filename=='predictions.json':data=orbitrows
        filename='orbital-'+filename
    if args.verify_resolution:filename=filename.replace('.json','-refined.json')
    (OUT/filename).write_text(json.dumps(data,indent=2)+'\n',newline='\n')
print(json.dumps(summary['scores'],indent=2))

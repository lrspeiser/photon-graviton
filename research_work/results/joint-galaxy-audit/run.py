"""Joint available-data audit; see protocol.json for scope and predeclared choices."""
import csv, hashlib, io, json, platform, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares, minimize, minimize_scalar
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]
C=299792.458; KPC=3.085677581491367e19; G=4.30091727003628e-6
INPUTS=['redshift_paper/all_164_groups.csv','temporal_candidate_audit/data/Rotmod_LTG.zip',
        'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt','temporal_candidate_audit/data/sparc_frozen.json',
        'research_work/results/milky-way-capture/inputs.json','research_work/results/electromagnetic-audit/results.json']

def save(name,obj):
    HERE.joinpath(name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rmse(x): return float(np.sqrt(np.mean(np.asarray(x)**2)))
def load(path): return json.loads((ROOT/path).read_text(encoding='utf-8'))

# A. Same old group splits, but train-only estimation of a history modification.
with (ROOT/INPUTS[0]).open(encoding='utf-8',newline='') as f: rr=list(csv.DictReader(f))
D=np.array([float(r['catalog_distance_mpc']) for r in rr]); x=D/100
y=np.array([float(r['observed_cmb_cz_kms']) for r in rr])
split=np.array([r['split'] for r in rr]); train=split=='train'
sd=D*np.log(10)/5*np.array([float(r['sbf_modulus_error_mag']) for r in rr])
alpha_fixed=.0002488993286382367
fit=minimize_scalar(lambda a:np.mean((C*np.expm1(a*D[train])-y[train])**2),bounds=(0,.001),method='bounded',options={'xatol':1e-15})
alpha=float(fit.x); acceleration_scale=(C*1000)**2*alpha/(1000*KPC)
# Use endpoint slopes a and a+2*b*xmax, both positive, as fit coordinates.
xmax=float(x.max())
def history(par):
    a,end=par; b=(end-a)/(2*xmax)
    return C*np.expm1(a*x+b*x*x)
history_fit=least_squares(lambda p:(history(p)[train]-y[train])/300,[100*alpha,100*alpha],bounds=([0,0],[.2,.2]),xtol=1e-12,ftol=1e-12,gtol=1e-10)
red={}; redrows=[]
preds={'fixed_alpha':C*np.expm1(alpha_fixed*D),'refit_alpha':C*np.expm1(alpha*D),'quadratic_history':history(history_fit.x)}
for name,yp in preds.items():
    if name=='quadratic_history':
        a,end=history_fit.x; b=(end-a)/(2*xmax); derivative=C/100*np.exp(a*x+b*x*x)*(a+2*b*x)
    else: derivative=C*(alpha_fixed if name=='fixed_alpha' else alpha)*np.exp((alpha_fixed if name=='fixed_alpha' else alpha)*D)
    sensitivity_sigma=np.sqrt(300**2+(derivative*sd)**2)
    red[name]={'scores':{}}
    for s in ['train','validation','test']:
        ii=split==s; e=yp[ii]-y[ii]
        red[name]['scores'][s]=dict(n=int(ii.sum()),RMSE_kms=rmse(e),bias_kms=float(np.mean(e)),
            within_300kms=int(np.sum(abs(e)<=300)),sensitivity_rms_standardized=rmse(e/sensitivity_sigma[ii]))
    for i,r in enumerate(rr): redrows.append(dict(model=name,group_pgc=r['group_pgc'],split=r['split'],D_Mpc=float(D[i]),observed_z=float(y[i]/C),predicted_z=float(yp[i]/C),residual_kms=float(yp[i]-y[i])))
red['parameters']={'alpha_per_Mpc':alpha,'alpha_fixed_per_Mpc':alpha_fixed,'a_star_m_s2':acceleration_scale,
    'history_a':float(history_fit.x[0]),'history_b':float((history_fit.x[1]-history_fit.x[0])/(2*xmax)),
    'history_optimizer_success':bool(history_fit.success),'max_distance_Mpc':float(D.max())}
save('redshift-predictions.json',redrows)
print('Redshift complete',flush=True)

# B. Refit the existing empirical extra-acceleration family on SPARC train only.
base=ROOT/'temporal_candidate_audit/data'; meta={}
for line in (base/'SPARC_Lelli2016c.mrt').read_text().splitlines():
    f=line.split()
    if len(f)!=19: continue
    try: meta[f[0]]=dict(method=int(f[4]),inc=float(f[5]),rd=float(f[11]),quality=int(f[17]))
    except ValueError: pass
galaxies={}
with zipfile.ZipFile(base/'Rotmod_LTG.zip') as zf:
    for fn in sorted(zf.namelist()):
        name=fn.replace('_rotmod.dat',''); m=meta[name]
        if m['quality']>2 or m['inc']<30 or m['rd']<=0: continue
        arr=np.atleast_2d(np.loadtxt(io.BytesIO(zf.read(fn))))
        vb2=arr[:,3]*abs(arr[:,3])+.5*arr[:,4]*abs(arr[:,4])+.7*arr[:,5]*abs(arr[:,5])
        good=np.isfinite(arr).all(axis=1)&(arr[:,0]>0)&(arr[:,1]>0)&(arr[:,2]>0)&(vb2>0)
        if good.sum()<5: continue
        arr=arr[good]; vb2=vb2[good]
        galaxies[name]=dict(name=name,**m,R=arr[:,0],v=arr[:,1],err=arr[:,2],vb2=vb2,gb=vb2*1e6/(arr[:,0]*KPC))
frozen=load('temporal_candidate_audit/data/sparc_frozen.json'); splits=frozen['split']
assert len(galaxies)==149 and sum(len(d['v']) for d in galaxies.values())==3150
def extra(gb,par): return 10**par[0]*acceleration_scale*(gb/acceleration_scale)**par[1]
def velocity(d,par): return np.sqrt(d['vb2']+extra(d['gb'],par)*d['R']*KPC/1e6)
def loss(par): return np.mean([np.mean(np.log10(velocity(galaxies[n],par)/galaxies[n]['v'])**2) for n in splits['train']])
opts=[minimize(loss,p,method='L-BFGS-B',bounds=[(-4,2),(0,1)],options={'ftol':1e-14,'gtol':1e-9}) for p in [[-.6,.5],[-1,0],[-1,1]]]
opt=min(opts,key=lambda o:o.fun); par=opt.x
sparc={'parameters':{'A':float(10**par[0]),'p':float(par[1]),'a_star_m_s2':acceleration_scale,'optimizer_success':bool(opt.success)},'scores':{}}; sprows=[]
for model in ['baryons','power']:
    sparc['scores'][model]={}
    for s,names in splits.items():
        errors=[]; logs=[]
        for n in names:
            d=galaxies[n]; yp=np.sqrt(d['vb2']) if model=='baryons' else velocity(d,par)
            errors.append(np.mean((yp-d['v'])**2)); logs.append(np.mean(np.log10(yp/d['v'])**2))
            for i in range(len(yp)): sprows.append(dict(model=model,galaxy=n,split=s,R_kpc=float(d['R'][i]),observed_kms=float(d['v'][i]),predicted_kms=float(yp[i]),gas_included=True))
        sparc['scores'][model][s]=dict(n_galaxies=len(names),n_rows=sum(len(galaxies[n]['v']) for n in names),galaxy_weighted_RMSE_kms=float(np.sqrt(np.mean(errors))),galaxy_weighted_log_RMSE=float(np.sqrt(np.mean(logs))))
save('galaxy-rotation-predictions.json',sprows)
print('149-galaxy transfer calibration complete',flush=True)

# C. Gas-only numerical potentials with a stored, hashed coefficient/version audit.
mw=load('research_work/results/milky-way-capture/inputs.json'); er=mw['eilers']['rows']; br=mw['bovy']['rows']
Robs=np.array([r['R_kpc'] for r in er]); Vobs=np.array([r['vc_kms'] for r in er]); Rz=np.array([r['R_kpc'] for r in br]); Zobs=np.array([r['Kz_over_2piG_Msun_pc2'] for r in br]); Ze=np.array([r['Kz_error'] for r in br])
from gas import make_gas,force,PARAMS
import galpy
rad=np.unique(np.r_[np.geomspace(.1,120,150),Robs,np.sqrt(Rz**2+1.1**2)])
cachefile=HERE/'gas-cache.json'
gas_source_hash=hashlib.sha256((HERE/'gas.py').read_bytes()).hexdigest()
if cachefile.exists():
    gc=json.loads(cachefile.read_text()); assert gc['gas_source_sha256']==gas_source_hash and gc['galpy_version']==galpy.__version__
    assert np.allclose(rad,gc['R_grid'])
else:
    with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
        gas30=make_gas(30); gas40=make_gas(40)
    gf=np.array([force(gas40,r,0)[0] for r in rad]); gz=np.array([force(gas40,r,1.1)[1] for r in Rz])
    checkR=np.array([5.27,8.19,12.25,20.27,24.82])
    lower=np.array([force(gas30,r,z) for r in checkR for z in [0,1.1]])
    higher=np.array([force(gas40,r,z) for r in checkR for z in [0,1.1]])
    absdiff=abs(lower-higher); rel=absdiff/np.maximum(abs(higher),1.)
    assert np.isfinite(gf).all() and np.isfinite(gz).all()
    gc=dict(gas_source_sha256=gas_source_hash,galpy_version=galpy.__version__,R_grid=rad.tolist(),radial_force=gf.tolist(),vertical_force=gz.tolist(),parameters=PARAMS,
            check_max_relative=float(rel.max()),check_max_absolute_force=float(absdiff.max()),check_R=checkR.tolist(),order30=lower.tolist(),order40=higher.tolist())
    save('gas-cache.json',gc)
gasR=PchipInterpolator(rad,np.array(gc['radial_force'])); gasZ=np.array(gc['vertical_force'])
print('Gas force calculation complete',gc['check_max_relative'],flush=True)

def stars(R,z,variant):
    # Known Miyamoto-Nagai formula; a=0 gives the Plummer bulge.
    components={'I':[(460*2.32e7,0,.3),(1700*2.32e7,5.3,.25),(1700*2.32e7,2.6,.8)],
                'II':[(1600*2.32e7,4.8,.25),(1700*2.32e7,2.,.8)]}[variant]
    fr=np.zeros_like(np.asarray(R,dtype=float)); fz=fr.copy()
    for mass,a,b in components:
        B=np.sqrt(z*z+b*b); den=(R*R+(a+B)**2)**1.5
        fr+=G*mass*R/den; fz+=G*mass*(a+B)*z/(B*den)
    return fr,fz

# Exact all-angle optical depth for beta=beta0*rs^2/(rs^2+r^2).
# Constant beta is a separate control, not an extremely large fitted rs.
def capture_profile(theta,kind,nr=900,nmu=96):
    mass,beta0=10**theta[0],10**theta[1]; rs=10**theta[2] if kind=='well' else None
    radius=np.r_[0,np.geomspace(.0001,100,nr-1)]
    mu,w=leggauss(nmu); r=radius[:,None]; impact2=r*r*(1-mu*mu); line=r*mu
    entry=np.sqrt(np.maximum(0,100**2-impact2))
    if kind=='uniform':
        tau=beta0*(line+entry); beta=np.full(len(radius),beta0)
    else:
        B=np.sqrt(rs*rs+impact2)
        tau=beta0*rs*rs/B*(np.arctan(line/B)+np.arctan(entry/B))
        beta=beta0/(1+(radius/rs)**2)
    J=.5*np.sum(w*np.exp(-np.maximum(0,tau)),axis=1)
    shape=beta*J; integral=cumulative_trapezoid(4*np.pi*radius**2*shape,radius,initial=0)
    assert integral[-1]>0
    enclosed=mass*integral/integral[-1]
    return radius,enclosed,shape/integral[-1]*mass

mwresults={}; mwrows=[]; profiles=[]
cal=Robs<=15; sigma=np.sqrt(np.array([.5*(r['err_minus_kms']+r['err_plus_kms']) for r in er])**2+(.03*Vobs)**2)
for variant in ['I','II']:
    baryR=stars(Robs,0,variant)[0]+gasR(Robs); baryZ=stars(Rz,1.1,variant)[1]+gasZ
    assert np.all(baryR>0)
    # Conservative spherical completion of the empirical midplane extra force.
    gb=baryR*1e6/KPC; gp=extra(gb,par)*KPC/1e6
    rz=np.sqrt(Rz**2+1.1**2); gbz=(stars(rz,0,variant)[0]+gasR(rz))*1e6/KPC
    transferZ=extra(gbz,par)*KPC/1e6*1.1/rz
    predictions={'baryons':(np.sqrt(Robs*baryR),baryZ),'sparc_transfer':(np.sqrt(Robs*(baryR+gp)),baryZ+transferZ)}
    mwresults[variant]={'capture_fits':{},'scores':{}}
    for kind in ['uniform','well']:
        def pred(theta,nr=900,nmu=96):
            rg,mg,rho=capture_profile(theta,kind,nr,nmu); interp=PchipInterpolator(rg,mg)
            return np.sqrt(Robs*baryR+G*interp(Robs)/Robs),baryZ+G*interp(rz)*1.1/rz**3
        lo=[9,-5]+([-.3] if kind=='well' else []); hi=[13,1]+([1.7] if kind=='well' else [])
        starts=[[11.5,-3],[11.5,-1]] if kind=='uniform' else [[11.5,-3,.5],[11.5,-1,.5],[11.5,-2,1.3]]
        fits=[least_squares(lambda p:(pred(p)[0][cal]-Vobs[cal])/sigma[cal],p,bounds=(lo,hi),max_nfev=150,ftol=1e-9,xtol=1e-9,gtol=1e-8) for p in starts]
        f=min(fits,key=lambda q:np.sum(q.fun*q.fun)); lowpred=pred(f.x); highpred=pred(f.x,1800,192)
        numerical_delta=max(float(np.max(abs(lowpred[0]-highpred[0]))),float(np.max(abs(lowpred[1]-highpred[1]))/(2*np.pi*G*1e6)))
        predictions['capture_'+kind]=highpred
        mwresults[variant]['capture_fits'][kind]=dict(log_parameters=f.x.tolist(),mass_100kpc_Msun=float(10**f.x[0]),beta0_per_kpc=float(10**f.x[1]),rs_kpc=float(10**f.x[2]) if kind=='well' else None,
            success=bool(f.success),calibration_score=float(np.sum(f.fun*f.fun)),at_bounds=[bool(abs(p-l)<.001 or abs(p-h)<.001) for p,l,h in zip(f.x,lo,hi)],
            max_resolution_change_kms_or_surface_units=numerical_delta)
        rg,mg,rho=capture_profile(f.x,kind,1800,192)
        profiles.append(dict(baryons=variant,kind=kind,R_kpc=rg.tolist(),mass_enclosed_Msun=mg.tolist(),rho_Msun_kpc3=rho.tolist()))
    for name,(vr,kz) in predictions.items():
        kzsurface=kz/(2*np.pi*G*1e6); err=vr-Vobs
        mwresults[variant]['scores'][name]=dict(radial_calibration_RMSE_kms=rmse(err[cal]),radial_outer_RMSE_kms=rmse(err[~cal]),radial_all_RMSE_kms=rmse(err),
            vertical_RMSE_surface_units=rmse(kzsurface-Zobs),vertical_diagonal_score=float(np.sum(((kzsurface-Zobs)/Ze)**2)),vertical_n=43)
        for i in range(len(Robs)): mwrows.append(dict(baryons=variant,model=name,observable='vc_kms',split='inner' if cal[i] else 'outer',R_kpc=float(Robs[i]),z_kpc=0,observed=float(Vobs[i]),predicted=float(vr[i])))
        for i in range(len(Rz)): mwrows.append(dict(baryons=variant,model=name,observable='abs_Kz_over_2piG_Msun_pc2',split='vertical_provisional',R_kpc=float(Rz[i]),z_kpc=1.1,observed=float(Zobs[i]),predicted=float(kzsurface[i])))
    print('Milky Way model',variant,'complete',flush=True)

checks={}
# Weak uniform absorption must approach uniform density and M proportional to r^3.
rg,mg,rho=capture_profile([11,-8],'uniform',1800,192)
checks['transparent_mass_half_ratio']=float(np.interp(50,rg,mg)/mg[-1])
assert abs(checks['transparent_mass_half_ratio']-.125)<2e-5
# The total incoming-minus-outgoing energy equals integrated capture, in a sphere.
# Independent ray chords for the variable beta case vs volume integration.
beta=.03; rs=7.; radius=np.r_[0,np.linspace(.0001,100,1800)]; mu,w=leggauss(192)
b=100*(mu+1)/2; wb=w*50
B=np.sqrt(rs*rs+b*b); chordtau=2*beta*rs*rs/B*np.arctan(np.sqrt(100**2-b*b)/B)
captured_ray=np.sum(2*b/100**2*(-np.expm1(-chordtau))*wb)
r=radius[:,None]; impact2=r*r*(1-mu*mu); BB=np.sqrt(rs*rs+impact2)
tau=beta*rs*rs/BB*(np.arctan(r*mu/BB)+np.arctan(np.sqrt(np.maximum(0,100**2-impact2))/BB))
J=.5*np.sum(w*np.exp(-tau),axis=1); shape=beta/(1+(radius/rs)**2)*J
captured_volume=4/100**2*np.trapezoid(radius**2*shape,radius)
checks['capture_energy_relative_error']=float(abs(captured_volume/captured_ray-1))
assert checks['capture_energy_relative_error']<5e-5
checks['redshift_energy_fraction_error']=float(np.max(abs(np.exp(-alpha*D)+(-np.expm1(-alpha*D))-1)))
assert checks['redshift_energy_fraction_error']<1e-14
checks['alpha_degeneracy_prediction_error']=float(np.max(abs(extra(np.array([1e-12,1e-10,1e-8]),par)-10**(par[0]-(1-par[1])*np.log10(2))*(2*acceleration_scale)*(np.array([1e-12,1e-10,1e-8])/(2*acceleration_scale))**par[1])))
assert checks['alpha_degeneracy_prediction_error']<1e-22
save('milky-way-predictions.json',mwrows); save('capture-profiles.json',profiles)
summary=dict(status='Conditional available-data audit, not a completed theory',redshift=red,sparc=sparc,milky_way=mwresults,checks=checks,
             gas_numerics=gc['check_max_relative'],software=dict(python=platform.python_version(),numpy=np.__version__,galpy=galpy.__version__))
save('results.json',summary)
save('input-manifest.json',{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS})
print(json.dumps(summary,indent=2),flush=True)

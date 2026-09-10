"""Training-only spherical Jeans/lensing pilot; not a completed theory."""
from pathlib import Path
import json,hashlib,argparse
import numpy as np
from scipy.integrate import quad,cumulative_trapezoid
from scipy.optimize import brentq
from scipy.stats import ncx2
from numpy.polynomial.legendre import leggauss
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
parser=argparse.ArgumentParser();parser.add_argument('--refined',action='store_true');parser.add_argument('--updated-profile',action='store_true');args=parser.parse_args()
datafile=H.parent/'lensing-data-readiness/lens-observations-and-image-models.json'
geomfile=H.parent/'lensing-data-readiness/conditional-geometry.json'
fitfile=H.parent/'joint-galaxy-audit/results.json'
data=json.loads(datafile.read_text());geo={r['Name']:r for r in json.loads(geomfile.read_text())}
data=[r for r in data if r['role']=='training' and r['Mph']=='E' and r['spectroscopic_dispersion_available']]
assert len(data)==33
profilefile=H.parent/'lens-photometric-audit/training-photometry.json'
if args.updated_profile:
    profiles={r['SDSS']:r for r in json.loads(profilefile.read_text())}
    for r in data:
        assert profiles[r['Name']]['Re(I)'] is not None
        r['Re']=profiles[r['Name']]['Re(I)']
pars=json.loads(fitfile.read_text())['sparc']['parameters'];G=4.30091727003628e-6;C=299792.458;KPC=3.085677581491367e19;RAD=np.pi/(180*3600)
A,p,astar=pars['A'],pars['p'],pars['a_star_m_s2']*KPC/1e6;MREF=1e11
x=np.geomspace(1e-6,1e5,12000 if args.refined else 6000);j=1/(x*(1+x)**3)
nq=128 if args.refined else 64;t,w=leggauss(nq)
def weights(ap,seeing):
    if seeing==0:
        u=np.minimum(1,(ap/x)**2);W=u/(1+np.sqrt(1-u))
        K=x**3/3*np.where(u==1,1,-np.expm1(1.5*np.log1p(-np.minimum(u,1-1e-16))))
    else:
        s=seeing/2.354820045
        top=np.arcsin(np.minimum(1,(ap+8*s)/x))
        theta=(t[None,:]+1)*top[:,None]/2
        P=ncx2.cdf((ap/s)**2,2,(x[:,None]*np.sin(theta)/s)**2)
        W=np.sum(P*np.sin(theta)*w[None,:],axis=1)*top/2
        K=cumulative_trapezoid(x*x*W,x,initial=0)
    return W,K
rows=[]
for r in data:
    dl=geo[r['Name']]['conditional_Dl_Mpc']*1000;ratio=geo[r['Name']]['conditional_Dls_over_Ds']
    a=dl*r['Re']*RAD/1.8153;ap=1.5/r['Re']*1.8153
    gb=G*MREF/a**2/(1+x)**2
    for fwhm in [0.,1.5]:
        W,K=weights(ap,fwhm/r['Re']*1.8153)
        denom=np.trapezoid(j*x*x*W,x)
        kb=a*np.trapezoid(j*gb*K,x)/denom
        for cutRe in [5.,20.,100.]:
            xt=cutRe*1.8153;gt=A*astar*((G*MREF/a**2/(1+xt)**2)/astar)**p
            gc=np.where(x<=xt,A*astar*(gb/astar)**p,gt*(xt/x)**2)
            kc=a*np.trapezoid(j*gc*K,x)/denom
            for model in ['baryons','empirical_companion']:
                if model=='baryons' and cutRe!=20:continue
                def mass(sig):return np.exp(brentq(lambda l:kb*np.exp(l)+(kc*np.exp(p*l) if model=='empirical_companion' else 0)-sig**2,-20,20))
                def deflection(theta,m):
                    b=dl*theta*RAD;u=b/a
                    # Fixed quadrature split at companion boundary; pointlike
                    # exterior tail analytic, baryon Hernquist extends to infinity.
                    bary=quad(lambda t:G*MREF*m/a**2/(1+u/np.cos(t))**2*b/np.cos(t),0,np.pi/2,epsabs=1e-7,epsrel=1e-9)[0]
                    comp=0.
                    if model=='empirical_companion':
                        if u<xt:
                            tm=np.arccos(u/xt)
                            comp=quad(lambda t:A*astar*((G*MREF*m/a**2/(1+u/np.cos(t))**2)/astar)**p*b/np.cos(t),0,tm,epsabs=1e-7,epsrel=1e-9)[0]
                            comp+=gt*m**p*(a*xt)**2/b*(u/xt)**2/(1+np.sin(tm))
                        else:comp=gt*m**p*(a*xt)**2/b
                    return 4*(bary+comp)/C**2*ratio/RAD
                predictions=[];masses=[]
                for sig in [r['sigma']-r['e_sigma'],r['sigma'],r['sigma']+r['e_sigma']]:
                    m=mass(sig);masses.append(m*MREF)
                    # Solve in log-angle; no Einstein angle used in mass inference.
                    rt=brentq(lambda l:deflection(np.exp(l),m)/np.exp(l)-1,-14,8)
                    predictions.append(float(np.exp(rt)))
                rows.append(dict(Name=r['Name'],role='training',model=model,seeing_fwhm_arcsec=fwhm,cutoff_over_Re=cutRe,
                    mass_from_sigma_Msun=masses[1],theta_pred_arcsec=predictions[1],theta_sigma_minus_arcsec=predictions[0],theta_sigma_plus_arcsec=predictions[2],
                    theta_SIE_arcsec=r['bSIE'],theta_LTM_arcsec=r['bLTM'],stellar_axis_ratio=r['b/a'],residual_arcsec=predictions[1]-r['bSIE']))
    print('computed',r['Name'],flush=True)
summary={'training_systems':len(data),'parameters':pars,'radial_grid_nodes':len(x),'angle_quadrature_nodes':nq,'scores':{},
    'holdout_scores_opened':False,'mass_fitted_to':'Spectroscopic aperture dispersion only, separately under each potential; not an independently measured stellar mass.',
    'assumptions':['Static Euclidean geometry using archived redshift distances','Spherical isotropic Hernquist light and ordinary-mass profile, Re=1.8153a','Uniform mass-to-light ratio; no separately modeled gas or central black hole','3-arcsecond diameter spectroscopic aperture','Seeing zero and Gaussian FWHM 1.5 arcsec are sensitivity cases, not per-object measurements','Equal temporal/spatial companion metric potentials','Source truncations 5,20,100 Re are prescribed sensitivity cases, not capture predictions','SIE intermediate-axis Einstein angle is an approximate circular target; no imaging likelihood or model covariance'],
    'source_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [datafile,geomfile,fitfile]}}
for model,seeing,cut in sorted(set((r['model'],r['seeing_fwhm_arcsec'],r['cutoff_over_Re']) for r in rows)):
    ss=[r for r in rows if (r['model'],r['seeing_fwhm_arcsec'],r['cutoff_over_Re'])==(model,seeing,cut)]
    res=np.array([r['residual_arcsec'] for r in ss]);rat=np.array([r['theta_pred_arcsec']/r['theta_SIE_arcsec'] for r in ss])
    summary['scores'][f'{model}/seeing{seeing}/cut{cut}']={'n':len(ss),'rmse_arcsec':float(np.sqrt(np.mean(res**2))),'median_predicted_over_SIE':float(np.median(rat)),'mean_residual_arcsec':float(res.mean())}
summary['updated_I_profile_used']=args.updated_profile
if args.updated_profile:summary['source_sha256'][str(profilefile.name)]=hashlib.sha256(profilefile.read_bytes()).hexdigest()
suffix=('-refined' if args.refined else '')+('-updated-profile' if args.updated_profile else '')
for name,obj in [('predictions',rows),('results',summary)]:
    (H/(name+suffix+'.json')).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(summary['scores'],indent=2))

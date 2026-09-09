"""Stationary thermal absorption/emission plus photon-energy conversion."""
from pathlib import Path
import hashlib,json,os
import numpy as np
from scipy.constants import h,k,c
from scipy.integrate import quad
from scipy.optimize import minimize_scalar,least_squares
from scipy.special import roots_laguerre

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'thermalizer'
Y,W=roots_laguerre(128)

def shape(x):
    x=np.asarray(x)
    return x**3*np.exp(-x)/(-np.expm1(-x))

def stationary_shape(x,a,nodes=Y,weights=W):
    return np.dot(shape(np.asarray(x)[...,None]*np.exp(nodes/a)),weights)

def planck(nu,t):
    return 2*(k*t)**3/(h*h*c*c)*shape(h*nu/(k*t))/1e-20

def spectrum(nu,t,a):
    return 2*(k*t)**3/(h*h*c*c)*stationary_shape(h*nu/(k*t),a)/1e-20

def main():
    protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    path=ROOT/protocol['input'];data=np.loadtxt(path)
    assert data.shape==(43,5) and np.all(data[:,3]>0)
    # Check the characteristic integral independently with adaptive quadrature.
    numerical=[]
    n0=quad(lambda x:float(shape(x))/x,0,100,epsabs=1e-10)[0]
    u0=quad(lambda x:float(shape(x)),0,100,epsabs=1e-10)[0]
    y2,w2=roots_laguerre(256)
    for a in protocol['absorption_to_conversion_ratios']:
        errs=[];diffs=[]
        for x in np.geomspace(.01,30,40):
            val=float(stationary_shape(x,a))
            adaptive=quad(lambda y:np.exp(-y)*float(shape(x*np.exp(y/a))),0,100,epsabs=1e-15,epsrel=1e-10)[0]
            errs.append(abs(val/adaptive-1))
            diffs.append(abs(val/float(stationary_shape(x,a,y2,w2))-1))
        n=quad(lambda x:float(stationary_shape(x,a))/x,0,100,epsabs=1e-9)[0]
        u=quad(lambda x:float(stationary_shape(x,a)),0,100,epsabs=1e-9)[0]
        assert max(errs)<1e-8 and max(diffs)<1e-8
        assert abs(n/n0-1)<1e-8 and abs(u/u0-a/(a+1))<1e-8
        # Steady energy-space equation: H nu u' + kappa(B-u)=0.
        residuals=[]
        for x in np.geomspace(.01,30,60):
            dx=x*1e-4;f=lambda z:float(stationary_shape(z,a))
            deriv=(-f(x+2*dx)+8*f(x+dx)-8*f(x-dx)+f(x-2*dx))/(12*dx)
            residuals.append(abs(x*deriv+a*(float(shape(x))-f(x)))/f(x))
        assert max(residuals)<1e-6
        numerical.append({'a':a,'number_ratio_to_material_Planck':n/n0,
                          'energy_ratio_to_material_Planck':u/u0,
                          'max_adaptive_relative_difference':max(errs),
                          'max_256_node_relative_difference':max(diffs),
                          'max_scaled_transport_residual':max(residuals)})
    nu=data[:,0]*100*c;obs=data[:,2];sig=data[:,3];gal=data[:,4]
    ref=planck(nu,protocol['reference_temperature_K']);wg=gal/sig
    lo,hi=protocol['temperature_bounds_K'];rows=[]
    for a in [None]+protocol['absorption_to_conversion_ratios']:
        model=lambda t:1000*((planck(nu,t) if a is None else spectrum(nu,t,a))-ref)
        def profile(t,full=False):
            r=(obs-model(t))/sig;g=float(np.dot(r,wg)/np.dot(wg,wg));r-=g*wg
            return (float(np.dot(r,r)),g) if full else float(np.dot(r,r))
        fit=minimize_scalar(profile,bounds=(lo,hi),method='bounded',options={'xatol':1e-12})
        chi,g=profile(fit.x,True)
        repeat=least_squares(lambda p:(model(p[0])+p[1]*gal-obs)/sig,[fit.x,g],
                             bounds=([lo,-np.inf],[hi,np.inf]),x_scale='jac',
                             ftol=1e-12,xtol=1e-12,gtol=1e-12)
        assert fit.success and repeat.success and lo<repeat.x[0]<hi
        chi2=float(np.dot(repeat.fun,repeat.fun));assert abs(chi2-chi)<1e-5*max(1,chi)
        rows.append({'absorption_to_conversion_ratio':a,'case':'unchanged Planck' if a is None else 'thermalizer with conversion',
                     'material_temperature_K':float(repeat.x[0]),'Galaxy_template_coefficient':float(repeat.x[1]),
                     'diagonal_chi_square':chi2,'weighted_residuals':repeat.fun.tolist()})
    result={'scope':protocol['scope'],'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'protocol_sha256':hashlib.sha256((HERE/'protocol.json').read_bytes()).hexdigest(),
            'numerical_checks':numerical,'fixed_ratio_profiles':rows,
            'energy_ledger':'U=a/(a+1) U_P(T_d), N=N_P(T_d). Required external material heating is H_c U; gross emission is kappa U_P(T_d), absorption is kappa U. No heating source or gravitating receiver is supplied.',
            'conditional_conclusion':'Finite gray thermalization produces a non-Planck stationary spectrum. Rapid thermalization approaches Planck while requiring external energy to replace conversion losses.',
            'unresolved':['Actual material opacity and temperature distribution','Source-funded heating, conversion vertex and receiving sectors','FIRAS covariance, calibration and foreground likelihood','Background angular structure, optical transparency and finite history']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'thermalizer-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'numerical_checks':numerical,'profiles':[{k:v for k,v in r.items() if k!='weighted_residuals'} for r in rows]},indent=2))

if __name__=='__main__':main()

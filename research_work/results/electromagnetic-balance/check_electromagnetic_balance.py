"""Constant electric, magnetic and parity-odd scalar couplings in a finite band."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'electromagnetic-balance'

def coefficients(v):
    qc=np.sqrt(1-v*v)
    w=lambda q:q*np.sqrt(v*v+q*q)
    A=lambda q:q*q-w(q)**2
    i0=qc*qc/2
    i1=quad(lambda q:q*w(q),0,qc)[0]
    i2=quad(lambda q:q*w(q)**2,0,qc)[0]
    ia=quad(lambda q:q*A(q),0,qc)[0]
    iaw=quad(lambda q:q*A(q)*w(q),0,qc)[0]
    ia2=quad(lambda q:q*A(q)**2,0,qc)[0]
    assert np.isclose(i2-ia/2,qc**4*(1+v*v)/8,rtol=1e-12)
    def plus(e):
        return (e*i0-2*i1+(i2-ia/2)/e+iaw/(2*e**2)+ia2/(16*e**3))/(16*np.pi)
    def minus(e):return ia2/(256*np.pi*e**3)
    return qc,plus,minus

def best_flatness(v,low_ratio,R=2.):
    qc,plus,minus=coefficients(v);e0=low_ratio*qc
    p0=plus(e0)
    p=lambda y:plus(e0*y)/p0
    # Positive monotone relative basis and convex total: minimax endpoints balance.
    t=(p(R)-1)/(1-R**-3)
    rate=lambda y:p(y)+t*y**-3
    fit=minimize_scalar(rate,bounds=(1,R),method='bounded',options={'xatol':1e-13})
    assert fit.success and 1<fit.x<R
    top=rate(1);bottom=rate(fit.x)
    assert abs(top-rate(R))<1e-12*top
    grid=np.linspace(1,R,4001);values=rate(grid)
    assert np.max(values)<=top*(1+1e-12) and np.min(values)>=bottom*(1-1e-12)
    err=(top-bottom)/(top+bottom)
    # Independent bounded minimization of mixing ratio using dense energy grid.
    def objective(logt):
        vals=p(grid)+np.exp(logt)*grid**-3
        return (max(vals)-min(vals))/(max(vals)+min(vals))
    alt=minimize_scalar(objective,bounds=(np.log(t)-3,np.log(t)+3),method='bounded',options={'xatol':1e-12})
    assert alt.success and abs(alt.fun-err)<1e-7
    return {'v':v,'low_energy_over_band_edge':low_ratio,'energy_interval_ratio':R,
            'best_low_endpoint_minus_to_plus_contribution':t,
            'required_coupling_weight_B_over_A':t*p0/minus(e0),
            'minimum_rate_energy_over_low_endpoint':float(fit.x),
            'max_rate_over_min_rate':top/bottom,
            'best_uniform_relative_error_with_free_normalization':err,
            'independent_grid_optimization_error':float(alt.fun)}

def main():
    rng=np.random.default_rng(30720);errors=[]
    for _ in range(300):
        mu=rng.uniform(-1,1);s=np.sqrt(1-mu*mu)
        ge,gb,go=rng.normal(size=3)
        ni=np.array([0.,0.,1.]);no=np.array([s,0.,mu])
        incoming=[np.array([1.,0.,0.]),np.array([0.,1.,0.])]
        outgoing=[np.array([mu,0.,-s]),np.array([0.,1.,0.])]
        total=0.
        for ei in incoming:
            for eo in outgoing:
                bi=np.cross(ni,ei);bo=np.cross(no,eo)
                amplitude=ge*np.dot(ei,eo)-gb*np.dot(bi,bo)+go*(np.dot(ei,bo)+np.dot(eo,bi))
                total+=amplitude**2/2
        decomposed=((ge-gb)**2*(1+mu)**2+((ge+gb)**2+4*go**2)*(1-mu)**2)/4
        errors.append(abs(total-decomposed)/max(1,total))
    assert max(errors)<1e-12
    integration=[]
    previous=json.loads((ROOT/'research_work/results/dispersive-companion/dispersive-companion-results.json').read_text(encoding='utf-8'))
    for v in [.1,.5,.9]:
        qc,plus,minus=coefficients(v)
        for ratio in [1,2,10,100]:
            e=ratio*qc
            def integrand(q,sgn):
                w=q*np.sqrt(v*v+q*q);ep=e-w
                mu=1-(q*q-w*w)/(2*e*ep)
                return q*ep*ep*(1+sgn*mu)**2/(64*np.pi*e)
            for sgn,fn in [(1,plus),(-1,minus)]:
                direct=quad(lambda q:integrand(q,sgn),0,qc,epsabs=1e-24,epsrel=1e-10)[0]
                integration.append(abs(direct/fn(e)-1))
            assert abs(minus(2*e)/minus(e)-1/8)<1e-14
        for row in previous['scenarios']:
            if row['v']==v and row['E_over_band_edge']>=1:
                e=row['E_over_M']
                assert np.isclose(plus(e)+minus(e),row['alpha_over_g_squared_M_cubed'],rtol=1e-10)
    assert max(integration)<1e-7
    rows=[best_flatness(v,r) for v in [.1,.5,.9] for r in [1,10,100,10000]]
    t=8/7;y=(3*t)**.25;top=1+t;bottom=y+t/y**3
    asym=(top-bottom)/(top+bottom)
    assert all(abs(r['best_uniform_relative_error_with_free_normalization']-asym)<1e-4 for r in rows if r['low_energy_over_band_edge']==10000)
    result={'scope':'Unpolarized spontaneous tree-level finite-band conversion with one scalar and constant local isotropic phi E^2, phi B^2 and phi E dot B couplings. No observational data fitted.',
            'polarization_factor':'[A(1+mu)^2+B(1-mu)^2]/4; A=(g_E-g_B)^2>=0, B=(g_E+g_B)^2+4g_O^2>=0',
            'rate_structure':'alpha(E)=A F_plus(E)+B D/E^3 for E>=Qc; F_plus has positive E leading coefficient and positive inverse-power coefficients, plus a constant term.',
            'octave_flatness_profiles':rows,
            'high_energy_octave_limit':{'best_low_endpoint_contribution_ratio':t,'minimum_energy_over_low_endpoint':y,'max_rate_over_min_rate':top/bottom,'best_uniform_relative_error':asym},
            'checks':{'polarization_decomposition_max_scaled_error':max(errors),'direct_vs_moment_integral_max_relative_error':max(integration)},
            'verdict':'Opposite color slopes permit a locally flat point but not a nonzero constant fractional loss over any full-band energy interval. At high energies the best unpolarized octave-wide uniform flatness error approaches the recorded positive floor.',
            'limits':['Finite-band and UV hierarchy assumptions inherited from dispersive-companion','Polarization-resolved transport, populated baths, frequency-dependent couplings and extra operators are not covered','No matter/clock, storage, gravity or cosmological completion']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'electromagnetic-balance-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'high_energy_octave_limit':result['high_energy_octave_limit'],'profiles':rows},indent=2))

if __name__=='__main__':main()

"""Independent mathematical checks; no astronomical inference."""
from pathlib import Path
import hashlib, json, os
import numpy as np
from scipy.integrate import quad
from scipy.special import log_ndtr, ndtr
from scipy.stats import multivariate_normal, norm
from likelihood import log_flux_likelihood, width_shape

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT/'research_work/generated'))/'timing-likelihood'


def dense(y, error, shape, am=1., ast=.5, bst=.1):
    covariance=np.diag(error**2)+ast**2*np.outer(shape,shape)+bst**2*np.ones((len(y),len(y)))
    solved=np.linalg.solve(covariance,y-am*shape)
    postmean=am+ast**2*shape@solved
    postvar=ast**2-ast**4*shape@np.linalg.solve(covariance,shape)
    return float(multivariate_normal.logpdf(y,mean=am*shape,cov=covariance)
                 +log_ndtr(postmean/np.sqrt(postvar))-log_ndtr(am/ast))


def numerical(y,error,shape,am=1.,ast=.5,bst=.1):
    # Integrate baseline independently; then integrate the positive amplitude.
    cov=np.diag(error**2)+bst**2*np.ones((len(y),len(y)))
    def integrand(a):
        return np.exp(multivariate_normal.logpdf(y,mean=a*shape,cov=cov)
                      +norm.logpdf(a,loc=am,scale=ast))/ndtr(am/ast)
    value,err=quad(integrand,0,np.inf,epsabs=1e-12,epsrel=1e-11,limit=300)
    assert value>0 and err < max(1e-11, value*1e-8)
    return float(np.log(value))


def main():
    protocol=json.loads((HERE/'protocol.json').read_text())
    tol=protocol['numerical_tolerances']
    rng=np.random.default_rng(20260909)
    dense_errors=[];quad_errors=[]
    for n in [1,2,5,12,30]:
        for level in [.1,.5,2.]:
            shape=rng.uniform(0,1,n);error=rng.uniform(.1,.5,n)
            y=level*shape+rng.normal(0,error)
            actual=log_flux_likelihood(y,error,shape)
            dense_errors.append(abs(actual-dense(y,error,shape)))
            if n<=5:quad_errors.append(abs(actual-numerical(y,error,shape)))
    # Negative measured flux probes the truncation term, not only A~1 cases.
    y=np.array([-.2,-.1]);error=np.array([.3,.4]);shape=np.array([1.,.8])
    quad_errors.append(abs(log_flux_likelihood(y,error,shape)-numerical(y,error,shape)))
    assert max(dense_errors)<tol['log_likelihood_absolute']
    assert max(quad_errors)<tol['log_likelihood_absolute']
    width_errors=[]
    for w in [2.,20.,160.]:
        for r in [.15,.3,.7]:
            for p,q in [(1.25,3.),(2.,2.),(3.,1.25)]:
                points=np.array([5-r*w,5.,5+(1-r)*w])
                width_errors.append(float(np.max(abs(width_shape(points,w,5.,r,p,q)-[.5,1.,.5]))))
    assert max(width_errors)<tol['half_maximum_absolute']
    integral,integral_error=quad(lambda f:np.exp(log_flux_likelihood([f],[.3],[.8])),
                                 -np.inf,np.inf,epsabs=1e-10,epsrel=1e-10)
    assert abs(integral-1)<tol['probability_integral_absolute']
    assert log_flux_likelihood([],[],[])==0
    collinear=log_flux_likelihood([1.,1.,1.],[.1,.1,.1],[1.,1.,1.])
    assert np.isfinite(collinear)
    invalid=0
    for y,e,s in [([1],[0],[1]),([1],[.1,.2],[1]),([float('nan')],[1],[1])]:
        try:log_flux_likelihood(y,e,s)
        except ValueError:invalid+=1
    assert invalid==3
    result={'scope':'Mathematical likelihood component only. No duration/population recovery or observational gate is tested.',
            'dense_covariance_comparisons':len(dense_errors),'max_dense_log_error':max(dense_errors),
            'independent_positive_amplitude_integrals':len(quad_errors),'max_quadrature_log_error':max(quad_errors),
            'half_maximum_cases':len(width_errors),'max_half_maximum_error':max(width_errors),
            'one_point_probability_integral':integral,'quadrature_error_estimate':integral_error,
            'empty_data_log_likelihood':0.,'collinear_shape_log_likelihood':collinear,
            'invalid_inputs_rejected':invalid,'mathematical_checks_pass':True,
            'revised_estimator_feasibility_established':False,'real_flux_used':False,
            'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'protocol.json',HERE/'likelihood.py',HERE/'check.py']}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'timing-likelihood-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()

"""Apply known constant-beta DF positivity criterion to the lens pilot."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
H=Path(__file__).resolve().parent
fit=H.parent/'joint-galaxy-audit/results.json'
p=json.loads(fit.read_text())['sparc']['parameters']['p'];assert 0<p<1
rows=[];roots=[]
for beta in [-.3,0.,.3]:
    k=1-2*beta
    coef=[(k+2)*(k+3),2*(k+3)*(k-1),k*(k-1)]
    rr=np.roots(coef)
    threshold=max([0.]+[float(v.real) for v in rr if abs(v.imag)<1e-10 and v.real>0])
    roots.append(dict(beta=beta,k=k,exterior_polynomial_coefficients=coef,sufficient_minimum_cutoff_over_a=threshold))
    for cutoff_Re in [5.,20.,100.]:
        xt=cutoff_Re*1.8153;assert xt>threshold
        x=np.unique(np.r_[np.geomspace(1e-7,1e7,4000),xt*(1-1e-8),xt,xt*(1+1e-8)])
        L=-k/x-3/(1+x);second=L*L+k/x**2+3/(1+x)**2
        inside_bound=k*(k+1)/x**2+4*k/(x*(1+x))+6/(1+x)**2
        outside_bound=np.polyval(coef,x)/(x*x*(1+x)**2)
        bound=np.where(x<=xt,inside_bound,outside_bound)
        assert np.all(bound>0)
        for strength in [0.,1e-6,.001,1.,1000.,1e6]:
            gb=1/(1+x)**2;db=-2*gb/(1+x)
            gc=np.where(x<=xt,strength/(1+x)**(2*p),strength/(1+xt)**(2*p)*(xt/x)**2)
            dc=np.where(x<=xt,-2*p*gc/(1+x),-2*gc/x)
            slope=(db+dc)/(gb+gc)
            actual=second-L*slope
            # This is h_PsiPsi*g²/h, so it has the sign of h_PsiPsi.
            assert np.all(actual>0)
            relative_margin=(actual-bound)/np.maximum(actual,bound)
            assert relative_margin.min()>-1e-12
            rows.append(dict(beta=beta,cutoff_over_Re=cutoff_Re,relative_companion_strength=strength,
                minimum_scaled_augmented_density_curvature=float(np.min(actual*x*x)),
                minimum_relative_margin_above_analytic_bound=float(relative_margin.min())))
# Check the inverse transform's normalization with independent quadratures for
# h(Psi)=Psi^n, which is an analytic transform benchmark, not a lens fit.
transform=[]
for beta in [-.3,0.,.3]:
    lam=1.5-beta;n=4-2*beta
    C=2**(1.5-beta)*np.pi**1.5*gamma(1-beta)/gamma(lam)
    D=2**(1.5-beta)*np.pi**1.5*gamma(1-beta)*gamma(.5+beta)
    def f(E):return gamma(n+1)/(C*gamma(lam)*gamma(n-lam+1))*E**(n-lam)
    for E in [.01,.2,.8]:
        numeric=quad(lambda psi:n*(n-1)*psi**(n-2),0,E,weight='alg',wvar=(0,beta-.5),epsabs=1e-14)[0]/D
        recovered=C*quad(f,0,E,weight='alg',wvar=(0,lam-1),epsabs=1e-14)[0]
        error=abs(numeric/f(E)-1);back=abs(recovered/E**n-1)
        assert error<1e-8 and back<1e-8
        transform.append(dict(beta=beta,E=E,inverse_relative_error=error,density_recovery_relative_error=back))
out={'p':p,'cutoff_over_a_minimum_tested':5*1.8153,'analytic_sufficient_cutoff_bounds':roots,
    'numeric_force_configurations':len(rows),'radial_probes_per_configuration':4003,
    'inverse_transform_checks':transform,
    'conclusion':'Known constant-anisotropy inversion admits a nonnegative stellar tracer DF for the specified Hernquist density, positive baryonic/companion amplitudes, tested betas and cutoffs; supported by the analytic bounds in report.md.',
    'not_established':['Dynamical stability','Distribution or support of the companion deposited source','Observed orbital anisotropy','A photon-conversion field equation','Observational improvement beyond already exposed training diagnostics'],
    'fit_sha256':hashlib.sha256(fit.read_bytes()).hexdigest(),'heldout_scores_opened':False}
for name,obj in [('results.json',out),('curvature-checks.json',rows)]:
    (H/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))

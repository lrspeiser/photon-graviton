"""Connect effective collecting area to a specified ordinary-matter receiver."""
from pathlib import Path
import json,math
from scipy.integrate import quad

OUT=Path(__file__).resolve().parent
def shape(T):
    return math.sqrt(math.pi*T)*math.erf(math.sqrt(T))+math.expm1(-T)

rows=[]
for T in [.001,.1,1.,10.,1000.]:
    exact=math.pi*shape(T) # area/a^2
    numerical=2*math.pi*quad(lambda b:b*(-math.expm1(-T/(1+b*b)**2)),0,math.inf,epsabs=1e-10,epsrel=1e-10)[0]
    error=abs(numerical/exact-1)
    assert error<1e-8 and 0<shape(T)/T<=1
    ratio=4*shape(T/4)/shape(T)
    rows.append(dict(central_optical_depth=T,area_over_a_squared=exact,
        area_over_sum_of_receiver_cross_sections=shape(T)/T,
        area_ratio_when_radius_doubles_at_fixed_mass=ratio,quadrature_relative_error=error))

# Independent integration of the Plummer mass along sightlines gives column density.
column_errors=[]
for b in [0.,.1,1.,10.]:
    column=2*quad(lambda z:3/(4*math.pi)*(1+b*b+z*z)**-2.5,0,math.inf,epsabs=1e-12)[0]
    analytic=1/(math.pi*(1+b*b)**2)
    column_errors.append(abs(column/analytic-1))
assert max(column_errors)<1e-8

# Implied opacity per ordinary mass for kappa=chi*(-Phi)^4.
implied=[]
for M,a in [(1.,1.),(10.,1.),(1.,2.)]:
    for x in [0.,1.,3.,10.]:
        r=x*a
        rho=3*M*a*a/(4*math.pi*(r*r+a*a)**2.5)
        depth=M/math.sqrt(r*r+a*a) # G=chi=1
        ratio=depth**4/rho
        expected=4*math.pi*M**3/(3*a)*math.sqrt(1+x*x)
        assert abs(ratio/expected-1)<1e-13
        implied.append(dict(M=M,a=a,radius_over_a=x,implied_mass_opacity=ratio,
            relative_to_same_well_center=math.sqrt(1+x*x)))

result=dict(scope='Fixed Plummer ordinary-matter receivers, straight rays, absorption with constant cross section per unit original receiver mass. Synthetic, no excitation feedback, source depletion or observed fit.',
    provenance='Known opacity definition, Plummer projection and absorption geometry. Constant receiver mass opacity is an optional comparison to the depth-to-fourth postulate.',
    formulas=['kappa=K*rho_b; tau(b)=T/(1+(b/a)^2)^2; T=K*M/(pi*a^2)',
    'sigma_eff=pi*a^2*[sqrt(pi*T)*erf(sqrt(T))+exp(-T)-1]',
    'sigma_eff<=K*M, from 1-exp(-tau)<=tau',
    'Depth-to-fourth implies K(r)=4*pi*chi*G^4*M^3/(3*a)*sqrt(1+(r/a)^2)'],
    constant_mass_opacity_cases=rows,max_column_relative_error=max(column_errors),depth_law_implied_mass_opacity=implied)
(OUT/'receiver-consistency-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2))

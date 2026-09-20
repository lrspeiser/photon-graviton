"""RB-1 finite-front force and energy-budget necessary conditions."""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parent
G=6.67430e-11
C=299792458.
YEAR=365.25*86400
TARGET=1e10*YEAR


def main():
    out=ROOT/'radiation-budget-v1';out.mkdir(exist_ok=False)
    rows=[]
    for s,v,mass,epsilon in itertools.product((0,1),(1e5,2e5,3e5),(2e40,2e41,2e42),(1.,1e-3,1e-6)):
        active=1+s;luminosity=v*v*C**3/(G*active)
        fuel=epsilon*mass*C*C;lifetime=fuel/luminosity
        required=v*v*C*TARGET/(epsilon*G*mass)
        rows.append(dict(s=s,speed_m_per_s=v,source_kg=mass,available_fraction=epsilon,active_coefficient=active,
                         luminosity_watt=luminosity,lifetime_years=lifetime/YEAR,required_active_coefficient=required,
                         target_lifetime_passed=lifetime>=TARGET,energy_identity_error=abs(luminosity*lifetime/fuel-1)))
    profiles=[]
    for s,v,R,fraction in itertools.product((0,1),(1e5,2e5,3e5),(6e20,6e21),(.1,.5,1.)):
        active=1+s;L=v*v*C**3/(G*active);r0=R/100;r=R*fraction
        def shell_energy(x):
            density=L/(4*np.pi*x*x*C)
            return 4*np.pi*x*x*density
        energy=quad(shell_energy,r0,r,epsabs=0,epsrel=1e-12)[0]
        active_mass=quad(lambda x:active*shell_energy(x)/(C*C),r0,r,epsabs=0,epsrel=1e-12)[0]
        actual_v2=G*active_mass/r
        expected_v2=v*v*(1-r0/r)
        expected_energy=v*v*C*C*(r-r0)/(active*G)
        ev=abs(actual_v2/expected_v2-1);ee=abs(energy/expected_energy-1)
        profiles.append(dict(s=s,speed_m_per_s=v,front_m=R,radius_m=r,emission_radius_m=r0,
                             wave_energy_joule=energy,wave_equivalent_kg=energy/C**2,active_mass_kg=active_mass,
                             circular_speed_m_per_s=np.sqrt(actual_v2),force_relative_error=ev,energy_relative_error=ee,passed=ev<1e-10 and ee<1e-10))
    manifest=dict(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),sources={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ('radiation-budget-protocol.md','radiation_budget.py')},constants=dict(G=G,c=C,Julian_year_seconds=YEAR),target_years=TARGET/YEAR)
    result=dict(budgets=rows,profiles=profiles,numerical_passed=all(x['passed'] for x in profiles) and max(r['energy_identity_error'] for r in rows)<1e-10,persistent_budget_pass_count=sum(r['target_lifetime_passed'] for r in rows))
    for name,obj in (('manifest.json',manifest),('results.json',result)):
        (out/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
    typical=next(r for r in rows if r['s']==1 and r['speed_m_per_s']==2e5 and r['source_kg']==2e41 and r['available_fraction']==1)
    stored=next(r for r in profiles if r['s']==1 and r['speed_m_per_s']==2e5 and r['front_m']==6e20 and r['radius_m']==6e20)
    report=f'''# RB-1: outgoing waves can supply the radial shape, but not the budget

20 September 2026. Analytic necessary-condition screen of the weak, massless, eta=lambda=0 SR-1 branch. These are illustrative inputs, not fitted observations. The derivative coupling of wave energy to the scalar potential has coefficient 1+s relative to slow matter; it is not a new adjustable charge.

The time-averaged outgoing energy density is u=L/(4 pi r^2 c). Its gravity has the required asymptotically flat radial shape inside the wave front:

v_wave^2(r) = (1+s) G L/c^3 * (1-r0/r).

That shape alone is insufficient. Sustained luminosity must be L=v_flat^2 c^3/[(1+s)G], and its maximum fuel-limited lifetime is epsilon (1+s) G M/(v_flat^2 c).

## Results

- {len(profiles)} independent shell-quadrature force/energy comparisons pass; maximum force relative error {max(r['force_relative_error'] for r in profiles):.6g}, energy error {max(r['energy_relative_error'] for r in profiles):.6g}.
- {len(rows)} source-budget fixtures evaluated; {result['persistent_budget_pass_count']} last the declared 10 billion years.
- Across fixtures, even the longest lifetime is {max(r['lifetime_years'] for r in rows):.9g} years, using all available rest energy.
- For s=1, v=200 km/s and M=2e41 kg: required luminosity {typical['luminosity_watt']:.9g} watts. Converting the entire source rest energy lasts only {typical['lifetime_years']:.9g} years. Allowing only 1e-6 of rest energy reduces that lifetime by 1e6.
- The same example would require an active-source coefficient {typical['required_active_coefficient']:.9g} to last 10 billion years with all rest energy, versus the derived coefficient 2. This is a diagnostic discrepancy, not permission to insert a fitted multiplier.
- Storing the wave population out to 6e20 m, with r0=6e18 m, requires {stored['wave_energy_joule']:.9g} joules, equivalent to {stored['wave_equivalent_kg']:.9g} kg or {stored['wave_equivalent_kg']/2e41:.6g} of the illustrative source rest mass. Trapping eliminates rapid escape but does not remove that positive-energy inventory.

## Decision and scope

Reject sustained, freely escaping weak-wave self-gravity as the long-lived explanation for every declared fixture. This does not reject all swirl models: coherent near-field constitutive response, nonperturbative states, or a different derived matter/field coupling are outside this calculation. Each would need a new Hamiltonian, stability and source-budget audit. Slower massive waves require rederiving their source coefficient and dispersion, not simply replacing c in this formula.

The conceptual advance is precise: an outward wave population can produce an extended 1/r acceleration, but the present Hamiltonian makes that route far too expensive for the declared targets. A useful next candidate must avoid depending on continuous relativistic energy escape and demonstrate an affordable stored or constitutive response. No such candidate is validated here. A finite front is essential; an infinite stationary 1/r^2 wave-energy distribution would have infinite total energy. Source depletion can invalidate the stationary approximation on the crossing timescale; that further weakens, rather than rescues, this route.

## Attribution and reproduction

Spherical flux conservation, Poisson integration and mass-energy accounting are established mathematics. The coefficient 1+s is derived from SR-1, whose shared metric optics is attributed in its protocol. Constants: [NIST G](https://www.nist.gov/how-do-you-measure-it/how-do-you-measure-strength-gravity) and [BIPM c](https://www.bipm.org/en/si-base-units/metre). Fixture masses/speeds are hypothetical. No cosmological distances, expansion or dark-matter component enter.

Run radiation_budget.py in a clean checkout without radiation-budget-v1 to reproduce; preserve the first archive. The manifest pins source and protocol hashes and source commit. Numerical identity checks do not constitute empirical validation.
'''
    (ROOT/'radiation-budget-report.md').write_text(report,encoding='utf-8')
    print(json.dumps(dict(numerical_passed=result['numerical_passed'],budget_passes=result['persistent_budget_pass_count'],typical=typical)))


if __name__=='__main__':main()

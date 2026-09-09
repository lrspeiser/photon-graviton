from pathlib import Path
import os
import sys, json, csv, shutil
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',Path(__file__).resolve().parents[3]/'research_work/generated'))/'transport'
OUT.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(ROOT/'research_work/energy'))
from transport_ledger import TransportRates, solve, permanent_capture_fraction
from ledger_solver import Rates, solve as old_solve

checks=[]
t=np.linspace(0,20,101)
for r in [TransportRates(1,.2),TransportRates(1,.2,2,3),TransportRates(1,.2,2,3,.4,.1)]:
    got=solve(t,r,initial_photons=1,fuel=2,source=.5,duration=100)
    def rhs(j):
        h,g,ep,ec,lc,ld=r.photon_conversion,r.companion_capture,r.photon_escape,r.companion_escape,r.companion_energy_loss,r.deposit_energy_loss
        return lambda _,y:[j-(h+ep)*y[0],h*y[0]-(g+ec+lc)*y[1],g*y[1]-ld*y[2],lc*y[1]+ld*y[2],ep*y[0]+ec*y[1],-j]
    a=solve_ivp(rhs(.5),(0,4),[1,0,0,0,0,2],dense_output=True,rtol=1e-11,atol=1e-13)
    b=solve_ivp(rhs(0),(4,20),a.y[:,-1],dense_output=True,rtol=1e-11,atol=1e-13)
    expected=np.array([a.sol(x) if x<=4 else b.sol(x) for x in t])
    err=float(abs(got-expected).max()); residual=float(abs(got.sum(axis=1)-3).max())
    assert err<1e-8 and residual<1e-11 and got.min()>-1e-12
    checks.append(dict(check='independent_ODE_finite_fuel_energy_conservation',rates=r.__dict__,error=err,residual=residual))
for r in [TransportRates(1,.2),TransportRates(1,1,2,3,.4),TransportRates(0,1),TransportRates(1,0)]:
    expected=permanent_capture_fraction(r)
    got=float(solve([1000],r)[0,2])
    assert abs(got-expected)<1e-11
    checks.append(dict(check='analytic_infinite_time_capture_fraction',rates=r.__dict__,expected=expected,actual=got))
r=TransportRates(1,.2,companion_energy_loss=.3,deposit_energy_loss=.1)
a=solve(t,r,fuel=2,source=.5,duration=3)
b=old_solve(t,Rates(1,.2,.3,.1),fuel=2,injection=.5,injection_duration=3)
assert abs(a[:,[0,1,2,3,5]]-b).max()<1e-11
checks.append(dict(check='zero_escape_reduces_to_previous_solver',passed=True))
# Same total emitted energy, different source durations: the eventual permanent
# deposit cannot be multiplied merely by stretching the source lifetime.
r=TransportRates(1,.2,2,3)
yields=[float(solve([1000],r,initial_photons=0,fuel=2,source=j,duration=100)[0,2]) for j in [.1,1,10]]
assert max(yields)-min(yields)<1e-11
assert abs(yields[0]-2*permanent_capture_fraction(r))<1e-11
checks.append(dict(check='fixed_fuel_longer_emission_does_not_create_energy',source_rates=[.1,1,10],deposits=yields))

source=OUT.parent/'energy/nearby-no-loss-comparison.csv'
rows=list(csv.DictReader(source.open(encoding='utf-8-sig')))
tradeoffs=[]
for years in [1e10,1e11,1e12,1e13,1e14]:
    for supply_multiplier in [1,100,10000]:
        for capture_fraction in [1,.1,.01]:
            per_galaxy=[dict(galaxy=row['galaxy'],required_response_multiplier=float(row['full_conversion_shortfall'])*1e10/(years*supply_multiplier*capture_fraction)) for row in rows]
            tradeoffs.append(dict(duration_years=years,integrated_source_multiplier=supply_multiplier,captured_fraction=capture_fraction,median_required_response_multiplier=float(np.median([x['required_response_multiplier'] for x in per_galaxy])),galaxies=per_galaxy))
result=dict(scope='Conditional one-zone transport; no microscopic capture or new gravitational law established',checks=checks,tradeoff_definition='eta_required = full_conversion_shortfall * 1e10 years / (T * B * f). B rescales integrated available photon energy relative to present luminosity times T; f includes conversion and capture. No stellar fuel or external-source history is yet established.',scenarios=tradeoffs)
(OUT/'transport-checks-and-tradeoffs.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
shutil.copy2(ROOT/'research_work/energy/transport_ledger.py',OUT/'transport_ledger.py')
shutil.copy2(__file__,OUT/'check_transport.py')
print(json.dumps(dict(checks_passed=len(checks),scenarios=len(tradeoffs),galaxies=len(rows),max_conservation_residual=max(x.get('residual',0) for x in checks)),indent=2))

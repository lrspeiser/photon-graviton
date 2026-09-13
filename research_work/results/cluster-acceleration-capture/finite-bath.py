"""Couple calculated collecting-area growth to finite emitted energy."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.linalg import expm
OUT=Path(__file__).resolve().parent;path=OUT/'feedback-results.json'
feedback=json.loads(path.read_text());run=feedback['runs'][-1]
mass=np.r_[0.,[v['deposit_mass_over_baryonic_mass'] for v in run['snapshots']]]
gain=np.r_[1.,[v['capture_power_over_initial'] for v in run['snapshots']]]
interp=PchipInterpolator(mass,gain)

def simulate(a,growing,interpolation='pchip'):
    def rhs(source):
        def f(t,y):
            photon,companion,deposit=y
            area=(float(interp(np.clip(deposit,0,1))) if interpolation=='pchip' else float(np.interp(deposit,mass,gain))) if growing else 1.
            return [source-a*photon,a*photon-area*companion,area*companion]
        return f
    first=solve_ivp(rhs(1.),[0,1],[0.,0.,0.],method='Radau',rtol=1e-10,atol=1e-12,dense_output=True)
    second=solve_ivp(rhs(0.),[1,300],first.y[:,-1],method='Radau',rtol=1e-10,atol=1e-12,dense_output=True)
    assert first.success and second.success
    def value(t):return first.sol(t) if t<=1 else second.sol(t)
    times=np.r_[np.linspace(0,1,51),np.geomspace(1.01,300,300)]
    states=np.array([value(t) for t in times])
    ledger=np.max(abs(states.sum(axis=1)-np.minimum(times,1.)))
    assert ledger<1e-9 and np.min(states)>-1e-10, (a,growing,interpolation,ledger,float(np.min(states)))
    assert np.min(np.diff(states[:,2]))>-1e-10
    assert abs(states[-1,2]-1)<1e-8
    if not growing:
        matrix=np.array([[-a,0,0,1],[a,-1,0,0],[0,1,0,0],[0,0,0,0]],float)
        for t in [0.1,1.,3.,30.]:
            exact=expm(matrix*min(t,1))@np.array([0.,0.,0.,1.])
            if t>1:exact[:3]=expm(matrix[:3,:3]*(t-1))@exact[:3]
            assert np.max(abs(exact[:3]-value(t)))<1e-9
    half=brentq(lambda t:value(t)[2]-.5,0,300)
    ninety=brentq(lambda t:value(t)[2]-.9,0,300)
    return dict(alpha_over_initial_beta=a,growing_area=growing,time_to_half_capture=half,time_to_90percent_capture=ninety,
      final_photon_companion_deposit=states[-1].tolist(),max_energy_ledger_error=float(ledger),
      samples=[dict(time=t,photon_companion_deposit=value(t).tolist()) for t in [1.,3.,10.,30.,100.,300.]])

rows=[]
for a in [.1,1.,10.]:
    fixed=simulate(a,False);growing=simulate(a,True)
    linear=simulate(a,True,'linear')
    assert growing['time_to_90percent_capture']<fixed['time_to_90percent_capture']
    rows.append(dict(alpha_over_initial_beta=a,fixed=fixed,growing=growing,
      fractional_reduction_in_time_to_90percent=1-growing['time_to_90percent_capture']/fixed['time_to_90percent_capture'],
      linear_area_interpolation_time_to_90percent=linear['time_to_90percent_capture'],
      interpolation_t90_relative=linear['time_to_90percent_capture']/growing['time_to_90percent_capture']-1))
result=dict(scope='Homogenized identical receivers and finite positive photon injection; shared radiation bath and growing-area interpolation from prior conditional immobile-storage model. No full gravitational energy, support or observed source population.',
    input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),time_unit='1/(beta_initial*c), beta_initial=n_receiver*sigma_initial',
    energy_unit='initial ordinary receiver rest energy per receiver; radiation expressed in the same per-receiver units',
    source='j=1 for dimensionless time 0 to 1, then zero; total injected energy=1',
    area_mass_knots=mass.tolist(),area_gain_knots=gain.tolist(),cases=rows)
(OUT/'finite-bath-results.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps([dict(a=v['alpha_over_initial_beta'],fixed_t90=v['fixed']['time_to_90percent_capture'],growing_t90=v['growing']['time_to_90percent_capture'],reduction=v['fractional_reduction_in_time_to_90percent']) for v in rows],indent=2))

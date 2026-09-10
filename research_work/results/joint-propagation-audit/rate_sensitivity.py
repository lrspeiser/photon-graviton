"""Conditional rate sensitivity requested after the original delay result."""
from pathlib import Path
import math,json
from scipy.constants import c,astronomical_unit
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
P=json.loads((HERE/'protocol.json').read_text())
ROOT=HERE.parents[2]
alpha=json.loads((ROOT/P['alpha_source']).read_text())['alpha_per_mpc']
D=P['distance_mpc'];MPC=1e6*astronomical_unit*648000/math.pi
T=D*MPC/c
def excess(A):return sum(A**n/math.factorial(n+1) for n in range(1,30))
def delay(a):return T*excess(a*D)
rates=[72.417/(c/1000),alpha,76.887/(c/1000)]
# Same illustrative +/-10-second emission scenario as the original audit.
lag_upper=P['gamma_minus_gw_seconds']+P['quoted_lag_error_seconds']
def mismatch(logA):
    A=math.exp(logA)
    return T*excess(A)-lag_upper-10*math.exp(A)
logA=brentq(mismatch,-60,-10,xtol=1e-12)
A=math.exp(logA);required=A/D
assert abs(mismatch(logA))<1e-9
smallA=2*(lag_upper+10)/T
assert abs(A/smallA-1)<1e-11
result=dict(status='Post-result conditional sensitivity; no refit or measured source-lag prior.',
    rate_range_source='conversion-first/report.md: exposed training sky-tile bootstrap for c*alpha; not a physical-model posterior',
    distance_mpc=D,rate_cases=[dict(alpha_per_mpc=a,transfer_z=math.expm1(a*D),
        delay_years=delay(a)/31557600) for a in rates],
    ten_second_intrinsic_scenario=dict(max_alpha_per_mpc=required,transfer_z=math.expm1(A),
        prior_alpha_divided_by_this=alpha/required,delay_seconds=delay(required),
        observed_lag_upper_seconds=lag_upper,source_lag_bound_seconds=10),
    scope='Uniform photon-only inverse-affine law, unaffected GW, nonnegative source field age. Not a general constraint on photon-to-gravity conversion.')
(HERE/'rate-sensitivity.json').write_text(json.dumps(result,indent=2)+'\n',newline='\n')
print(json.dumps(result,indent=2))

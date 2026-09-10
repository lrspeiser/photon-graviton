"""Fixed conditional screening diagnostics; no new fit or astronomical validation."""
from pathlib import Path
import hashlib
import itertools
import json
import math

import mpmath as mp
import sympy as sp

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'conversion-first/results.json'
ALPHA = json.loads(SOURCE.read_text(encoding='utf-8'))['fitted']['alpha_per_mpc']
C = 299792458.0
AU = 149597870700.0
YEAR = 31557600.0
MPC = 1e6 * (648000 / math.pi) * AU
LY = C * YEAR
GM_E = 3.986004e14
GM_S = 1.3271244e20
R_E = 6378100.0
EPSILONS = (1e-12, 1e-9, 1e-6, .01, .1)
STARS = (1e6, 1e8, 1e10)
NS = (1, 2, 4)
mp.mp.dps = 60
errors = {'path': [], 'clock_ratio': [], 'slope': []}


def relerr(value, reference):
    return float(abs(mp.mpf(value) - reference) / abs(reference))


def terms(w, star, n, eps, sign):
    x = (w / star) ** n
    screen = 1 / (1 + x)
    logf = math.log1p(sign * eps * screen)
    logbase = .5 * math.log1p(-2 * w / C**2)
    base_slope = -1 / (C**2 - 2 * w)
    f_slope = -sign * eps * n * x * screen**2 / (w * (1 + sign * eps * screen))
    return screen, logf, logbase, base_slope, f_slope


# Independent symbolic differentiation of both environmental signs and all n.
w, star, eps, cc = sp.symbols('w star eps cc', positive=True)
symbolic_checks = 0
for n, sign in itertools.product(NS, (-1, 1)):
    x = (w / star)**n
    f = 1 + sign * eps / (1 + x)
    expected = -sign * eps * n * x / (w * (1 + x)**2 * f)
    assert sp.simplify(sp.diff(sp.log(f), w) - expected) == 0
    assert sp.simplify(sp.diff(sp.log(sp.sqrt(1-2*w/cc**2)), w) + 1/(cc**2-2*w)) == 0
    symbolic_checks += 2

paths = []
for name, length in [('1 metre', 1.), ('20200 km', 20200e3), ('1 au', AU),
                     ('100 au', 100*AU), ('1 light-year', LY),
                     ('1 million light-years', 1e6*LY), ('100 million light-years', 1e8*LY)]:
    a = ALPHA * length / MPC
    z = math.expm1(a)
    loss = -math.expm1(-a)
    # Reconstruct the length conversion at higher precision from defining constants.
    mpc_high = mp.mpf(10)**6 * (648000/mp.pi) * mp.mpf(str(AU))
    a_high = mp.mpf(str(ALPHA))*mp.mpf(str(length))/mpc_high
    errors['path'].append(relerr(z, mp.expm1(a_high)))
    assert abs((math.exp(-a) + loss) - 1) < 1e-15
    paths.append(dict(path=name, metres=length, travel_seconds=length/C,
                      optical_depth=a, conversion_z=z, photon_energy_loss_fraction=loss))

records = []
for background, epsilon, wstar, n, sign in itertools.product(
        (0., 3e10), EPSILONS, STARS, NS, (-1, 1)):
    wells = dict(earth=GM_E/R_E + GM_S/AU + background,
                 gps_height=GM_E/(R_E+20200e3) + GM_S/AU + background,
                 fifty_au=GM_S/(50*AU) + GM_E/(49*AU) + background,
                 transition=wstar)
    t = {name: terms(value, wstar, n, epsilon, sign) for name, value in wells.items()}
    ratios = {}
    for target in ('gps_height', 'fifty_au'):
        extra = math.expm1(t[target][1] - t['earth'][1])
        baseline = math.expm1(t[target][2] - t['earth'][2])
        total = math.expm1(t[target][2] - t['earth'][2] + t[target][1] - t['earth'][1])
        wh, wl = mp.mpf(str(wells[target])), mp.mpf(str(wells['earth']))
        eh, sh = mp.mpf(str(epsilon)), mp.mpf(str(wstar))
        fh = lambda ww: 1 + sign*eh/(1+(ww/sh)**n)
        ref = fh(wh)/fh(wl)-1
        errors['clock_ratio'].append(relerr(extra, ref))
        ratios[target] = dict(extra_fractional_clock_ratio=extra,
                              baseline_fractional_clock_ratio=baseline,
                              total_fractional_clock_ratio=total,
                              extra_log_ratio_over_baseline_log_ratio=(t[target][1]-t['earth'][1])/(t[target][2]-t['earth'][2]),
                              seconds_per_day_extra_on_baseline_normalized_ratio=86400*extra)
    slopes = {}
    for name, value in wells.items():
        wh = mp.mpf(str(value))
        eh, sh = mp.mpf(str(epsilon)), mp.mpf(str(wstar))
        derivative = mp.diff(lambda ww: mp.log(1+sign*eh/(1+(ww/sh)**n)), wh)
        errors['slope'].append(relerr(t[name][4], derivative))
        ratio = t[name][4]/t[name][3]
        slopes[name] = dict(extra_log_lapse_slope_per_w=t[name][4],
                            extra_to_baseline_log_lapse_slope=ratio,
                            extra_to_baseline_radial_acceleration=ratio,
                            total_clock_rate_increases_outward=(1+ratio)>0)
        assert (ratio < 0) if sign == -1 else (ratio > 0)
    records.append(dict(background_w=background, epsilon=epsilon, wstar=wstar, n=n,
                        branch='slower_void' if sign == -1 else 'faster_void_control',
                        wells=wells, earth_screen=t['earth'][0],
                        earth_conversion_rate_per_mpc=ALPHA*t['earth'][0],
                        illustrative_kappa_per_mpc=ALPHA/epsilon,
                        clock_ratios=ratios, slopes=slopes))

# Uniform-F and zero-epsilon controls: exactly no extra clock ratios or forces.
for epsilon, sign in itertools.product(EPSILONS, (-1, 1)):
    factor = 1 + sign*epsilon
    assert factor > 0 and factor/factor == 1
for star_value, n, sign in itertools.product(STARS, NS, (-1, 1)):
    t0 = terms(1e9, star_value, n, 0., sign)
    assert t0[1] == 0 and t0[4] == 0

positive_cases = 0
for epsilon, star_value, n, sign in itertools.product(EPSILONS, STARS, NS, (-1, 1)):
    for w_value in [0.] + [10**(i/10) for i in range(121)]:
        f = 1 + sign*epsilon/(1+(w_value/star_value)**n)
        assert f > 0 and (1-2*w_value/C**2) > 0
        positive_cases += 1

critical = []
for star_value, n in itertools.product(STARS, NS):
    critical_epsilon = 4*star_value/(n*(C**2-2*star_value)+2*star_value)
    tt = terms(star_value, star_value, n, critical_epsilon, -1)
    assert abs(tt[4]/tt[3]+1) < 1e-12
    critical.append(dict(wstar=star_value, n=n,
                         epsilon_at_zero_total_gradient_at_transition=critical_epsilon))

maxerr = {name: max(values) for name, values in errors.items()}
assert maxerr['path'] < 1e-12
assert maxerr['clock_ratio'] < 1e-7
assert maxerr['slope'] < 1e-7
assert len(records) == 180
output = dict(status='conditional fixed-grid diagnostic, not an empirical screening validation',
              source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
              protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
              alpha_per_mpc=ALPHA, fractional_energy_loss_rate_per_second=ALPHA*C/MPC,
              constants=dict(c_m_s=C, au_m=AU, julian_year_s=YEAR, mpc_m=MPC,
                             nominal_gm_earth=GM_E, nominal_gm_sun=GM_S, nominal_earth_radius=R_E),
              path_cases=paths, screening_cases=records, transition_critical_amplitudes=critical,
              checks=dict(symbolic=symbolic_checks, high_precision_clock_ratios=len(errors['clock_ratio']),
                          high_precision_slopes=len(errors['slope']), positive_lapse=positive_cases,
                          max_relative_errors=maxerr, uniform_factor_cancellation=True, epsilon_zero_recovery=True))
(HERE/'results.json').write_text(json.dumps(output, indent=2, allow_nan=False)+'\n', encoding='utf-8')

lines = ['# Numerical tables: void-screening diagnostic', '',
         'Generated by `run.py`; interpretation and assumptions are in [report.md](report.md) and [protocol.md](protocol.md).', '',
         '## Unscreened photon conversion at the reused exploratory coefficient', '',
         '| Path | Conversion redshift | Photon energy fraction transferred |',
         '|---|---:|---:|']
for row in paths:
    lines.append(f"| {row['path']} | {row['conversion_z']:.8g} | {row['photon_energy_loss_fraction']:.8g} |")
lines += ['', '## All fixed screening cases', '',
          'Ratios are hypothetical stationary-clock comparisons, not satellite fits. R is the extra-to-baseline log-lapse gradient ratio; R below -1 reverses the baseline gradient in the stated coupling.', '',
          '| Background W | Branch | epsilon | W_star | n | Extra Earth-to-GPS-height clock ratio | R at Earth | R at transition |',
          '|---:|---|---:|---:|---:|---:|---:|---:|']
for r in records:
    lines.append(f"| {r['background_w']:.3g} | {r['branch']} | {r['epsilon']:.3g} | {r['wstar']:.3g} | {r['n']} | {r['clock_ratios']['gps_height']['extra_fractional_clock_ratio']:.6g} | {r['slopes']['earth']['extra_to_baseline_log_lapse_slope']:.6g} | {r['slopes']['transition']['extra_to_baseline_log_lapse_slope']:.6g} |")
lines += ['', '## Verification', '', '```json', json.dumps(output['checks'], indent=2), '```', '']
(HERE/'tables.md').write_text('\n'.join(lines), encoding='utf-8')
print(json.dumps(dict(paths=paths, checks=output['checks'], cases=len(records)), indent=2))

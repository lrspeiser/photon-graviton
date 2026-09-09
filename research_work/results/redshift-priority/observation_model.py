"""Factor accounting and sensitivity only; no per-object corrections are fitted."""
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
C = 299792.458
ALPHA = 0.0002488993286382367  # Previously calibrated, unchanged.


def velocity_to_z(value, convention):
    b = value / C
    if convention == 'optical':
        if b <= -1:
            raise ValueError('Nonpositive spectral factor')
        return b
    if convention == 'radio':
        if b >= 1:
            raise ValueError('Nonpositive spectral factor')
        return b / (1-b)
    if convention == 'relativistic_radial':
        if abs(b) >= 1:
            raise ValueError('Radial speed must be subluminal')
        return math.expm1(math.atanh(b))
    raise ValueError('Explicit known velocity convention required')


def z_to_velocity(z, convention):
    if z <= -1:
        raise ValueError('Nonpositive spectral factor')
    if convention == 'optical':
        return C*z
    if convention == 'radio':
        return C*z/(1+z)
    if convention == 'relativistic_radial':
        return C*math.tanh(math.log1p(z))
    raise ValueError('Explicit known velocity convention required')


def combine(time_depth, motion_depth=0., endpoint_depth=0., measurement_depth=0.):
    return math.expm1(time_depth + motion_depth + endpoint_depth + measurement_depth)


def main():
    source = ROOT/'redshift_paper/all_164_groups.csv'
    rows = list(csv.DictReader(source.open(encoding='utf-8-sig')))
    result = []
    for row in rows:
        d, z = float(row['catalog_distance_mpc']), float(row['observed_cmb_z'])
        predicted_depth = ALPHA*d
        predicted_z = combine(predicted_depth)
        excess = math.log1p(z)-predicted_depth
        record = dict(pgc=row['pgc'], group_pgc=row['group_pgc'], distance_mpc=d,
            observed_z=z, fixed_baseline_z=predicted_z,
            observed_minus_baseline_log_factor=excess,
            equivalent_radial_motion_kms=C*math.tanh(excess))
        # Hypothetical sensitivity scales only. Not priors or confidence intervals.
        for speed in (100, 300, 1000):
            motion_depth = math.atanh(speed/C)
            record[f'scenario_plus_{speed}_kms_z'] = combine(predicted_depth, motion_depth)
            record[f'scenario_minus_{speed}_kms_z'] = combine(predicted_depth, -motion_depth)
            record[f'plus_{speed}_additive_shortcut_error_cdz_kms'] = C*(combine(predicted_depth,motion_depth)-predicted_z-speed/C)
        result.append(record)
    with (HERE/'observation-factor-diagnostics.csv').open('w', newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(result[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(result)
    checks=0
    for convention in ('optical','radio','relativistic_radial'):
        for z in (-.5,-.01,0,.01,.1,1,5):
            recovered=velocity_to_z(z_to_velocity(z,convention),convention)
            assert math.isclose(recovered,z,rel_tol=2e-14,abs_tol=1e-15)
            checks+=1
    # Factor cross term and round-trip inversion, including blueshift and cancellation.
    for t in (0.,.001,.1,1.):
        for m in (-.1,0.,.1):
            z=combine(t,m,.002,-.001)
            assert math.isclose(math.log1p(z)-m-.002+.001,t,abs_tol=1e-14)
            assert math.isclose(1+z,math.exp(t)*math.exp(m)*math.exp(.002)*math.exp(-.001),rel_tol=1e-14)
            checks+=1
    for v, convention in [(C,'radio'),(C,'relativistic_radial'),(-C,'relativistic_radial'),(-C,'optical'),(0,'unspecified')]:
        try:
            velocity_to_z(v, convention)
        except ValueError:
            checks+=1
        else:
            raise AssertionError('Invalid convention/input accepted')
    summary=dict(rows=len(rows), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        baseline_alpha_per_mpc=ALPHA, numerical_checks=checks,
        convention_status='CF4 CMB frame confirmed; exact original velocity convention and frame-transform implementation not certified by audited appendix.',
        residual_status='Equivalent radial speeds are residual diagnostics, never measured motions or fitted corrections.',
        scenario_status='100/300/1000 km/s are illustrative fixed sensitivities, not independently justified population widths or predictive intervals.',
        positive_300_additive_shortcut_error_cdz_kms_range=[min(r['plus_300_additive_shortcut_error_cdz_kms'] for r in result),max(r['plus_300_additive_shortcut_error_cdz_kms'] for r in result)])
    (HERE/'observation-factor-results.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()

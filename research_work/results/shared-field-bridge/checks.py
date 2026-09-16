"""Suite job for BRIDGE-1B (see protocol.md). At reduced size it checks:
- the crossing spectrum against the Landau-Zener formula (V1);
- energy conservation and the back-reaction's stopping point (V2, V3);
- the owner's four controls (V5);
- the exponent identity behind the atomic-reference factor;
- the physical mapping's internal relations at the reference point;
- the driver's short-run regression anchor against the archived results.

The Pantheon+ covariance is not loaded here; the brightness numbers belong to the driver."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bridge as B  # noqa: E402  (also puts modes on the path)


def main():
    v1 = B.trial_crossing_spectrum()
    trap = B.trial_trapping(4.)
    controls, control_checks, _ = B.trial_reference_and_controls()
    ident = B.exponent_identity()
    ref = B.scenario(B.REF_NSTAR, B.REF_RE)
    h = B.history(B.REF_NSTAR, B.REF_RE)
    rel = lambda a, b: abs(a - b)/max(abs(b), 1e-300)
    mapping = dict(
        k_star_squared_over_g_ndot=rel(ref['k_star_eV']**2, ref['g_eV']*h['ndot_star']*B.HBAR_EVS),
        abundance=rel(ref['number_density_eV3']*ref['g_eV']*(1 - B.REF_NSTAR), B.RHO_C0*B.EV4_PER_J_M3),
        speed_law=rel(ref['v_rms_today'], ref['k_rms_eV']/ref['mass_today_eV']),
        rolls_today=bool(ref['trapping_distance_in_n'] > ref['needed_span_in_n']))
    mapping['passed'] = bool(max(mapping['k_star_squared_over_g_ndot'], mapping['abundance'],
                                 mapping['speed_law']) < 1e-9 and mapping['rolls_today'])
    energy = max([t['energy_error'] for t in (trap,)] + [c['energy_error'] for c in controls])
    out = dict(V1=dict(max_absolute_error=v1['max_absolute_error'],
                       max_relative_error=v1['max_relative_error_above_0p01'], passed=v1['passed']),
               V2=dict(worst_relative_energy_error=energy, passed=bool(energy <= 1e-8)),
               V3=dict(stop_relative_error=trap['stop_relative_error'], passed=trap['passed']),
               V5=control_checks, identity=dict(max_residual=ident['max_residual'], passed=ident['holds']),
               mapping=mapping)
    out['V5']['passed'] = bool(all(v for k, v in control_checks.items() if isinstance(v, bool)))
    arch = json.loads((HERE/'bridge-results.json').read_text(encoding='utf-8')).get('checks_short_run', {})
    recomputed = dict(v1_absolute_error=v1['max_absolute_error'], v1_number=v1['number_density'],
                      trap_lambda4_stop=trap['stop_n'], reference_companion_energy=controls[0]['companion_energy'],
                      reference_field_kinetic=controls[0]['field_kinetic_energy'],
                      reference_k_star_eV=ref['k_star_eV'], reference_mass_today_eV=ref['mass_today_eV'],
                      identity_max_residual=ident['max_residual'])
    same = bool(arch) and all(abs(recomputed[k] - arch[k]) <= 1e-9*max(1., abs(arch[k])) for k in recomputed)
    out.update(short_run=recomputed, short_run_matches_archive=same)
    out['passed'] = bool(all(out[k]['passed'] for k in ('V1', 'V2', 'V3', 'V5', 'identity', 'mapping')) and same)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())

"""Suite job for PM-1 (see protocol.md). At reduced size it checks:
- candidate A's control: a rigidly rotating axisymmetric source is the static one, so delay sees no orbit;
- candidate B's fitted constant and its training score, and that the banded form agrees with its limit;
- the sourcing rules' mass-speed slopes, which is the gate that separates them;
- candidate C's cubic against the owner's root, and that it has a growing mode;
- the driver's short-run regression anchor against the archived results."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pm1 as P  # noqa: E402


def main():
    gals = P.sparc()
    mw = P.milky_way()
    a = P.candidate_a(gals[0])
    root = P.candidate_b(gals, mw, 'root', (1e-2, 1e-1))
    linear = P.candidate_b(gals, mw, 'linear', (1e-7, 1e-6))
    collective = P.candidate_b(gals, mw, 'collective', (1e-2, 1e-1))
    ring = P.ring_audit()
    rep = P.ring_representation(max(gals, key=lambda g: len(g['R'])))
    c = P.candidate_c(gals)
    out = dict(
        A=dict(rotating_minus_static=a['rotating_minus_static_max'], passed=a['passed']),
        B=dict(constant=root['constant'], train=root['train'],
               a_star_m_s2=root['implied_a_star_m_s2'],
               beats_baryons=root['beats_baryons'],
               passed=bool(root['beats_baryons'] and root['train'] < P.BARYONS['train'])),
        mass_speed=dict(observed=root['mass_speed_slope_observed'], root=root['mass_speed_slope_model'],
                        linear=linear['mass_speed_slope_model'],
                        passed=bool(abs(root['mass_speed_slope_model'] - root['mass_speed_slope_observed'])
                                    < abs(linear['mass_speed_slope_model'] - root['mass_speed_slope_observed']))),
        C=dict(owner_root=c['owner_root_check'], worst_growth=c['worst_growth_rate'],
               passed=bool(abs(c['owner_root_check']['got'][0] - .1533) < 1e-3
                           and abs(c['owner_root_check']['got'][1] - 1.2858) < 1e-3
                           and not c['stable_anywhere'])))
    out['ring'] = dict(source_integrated=ring['owner_example']['source_integrated'],
                       on_ring=ring['owner_example']['on_ring_at_r'],
                       reproduces_owner=ring['reproduces_owner'],
                       representation_growth=rep['growth_over_8x_refinement'],
                       expected_sqrt_8=rep['expected_sqrt_8'],
                       passed=bool(ring['reproduces_owner']
                                   and abs(rep['growth_over_8x_refinement']/rep['expected_sqrt_8'] - 1) < 1e-6))
    arch = json.loads((HERE/'pm1-results.json').read_text(encoding='utf-8')).get('checks_short_run', {})
    now = dict(root_constant=root['constant'], root_train=root['train'], root_test=root['test'],
               linear_train=linear['train'], collective_train=collective['train'],
               c_worst_growth=c['worst_growth_rate'],
               ring_source_integrated=ring['owner_example']['source_integrated'])
    same = bool(arch) and all(abs(now[k] - arch[k]) <= 1e-9*max(1., abs(arch[k])) for k in now)
    out.update(short_run=now, short_run_matches_archive=same)
    out['passed'] = bool(all(out[k]['passed'] for k in ('A', 'B', 'mass_speed', 'C', 'ring')) and same)
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())

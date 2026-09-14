"""Quick RPG-1 checks for run_checks.py: the analytic solver validations (V1 a-c) and the three
convergence galaxies against the archived rows in rpg1-results.json.

    python checks.py
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rpg1  # noqa: E402


def main():
    archived = json.loads((HERE/'rpg1-results.json').read_text(encoding='utf-8'))
    v1 = rpg1.validation()
    diffs = rpg1.evidence_io.compare(v1, {k: v for k, v in archived['V1'].items() if k != 'sparc_convergence'},
                                     rtol=1e-6, atol=1e-12)
    rpg1._init()
    rows = {r['name']: r for r in archived['galaxies']}
    speed = {}
    for name in archived['convergence_selection'].values():
        new = rpg1.solve_galaxy((name, .1, 1, 1., False))
        for key in ('newtonian_kms', 'algebraic_kms', 'aqual_kms'):
            speed[f'{name}/{key}'] = float(np.max(np.abs(np.array(new[key]) - np.array(rows[name][key]))))
    ok = all(v['passed'] for v in v1.values()) and not diffs and max(speed.values()) < 1e-6
    print(json.dumps(dict(V1=v1, V1_regression_differences=diffs, max_speed_differences_kms=speed, passed=ok), indent=1))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())

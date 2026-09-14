"""Quick CR-2 checks for run_checks.py:
- the Thomas–Fermi solver against sin(x)/x (V1);
- the reproduction of the archived FLRW stars-only and NFW exact-lens fits for every lens (V2);
- one archived condensate fit per geometry, recomputed at its archived parameters.

    python checks.py
"""
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cr2  # noqa: E402
import tf  # noqa: E402


def main():
    np.seterr(over='raise', invalid='raise', divide='raise')
    archived = json.loads((HERE/'cr2-results.json').read_text(encoding='utf-8'))
    x = np.geomspace(1e-5, 1e5, 101)
    s0 = tf.Family(x, x**3/(1 + x**3), 1.).solve(0.)
    v1 = bool(abs(s0['x_edge'] - math.pi) < 1e-8 and abs(s0['m_edge'] - math.pi) < 1e-8)
    v2 = {}
    for name in cr2.LENSES:
        lens = cr2.Lens(name, 'G1_flat_FLRW')
        for kind, value in (('stars_only_exact_lens', lens.stars_only()['chi2']), ('NFW_exact_lens', lens.nfw()['chi2'])):
            arc = next(r for r in cr2.ARCHIVE if r['Name'] == name and r['model'] == kind)['stellar_chi2']
            v2[f'{name} {kind}'] = abs(value/arc - 1)
    refits = {}
    for which, g in archived['geometries'].items():
        row = g['M1']['rows'][0]
        lens = cr2.Lens(row['name'], which)
        q = lens.m1_at(tf.Family(lens.model.r, lens.frac, tf.k_of(g['M1']['best_R_tf_kpc'])), row['lam'])
        refits[which] = abs(q['chi2']/row['chi2'] - 1)
    ok = v1 and max(v2.values()) < 1e-4 and max(refits.values()) < 1e-6
    print(json.dumps(dict(V1=v1, V2_max_relative_difference=max(v2.values()), refit_relative_differences=refits, passed=ok), indent=1))
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())

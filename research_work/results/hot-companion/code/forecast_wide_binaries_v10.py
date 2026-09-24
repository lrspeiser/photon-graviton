"""A locked forecast (proposal 2): wide binary stars near the Sun, for Gaia's fourth data release.

Written and committed before those data are examined. The law is frozen at the adopted round-9 form
(regression/law_config.py, 'round9': the round-3 constants with the companion released over
30,000 AU), and the Galaxy's Newtonian pull at the Sun comes from turning the law inside out for
230 km/s at 8.2 kpc, as in the regression suite (t_precision.galactic_pull).

Forecast quantity: the orientation-averaged pull between the two stars divided by Newton's,
gamma(s) = 1 + M_ph(<s)/M, at 3D separation s, for total masses 1.0, 1.5 and 2.0 Msun; and the
corresponding ratio of orbital speeds, sqrt(gamma). The candidate 'no_hold' route for the dwarfs
(round 9, not adopted) is given alongside, labelled, because the binaries separate the two.

    python code/forecast_wide_binaries_v10.py --output forecasts/wide_binaries_gaia_dr4_v10.json
"""
from __future__ import annotations
import argparse, hashlib, json, subprocess, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'regression'))
import common  # noqa: E402,F401  (paths)
import t_precision as P  # noqa: E402
from law_config import load_law  # noqa: E402

SEPS = (1000, 2000, 3000, 5000, 7000, 10000, 15000, 20000, 30000, 50000)
MASSES = (1.0, 1.5, 2.0)


def profile(law, hold=None, L=None):
    a, gd = law['a_SI'], law['g_d_SI']
    c = law['external_hold'] if hold is None else hold
    Lau = law['release_length_au'] if L is None else L
    ge, gobs = P.galactic_pull(law)
    out = {}
    for m in MASSES:
        _, b = P.sun_in_galaxy(m * P.MSUN, c * ge, 0.0, 0.0, Lau, a, gd, seps_au=SEPS)
        out[f'{m:g}'] = {str(s): dict(gamma=round(b[s], 4), speed_ratio=round(b[s] ** 0.5, 4)) for s in SEPS}
    return out, ge, gobs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    law = load_law('round9')
    adopted, ge, gobs = profile(law)
    cand, _, _ = profile(law, hold=0.0, L=200000.0)
    body = dict(
        forecast='wide binary stars near the Sun (3D separation s), for Gaia DR4',
        law=dict(name=law['name'], description=law['description'], a_SI=law['a_SI'], g_d_SI=law['g_d_SI'],
                 u_kms=law['u_kms'], release_length_au=law['release_length_au'], external_hold=law['external_hold']),
        galaxy_newtonian_pull_at_sun_SI=ge, galaxy_observed_pull_at_sun_SI=gobs,
        quantity='gamma = pull between the stars / Newton (orientation-averaged); speed_ratio = sqrt(gamma)',
        adopted_law=adopted,
        candidate_no_hold_not_adopted=dict(note='external hold 0, release over 200,000 AU (round 9 candidate for the dwarfs)', values=cand),
        what_would_count=('the adopted law is supported if the measured gamma is 1.03-1.05 at s = 7,000 AU and 1.08-1.10 at '
                          '20,000 AU for solar-mass pairs (within the analysis errors); it is refuted by gamma < 1.02 or > 1.2 '
                          'at 20,000 AU with errors below 0.03, or by a rise before 3,000 AU'),
        comparison=dict(newton='gamma = 1', mond_simple_quoted_in_round_6='about 1.43 at 10,000-30,000 AU'))
    blob = json.dumps(body, sort_keys=True).encode()
    body['sha256_of_forecast'] = hashlib.sha256(blob).hexdigest()
    try:
        body['written_at_commit'] = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=HERE, capture_output=True, text=True).stdout.strip()
    except Exception:
        pass
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(body, indent=1) + '\n')
    print(json.dumps({m: {s: v['gamma'] for s, v in d.items()} for m, d in adopted.items()}, indent=0))
    print('candidate no_hold:', {m: {s: v['gamma'] for s, v in d.items()} for m, d in cand.items()}['1'])
    print('sha256', body['sha256_of_forecast'])


if __name__ == '__main__':
    main()

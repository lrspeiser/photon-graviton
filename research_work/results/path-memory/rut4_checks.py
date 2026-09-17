"""Suite job for RUT-1 stage 4 (protocol-rut1.md, stage 4). rut4.py is the science driver and takes tens
of minutes on many cores; this job checks, in about a minute, that its archive is what the committed code
produces:

- the committed formation.py, longrun.py, rut4.py, rut3.py and rut1.py hash to the values recorded when the
  campaign was LAUNCHED, so the code on the branch is the code that ran;
- every committed series file hashes to the value recorded in rut4-results.json;
- a deterministic 5 T0 prefix of the representative two-stage run is replayed and must reproduce the
  archived 5 T0 checkpoint -- positions, velocities, the eps residual, and the every-step and per-band work
  and torque;
- the whole analysis -- trends, numerical uncertainty, classification, common-time tables, the primed
  challenge and the gates -- is recomputed from the committed series and must reproduce the archive.
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import rut4 as R4  # noqa: E402

TOL = 1e-9


def _load(name):
    with gzip.open(HERE/'rut4-series'/f'{name}.json.gz', 'rt', encoding='utf-8') as f:
        return json.load(f)


def _close(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return bool(a.shape == b.shape and np.all(np.abs(a - b) <= TOL*np.maximum(1., np.abs(b))))


def main():
    archive = json.loads((HERE/'rut4-results.json').read_text(encoding='utf-8'))
    out = {}

    code_now = {f: hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in archive['code_sha256_at_launch']}
    mismatched = [f for f, h in code_now.items() if h != archive['code_sha256_at_launch'][f]]
    out['code_is_what_ran'] = dict(mismatched=mismatched, passed=not mismatched)

    series_bad = [n for n, h in archive['series_sha256'].items()
                  if hashlib.sha256((HERE/'rut4-series'/f'{n}.json.gz').read_bytes()).hexdigest() != h]
    out['series_hashes'] = dict(mismatched=series_bad, passed=not series_bad)

    name = 'two_stage_32_w0.2_s1'
    spec = next(s for s in archive['runs'] if s['name'] == name)
    replay = R4.execute(spec, horizon=5., checkpoints=(5,))
    stored = _load(name)['checkpoints']['5']
    fresh = replay['checkpoints']['5']
    keys = ('positions', 'velocities', 'radius', 'L_ratio', 'eps_residual', 'W_total', 'J_total')
    prefix = {k: _close(fresh[k], stored[k]) for k in keys}
    prefix.update({f'W_band:{b}': _close(fresh['W_band'][b], stored['W_band'][b]) for b in stored['W_band']})
    prefix.update({f'J_band:{b}': _close(fresh['J_band'][b], stored['J_band'][b]) for b in stored['J_band']})
    out['deterministic_prefix'] = dict(run=name, horizon_T0=5., checked=len(prefix),
                                       failed=[k for k, v in prefix.items() if not v], passed=all(prefix.values()))

    series = {s['name']: _load(s['name']) for s in archive['runs']}
    analyses, unc, classes, common, primed, rotation = R4.analyze(series)
    gate = R4.gates(series, analyses, unc, rotation)
    for a in analyses.values():
        a.pop('wall_seconds', None)
    recomputed = json.loads(json.dumps(dict(analyses=analyses, numerical_uncertainty=unc, classification=classes,
                                            common_time_comparison=common, primed_challenge=primed, gates=gate),
                                       default=float))
    stored_part = {k: archive[k] for k in recomputed}
    diffs = evidence_io.compare(recomputed, stored_part)
    out['analysis_reproduces'] = dict(differences=diffs[:10], count=len(diffs), passed=not diffs)

    out['archive_gates'] = {k: (v or {}).get('passed') for k, v in archive['gates'].items()}
    out['passed'] = bool(all(v['passed'] for k, v in out.items() if isinstance(v, dict) and 'passed' in v)
                         and archive['passed'])
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())

"""Audit corrected sampling against existing annulus constructions; no databases needed.

    python -m research_work.annulus_sampling.audit --all \
        --output research_work/generated/annulus-sampling-v2.json

Writes a new report, never a canonical historical archive. A numerical or
construction failure exits nonzero. Neither a good sampler nor an equilibrium
fixed-point residual is a stability result.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any

import numpy as np
import scipy

from .adapter import ARCHIVED_PATH, WarmAnnulus
from .sampling import VERSION, verify_distribution

ROOT = Path(__file__).resolve().parents[2]
CASES = {'narrow_cold': (.06, .010), 'narrow_warm': (.06, .030),
         'wide_cold': (.12, .010), 'wide_warm': (.12, .030)}


def new_output_path(path: Path) -> Path:
    path = path.resolve()
    if path.exists():
        raise FileExistsError(f'Refusing to overwrite {path}')
    if path.is_relative_to((ROOT/'research_work/results').resolve()):
        raise ValueError('Write to research_work/generated or another scratch directory, not frozen results')
    if path.suffix.lower() != '.json':
        raise ValueError('Report output must have a .json suffix')
    return path


def write_new_report(path: Path, report: dict[str, Any]) -> None:
    path = new_output_path(path)
    # Serialization precedes opening the file, so NaN cannot leave a partial scientific report.
    text = json.dumps(report, indent=2, allow_nan=False) + '\n'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        f.write(text)


def source_fingerprint(a: WarmAnnulus) -> str:
    h = hashlib.sha256()
    for name in ('r', 'sigma', 'C'):
        arr = np.asarray(getattr(a, name), dtype='<f8')
        h.update(name.encode()); h.update(arr.tobytes())
    h.update(json.dumps({k: float(getattr(a, k)) for k in
                         ('L0', 'dL', 'dE', 'w', 'tau_keep', 'alpha')}, sort_keys=True).encode())
    return h.hexdigest()


def audit_case(name: str, n_r: int = 280) -> dict[str, Any]:
    dL, dE = CASES[name]
    # This deliberately repeats the archived matched-support construction. It is
    # NOT a fixed-mass/fixed-alpha temperature experiment and does not rescale a draw.
    a = WarmAnnulus(L0=1., dL=dL, dE=dE, w=.2, n_r=n_r).solve(target_support=.1)
    summary = a.summary()
    consistency = a.consistency_residual()
    if not all(np.isfinite(v) for v in consistency.values()) or max(consistency.values()) > 1e-6:
        raise RuntimeError(f'{name}: equilibrium fixed-point check failed: {consistency}')
    if not np.isfinite(summary['mass']) or summary['mass'] <= 0:
        raise RuntimeError(f'{name}: invalid constructed mass')
    before = source_fingerprint(a)
    nodes = a.sampling_distribution()
    verification = verify_distribution(a, nodes)
    biased = replace(nodes, node_mass=nodes.node_mass*nodes.node_radius)
    bias_verification = verify_distribution(a, biased)
    moments = nodes.moments()['radius']
    shift = biased.moments()['radius']['mean'] - moments['mean']
    expected_shift = moments['variance']/moments['mean']
    unchanged = before == source_fingerprint(a)
    passed = (verification['numerical_verification_passed'] and unchanged
              and not bias_verification['numerical_verification_passed']
              and abs(shift-expected_shift) <= 1e-12)
    return {'case': name, 'numerical_verification_passed': bool(passed),
            'constructed_state': summary, 'fixed_point_residual': consistency,
            'total_writing_rate': float(a.alpha*summary['mass']),
            'requested_support_at_density_peak': .1,
            'achieved_support_at_density_peak': summary['support_fraction_at_peak'],
            'sampler': verification, 'source_unchanged': unchanged,
            'source_fingerprint': before,
            'deliberate_radius_bias_control': {
                'rejected': not bias_verification['numerical_verification_passed'],
                'biased_mean_radius': biased.moments()['radius']['mean'],
                'mean_radius_shift': shift, 'variance_over_mean': expected_shift,
                'radial_marginal_error': bias_verification['radial_marginal_max_absolute_error']},
            'quadrature': {'radial_nodes': n_r, 'L_nodes': 160, 'positive_vr_nodes': 96,
                           'reach_in_widths': 5., 'radial_rule': 'trapezoid',
                           'L_rule': 'archived equal rectangular weights'},
            'scope': 'Discrete sampling verified; continuum/domain convergence and stability not evaluated.'}


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--all', action='store_true', help='Audit all four archived construction families')
    p.add_argument('--cases', nargs='+', choices=tuple(CASES), default=['narrow_cold', 'narrow_warm'])
    p.add_argument('--n-r', type=int, default=280, help='Radius-grid size; defaults to the archived 280')
    p.add_argument('--output', type=Path, required=True, help='New JSON path; existing files are refused')
    args = p.parse_args(argv)
    try:
        output = new_output_path(args.output)
    except (ValueError, OSError) as exc:
        p.error(str(exc))
    if args.n_r < 3:
        p.error('--n-r must be >= 3')
    names = list(CASES) if args.all else list(dict.fromkeys(args.cases))
    rows, failures = [], []
    for name in names:
        try:
            row = audit_case(name, args.n_r)
            rows.append(row)
            print(f"{name}: {'PASS' if row['numerical_verification_passed'] else 'FAIL'}", flush=True)
        except Exception as exc:
            failures.append({'case': name, 'error_type': type(exc).__name__, 'message': str(exc)})
            print(f'{name}: FAIL {exc}', file=sys.stderr, flush=True)
    passed = (not failures and len(rows) == len(names)
              and all(r['numerical_verification_passed'] for r in rows))
    report = {'sampler_version': VERSION,
              'status': {'numerical_verification': 'passed' if passed else 'failed',
                         'historical_archive_reproduction': 'not_run',
                         'historical_H6': 'failed_preserved_unchanged',
                         'physical_stability': 'not_evaluated',
                         'continuum_convergence': 'not_evaluated'},
              'requested_cases': names, 'completed_cases': len(rows),
              'failures': failures, 'cases': rows,
              'provenance': {'python': platform.python_version(), 'numpy': np.__version__,
                             'scipy': scipy.__version__,
                             'archived_equilibrium_sha256': hashlib.sha256(ARCHIVED_PATH.read_bytes()).hexdigest(),
                             'implementation_sha256': {
                                 name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                                 for name in ('sampling.py', 'adapter.py', 'audit.py')}}}
    write_new_report(output, report)
    print(f'Report: {output}')
    return 0 if passed else 1


if __name__ == '__main__':
    raise SystemExit(main())

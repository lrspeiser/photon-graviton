"""Entry-point handling shared by diagnostic scripts.

Regenerated results go to a fresh directory (research_work/generated/ by
default, or the run_checks.py output directory); archived canonical files are
only overwritten with an explicit --canonical. Every run compares its numbers
with the archived file, so the preserved assertions and outputs act as a
regression check.
"""
from pathlib import Path
import argparse
import json
import math
import os
import re
import time
import uuid

ROOT = Path(__file__).resolve().parents[3]


def parse(description):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument('--output-dir', type=Path,
                    help='Directory for regenerated results (default: $PHOTON_GRAVITON_RESULTS/<name> or research_work/generated/<run>).')
    ap.add_argument('--canonical', action='store_true', help='Also overwrite the archived canonical result file.')
    return ap.parse_args()


def output_dir(args, name):
    if args.output_dir:
        d = args.output_dir
    elif os.environ.get('PHOTON_GRAVITON_RESULTS'):
        d = Path(os.environ['PHOTON_GRAVITON_RESULTS'])/name
    else:
        d = ROOT/'research_work/generated'/(time.strftime('%Y%m%d-%H%M%S') + '-' + name + '-' + uuid.uuid4().hex[:6])
    d.mkdir(parents=True, exist_ok=True)
    return d


def compare(new, old, path='', ignore=(), rules=(), rtol=1e-9, atol=1e-12):
    """Paths at which two parsed JSON documents differ numerically or structurally.

    ignore: exact paths skipped entirely. rules: (regex, rtol, atol) tuples; the
    first pattern matching a path sets the tolerance for it and its children,
    and rtol=None skips it. Rules are only for quantities an optimizer does not
    determine uniquely; each use must say why.
    """
    if path in ignore:
        return []
    for pattern, r, a in rules:
        if re.match(pattern, path):
            if r is None:
                return []
            rtol, atol = r, a
            break
    if isinstance(new, dict) and isinstance(old, dict):
        diffs = [f'{path}/{k}: missing' for k in set(new) ^ set(old) if f'{path}/{k}' not in ignore]
        for k in set(new) & set(old):
            diffs += compare(new[k], old[k], f'{path}/{k}', ignore, rules, rtol, atol)
        return diffs
    if isinstance(new, list) and isinstance(old, list):
        if len(new) != len(old):
            return [f'{path}: length {len(new)} vs {len(old)}']
        return [d for i, (a, b) in enumerate(zip(new, old))
                for d in compare(a, b, f'{path}[{i}]', ignore, rules, rtol, atol)]
    if isinstance(new, bool) or isinstance(old, bool) or not isinstance(new, (int, float)) or not isinstance(old, (int, float)):
        return [] if new == old else [f'{path}: {new!r} vs {old!r}']
    return [] if math.isclose(new, old, rel_tol=rtol, abs_tol=atol) else [f'{path}: {new!r} vs {old!r}']


def finish(args, name, text, canonical, ignore=(), rules=()):
    """Write regenerated results, compare with the archived file, optionally replace it."""
    out = output_dir(args, name)/canonical.name
    out.write_text(text, encoding='utf-8', newline='\n')
    diffs = compare(json.loads(text), json.loads(canonical.read_text(encoding='utf-8')), ignore=set(ignore), rules=rules)
    if args.canonical:
        canonical.write_text(text, encoding='utf-8', newline='\n')
    print(json.dumps(dict(output=str(out), canonical=str(canonical), canonical_rewritten=args.canonical,
                          n_regression_differences=len(diffs), first_differences=diffs[:20]), indent=1))
    return 0 if args.canonical or not diffs else 1

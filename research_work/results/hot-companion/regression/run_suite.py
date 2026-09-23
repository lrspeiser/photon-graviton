#!/usr/bin/env python3
"""Run every test of the law, score it, and compare with the saved baseline.

    python run_suite.py                               # the round-3 law, quick tier (about 1 minute)
    python run_suite.py --tier full                   # + five colliding clusters (about 20 minutes)
    python run_suite.py --law gd_x1p5                 # a candidate from candidates/
    python run_suite.py --only dwarfs,precision       # some groups only
    python run_suite.py --tier full --save-baseline   # make this run the new baseline

Exit code 1 when any check got worse in status (pass -> close/fail, close -> fail) or crashed,
so the suite can guard every change. See README.md in this folder.
"""
from __future__ import annotations
import argparse, datetime, importlib, json, subprocess, sys, time, traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import common as C                    # noqa: E402  (sets up the import paths)
from checks import Check, compare     # noqa: E402
from law_config import load_law, describe   # noqa: E402

MODULES = ['t_machinery', 't_galaxies', 't_clusters', 't_lensing', 't_milky_way', 't_dwarfs', 't_precision', 't_collisions']
MARK = {'pass': 'PASS', 'close': 'close', 'fail': 'FAIL', 'info': 'info', 'error': 'ERROR'}
MOVE = dict(regressed='REGRESSED', improved='improved', worse='worse', better='better', changed='changed', same='', new='new', known='known', error='CRASH')


def fmt(v):
    if v is None: return '-'
    a = abs(v)
    if a == 0: return '0'
    if a >= 1e4 or a < 1e-3: return f'{v:.3e}'
    if a >= 100: return f'{v:.1f}'
    return f'{v:.4g}'


def git_commit():
    try:
        return subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=HERE, capture_output=True, text=True).stdout.strip()
    except Exception:
        return ''


def run_all(law, ctx, only=None):
    checks, seconds = [], {}
    for name in MODULES:
        mod = importlib.import_module(name)
        if only and mod.GROUP not in only and name[2:] not in only:
            continue
        t0 = time.monotonic()
        print(f'== {mod.GROUP}', flush=True)
        try:
            got = mod.run(law, ctx)
        except Exception as e:
            tb = traceback.format_exc()
            print(tb, flush=True)
            got = [Check(id=f'{mod.GROUP}.crash', group=mod.GROUP, title=f'{name} crashed', value=None, status='error',
                         note=f'{type(e).__name__}: {e}', detail=dict(traceback=tb[-2000:]))]
        seconds[mod.GROUP] = round(time.monotonic() - t0, 1)
        checks += got
    return checks, seconds


def scoreboard(rows, base):
    lines = []
    group = None
    for r in rows:
        if r['group'] != group:
            group = r['group']; lines.append(f'\n[{group}]')
        b = base.get(r['id']) if base else None
        was = f"was {fmt(b['value'])} {b['status']}" if b else ''
        mv = MOVE.get(r['move'], r['move'])
        st = MARK.get(r['status'], r['status'])
        lines.append(f"  {st:5s} {r['title'][:78]:78s} {fmt(r['value']):>11s} {r['unit']:9s} | {r['target'][:60]:60s} {mv:9s} {was}")
    return '\n'.join(lines)


def report_md(payload):
    L = [f"# Regression suite: {payload['law']['name']} ({payload['tier']} tier)", '',
         f"Law: {describe(payload['law'])}", '', f"{payload['law']['description']}", '',
         f"Run {payload['date']} on commit {payload['commit']}, {payload['seconds']:.0f} s. "
         f"Baseline: {payload['baseline'] or 'none'}.", '',
         '| | count |', '|---|---:|'] + [f'| {k} | {v} |' for k, v in payload['summary']['status'].items()] + \
        [f"| vs baseline: {k} | {v} |" for k, v in payload['summary']['movement'].items() if k != 'same'] + ['']
    group = None
    for r in payload['checks']:
        if r['group'] != group:
            group = r['group']
            L += ['', f'## {group}', '', '| Check | Ours | Measured / required | Status | vs baseline |', '|---|---:|---|---|---|']
        b = r.get('baseline')
        was = (f"{MOVE.get(r['move'], r['move'])} (was {fmt(b['value'])}, {b['status']})" if b else MOVE.get(r['move'], r['move'])).strip()
        L.append(f"| {r['title']} | {fmt(r['value'])} {r['unit']} | {r['target']} | {r['status']} | {was} |")
    return '\n'.join(L) + '\n'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--law', default=None, help='candidate name in candidates/ or a JSON path (default: the round-3 law)')
    ap.add_argument('--tier', choices=('quick', 'full'), default='quick')
    ap.add_argument('--only', default='', help='comma-separated groups (machinery, galaxies, clusters, lensing, milky_way, dwarfs, precision, collisions)')
    ap.add_argument('--baseline', type=Path, default=HERE / 'baseline.json')
    ap.add_argument('--save-baseline', action='store_true', help='write this run as the baseline')
    ap.add_argument('--output-dir', type=Path, default=None, help='default: runs/<law>-<tier>/ in this folder')
    ap.add_argument('--no-fail', action='store_true', help='exit 0 even when something regressed')
    args = ap.parse_args()
    t0 = time.monotonic()
    law = load_law(args.law)
    ctx = C.Context(tier=args.tier)
    print(f"law: {law['name']} -- {describe(law)}", flush=True)
    if law['refit']:
        print(f"refitting {', '.join(law['refit'])} on the home data...", flush=True)
        law = C.refit_constants(law, ctx)
        print(f"  refitted: {describe(law)}", flush=True)
    only = {s.strip() for s in args.only.split(',') if s.strip()} or None
    checks, seconds = run_all(law, ctx, only)

    base_all = json.loads(args.baseline.read_text()) if args.baseline.exists() else None
    base = {c['id']: c for c in base_all['checks']} if base_all else {}
    rows = []
    for c in checks:
        d = c.as_dict()
        b = base.get(c.id)
        d['move'] = compare(d, b) if base_all else ('error' if c.status == 'error' else '')
        if b: d['baseline'] = dict(value=b['value'], status=b['status'], score=b.get('score'))
        rows.append(d)
    count = lambda key, vals: {v: sum(1 for r in rows if r[key] == v) for v in vals}
    summary = dict(status=count('status', ('pass', 'close', 'fail', 'info', 'error')),
                   movement=count('move', ('regressed', 'improved', 'worse', 'better', 'changed', 'new', 'known', 'same', 'error')) if base_all else {},
                   not_run=sorted(set(base) - {r['id'] for r in rows}))
    payload = dict(suite='hot-companion regression suite', law=law, tier=args.tier, only=sorted(only) if only else None,
                   date=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC'), commit=git_commit(),
                   baseline=(f"{args.baseline.name} ({base_all['law']['name']}, {base_all['tier']}, {base_all['date']})" if base_all else None),
                   seconds=time.monotonic() - t0, seconds_by_group=seconds, summary=summary, checks=rows)
    print(scoreboard(rows, base), flush=True)
    print(f"\nstatus: {summary['status']}", flush=True)
    if base_all:
        print(f"against the baseline ({payload['baseline']}): {summary['movement']}", flush=True)
        for key in ('regressed', 'error', 'improved', 'worse', 'better', 'changed'):
            for r in rows:
                if r['move'] == key:
                    print(f"  {MOVE[key]:9s} {r['id']}: {fmt(r['value'])} ({r['status']}), was {fmt(r['baseline']['value']) if r.get('baseline') else '-'}"
                          f" ({r['baseline']['status'] if r.get('baseline') else '-'})", flush=True)
    out = args.output_dir or HERE / 'runs' / f"{law['name']}-{args.tier}"
    out.mkdir(parents=True, exist_ok=True)
    (out / 'results.json').write_text(json.dumps(payload, indent=1, default=float) + '\n')
    (out / 'report.md').write_text(report_md(payload))
    print(f"\nwrote {out / 'results.json'} and report.md ({payload['seconds']:.0f} s)", flush=True)
    if args.save_baseline:
        keep = [dict(id=r['id'], group=r['group'], title=r['title'], value=r['value'], unit=r['unit'], target=r['target'],
                     status=r['status'], score=r['score']) for r in rows]
        if only and base_all:        # a partial run updates only its own checks
            ids = {r['id'] for r in keep}
            keep = [c for c in base_all['checks'] if c['id'] not in ids] + keep
        args.baseline.write_text(json.dumps(dict(suite=payload['suite'], law=law, tier=args.tier, date=payload['date'], commit=payload['commit'],
                                                 checks=keep), indent=1, default=float) + '\n')
        print(f'saved the baseline: {args.baseline}', flush=True)
    bad = summary['movement'].get('regressed', 0) + summary['status'].get('error', 0) if base_all else summary['status'].get('error', 0)
    sys.exit(1 if bad and not args.no_fail else 0)


if __name__ == '__main__':
    main()

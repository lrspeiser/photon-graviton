"""CC-2 stage 2B-F1 revision (see revision-protocol.md): a direct normalization search around each combination's best
sampled RMSE, repeats at higher resolution, F3 with more tracers, and universality if the best source changes.

    python revision.py [--canonical] [--output-dir DIR]

It reads the canonical 2B-F1 results (f1-results.json, or the file F1R_SOURCE names) for every combination's
completed runs. F1_WORKERS sets the pool size and F1_BUDGET each run's wall-clock budget; F1_SMOKE=1 runs against a
smoke output. Seeds come from the revision's own counter and are drawn in a fixed order, and results return in the
order the tasks were issued.
"""
import json
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import f1 as F  # noqa: E402  (also puts field, mc, formation and evidence_io on the path)
import evidence_io  # noqa: E402

SEEDS, HI_SEEDS, HI_FACTOR, F3_FACTOR = 3, 4, 4., 16
BRACKET = (-4, -2, 0, 2, 4)             # rates q_s 3^(e/8): five rates a step of 3^(1/4) apart
MAX_EXT = 2                             # outward extensions per side, one step of 3^(1/4) each
SCORES = ('rmse_38', 'rmse_inner20', 'slope_8_20', 'cost')
KEEP = ('status', 'mw', 'gates', 'M_T', 'B1_ratio', 'B2_ratio', 'growth_dlnM_dlnt_at_T', 'regime', 'energy',
        'test_radii', 'profile', 'ledger', 'mass_ledger_relative')
SOURCE = Path(os.environ.get('F1R_SOURCE', HERE/'f1-results.json'))
log = F.log


def targets(d):
    """Every combination with a completed run, and the 150 kpc sensitivity: the completed runs' rates and RMSE, the
    rates where a run stopped early, and the rescaling-root verification for comparison."""
    out = []
    for name, c in d['combos'].items():
        runs = [(r['q'], r['rmse_38']) for r in c['ladder'] if r['status'] == 'completed' and r.get('rmse_38') is not None]
        v = c.get('verify') or {}
        if v.get('status') == 'completed' and v.get('mw'):
            runs.append((c['q_star'], v['mw']['rmse_38']))
        if runs:
            out.append(dict(name=name, v_d=c['v_d_kms'], sm=c['sigma_over_m_cm2g'], gravity=c['bath_gravity'], R_b=None,
                            runs=runs, stopped=[r['q'] for r in c['ladder'] if r['status'] != 'completed'],
                            root=dict(q=c['q_star'], gates=v.get('gates'),
                                      **({k: v['mw'][k] for k in SCORES} if v.get('mw') else {}))))
    S = d.get('sensitivity_150kpc') or {}
    if d.get('best') and S.get('ladder'):
        b = d['combos'][d['best']['name']]
        runs = [(r['q'], r['rmse_38']) for r in S['ladder'] if r['status'] == 'completed' and r.get('rmse_38') is not None]
        if runs:
            sv = S.get('verify') or {}
            # the archive keeps the 150 kpc verification's scores but not its rate, so it cannot be a best sample
            out.append(dict(name='150 kpc zone | ' + d['best']['name'], v_d=b['v_d_kms'], sm=b['sigma_over_m_cm2g'],
                            gravity=b['bath_gravity'], R_b=150., runs=runs,
                            stopped=[r['q'] for r in S['ladder'] if r['status'] != 'completed'],
                            root=dict(gates=sv.get('gates'), **({k: sv['mw'][k] for k in SCORES} if sv.get('mw') else {}))))
    return out


def summary(rs, D):
    """Seed means and standard deviations of the scores, the gates on the means, and each seed's own gates."""
    vals = {k: [r['mw'][k] for r in rs] for k in SCORES}
    mean = {k: float(np.mean(v)) for k, v in vals.items()}
    sd = {k: (float(np.std(v, ddof=1)) if len(v) > 1 else None) for k, v in vals.items()}
    return dict(mean=mean, sd=sd, gates=F.gates(mean, D), seed_gates=[r['gates'] for r in rs])


def tier(g):
    return 0 if g['profile_gate_passed'] else 1 if (g['G1_fit'] and g['G2_shape']) else 2 if g['G1_fit'] else 3


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    D = F.mw_data()
    d = json.loads(SOURCE.read_text(encoding='utf-8'))
    T = targets(d)
    byname = {t['name']: t for t in T}
    cnt = [20000]

    def nxt():
        cnt[0] += 1
        return cnt[0]

    def rate(t, e):
        return t['q_s']*3**(e/8)

    def spec(t, e, role, n_factor=1., full=False):
        sp = dict(role=role, key='MW', v_d=t['v_d'], sm=t['sm'], gravity=t['gravity'], q=rate(t, e), rng=nxt(),
                  full=full, target=t['name'], step=e)
        if n_factor != 1.:
            sp['n_factor'] = n_factor
        if t['R_b']:
            sp['R_b'] = t['R_b']
        return sp

    def blocked(t, e):
        """At or beyond a rate where a run of this target stopped early, on the far side of q_s."""
        q = rate(t, e)
        return any((s >= t['q_s'] and q >= s) or (s <= t['q_s'] and q <= s) for s in t['stopped'])

    def complete(t):
        """Mean RMSE at every step where all seeds completed."""
        return {e: float(np.mean([r['mw']['rmse_38'] for r in rs])) for e, rs in t['evals'].items()
                if len(rs) >= SEEDS and all(r['status'] == 'completed' and r.get('mw') for r in rs)}

    def absorb(rs):
        for r in rs:
            if r['spec']['role'] in ('bracket', 'extension', 'refinement'):
                sp = r['spec']
                t = byname[sp['target']]
                t['evals'].setdefault(sp['step'], []).append(r)
                if r['status'] != 'completed':
                    t['stopped'].append(sp['q'])

    for t in T:
        t['q_s'], t['best_rmse'] = min(t['runs'], key=lambda p: p[1])
        t['evals'], t['skipped'] = {}, []
        log(f"target {t['name']}: best sampled RMSE {t['best_rmse']:.1f} km/s at q = {t['q_s']:.4g}")
    # F3 at higher resolution runs in a process of its own from the start and is collected at the end, so that its long
    # run does not hold up the rounds; a failure is recorded, not raised, so that the rest of the results survive it.
    # Its cap is four times its births, so the engine never thins it, and it has a budget of its own (amendment A2)
    n3 = int(64000*(.05 if F.SMOKE else 1.)*F3_FACTOR)
    f3spec = dict(role='validation', name='F3', key='MW', sm=0., q=1e3, gravity=False, n=n3, n_max=4*n3, rng=nxt(),
                  budget=float(os.environ.get('F3_BUDGET', '900' if F.SMOKE else '43200')))
    f3pool = Pool(1)
    f3job = f3pool.apply_async(F.task, (f3spec,))
    log(f'F3 with {n3} tracers started in its own process')
    # round A: the bracket
    specs = []
    for t in T:
        for e in BRACKET:
            if blocked(t, e):
                t['skipped'].append(e)
            else:
                specs += [spec(t, e, 'bracket') for _ in range(SEEDS)]
    absorb(F.run_tasks(specs, 'round A (bracket)'))
    # round B: while the lowest mean sits at an end of the completed steps, step outward (at most MAX_EXT per side)
    for ext in range(MAX_EXT):
        specs = []
        for t in T:
            m = complete(t)
            if not m:
                continue
            e_best = min(m, key=m.get)
            for side, edge in ((1, max(m)), (-1, min(m))):
                e_new = edge + 2*side
                if e_best == edge and abs(e_new) <= max(BRACKET) + 2*MAX_EXT and e_new not in t['evals']:
                    if blocked(t, e_new):
                        t['skipped'].append(e_new)
                    else:
                        specs += [spec(t, e_new, 'extension') for _ in range(SEEDS)]
        if not specs:
            break
        absorb(F.run_tasks(specs, f'round B.{ext + 1} (extension)'))
    # round C: the refinement, a step of 3^(1/8) either side of the lowest mean
    specs = []
    for t in T:
        m = complete(t)
        if not m:
            continue
        e_best = min(m, key=m.get)
        for e_new in (e_best - 1, e_best + 1):
            if e_new in t['evals']:
                continue
            if blocked(t, e_new):
                t['skipped'].append(e_new)
            else:
                specs += [spec(t, e_new, 'refinement') for _ in range(SEEDS)]
    if specs:
        absorb(F.run_tasks(specs, 'round C (refinement)'))
    for t in T:
        m = complete(t)
        t['sel'] = min(m, key=m.get) if m else None
        if t['sel'] is not None:
            s = summary(t['evals'][t['sel']], D)
            log(f"selected {t['name']}: q = {rate(t, t['sel']):.4g}, mean RMSE {s['mean']['rmse_38']:.1f}, slope "
                f"{s['mean']['slope_8_20']:.2f}, cost {s['mean']['cost']:.3g}, gates {s['gates']}")
    # round D: four more seeds at four times the tracers wherever the seed means pass G1 and G2
    specs = []
    for t in T:
        if t['sel'] is not None:
            g = summary(t['evals'][t['sel']], D)['gates']
            if g['G1_fit'] and g['G2_shape']:
                specs += [spec(t, t['sel'], 'high resolution', n_factor=HI_FACTOR, full=True) for _ in range(HI_SEEDS)]
    hi = {}
    for r in (F.run_tasks(specs, 'round D (four times the tracers)') if specs else []):
        hi.setdefault(r['spec']['target'], []).append(r)
    # round R7 (amendment A3): numerical controls beyond tracer count, for the reference case at its selected rate, one
    # change at a time: the timestep halved, the potential grid doubled, and the pools regenerated twice as often with
    # twice the draws. Round D measures tracer convergence only
    controls = {}
    ref = byname.get(f'v_d {F.V_D[0]:g} km/s | 0 cm2/g | bath gravity omitted')
    if ref is not None and ref['sel'] is not None:
        variants = dict(timestep_halved=dict(dmax_factor=.5), grid_doubled=dict(n_grid=2048),
                        pools_doubled=dict(regen=5, pool_factor=2.))
        specs = [dict(spec(ref, ref['sel'], 'numerical control'), control=name, **kw)
                 for name, kw in variants.items() for _ in range(SEEDS)]
        for r in F.run_tasks(specs, 'round R7 (numerical controls for the reference case)'):
            controls.setdefault(r['spec']['control'], []).append(r)
    # the revised best source; universality reruns only if it changed
    cands = [(t, summary(t['evals'][t['sel']], D)) for t in T if t['sel'] is not None and not t['R_b']]
    best_t, best_s = min(cands, key=lambda ts: (tier(ts[1]['gates']),
                                                ts[1]['mean']['cost'] if tier(ts[1]['gates']) < 3 else ts[1]['mean']['rmse_38']))
    changed = bool(d.get('best') and best_t['name'] != d['best']['name'])
    uni = None
    if changed:
        base = dict(v_d=best_t['v_d'], sm=best_t['sm'], gravity=best_t['gravity'], q=rate(best_t, best_t['sel']), full=True)
        uspecs = [dict(base, role='universality', key=k, rng=nxt()) for k in ('J1630', 'Coma low', 'Coma high')]
        uni = {r['spec']['key']: {k: x for k, x in r.items() if k != 'spec'}
               for r in F.run_tasks(uspecs, 'round E (universality)')}
    f3 = f3job.get()
    f3pool.close()
    f3pool.join()
    if f3['status'] == 'completed':
        log(f"F3 with {n3} tracers: passed={f3.get('passed')}, worst bin {f3.get('worst_over_tolerance')} of its tolerance")
    else:
        log(f"F3 with {n3} tracers did not complete: {f3['status'][:200]}")
    # output
    out = dict(scope='CC-2 stage 2B-F1 revision: a direct normalization search, repeats at higher resolution, and F3 with '
                     'more tracers', smoke=F.SMOKE, source=SOURCE.name, T_Gyr=F.T_REF, seeds_per_rate=SEEDS,
               high_resolution_seeds=HI_SEEDS, high_resolution_factor=HI_FACTOR, targets={})
    for t in T:
        rates = []
        for e in sorted(t['evals']):
            rs = t['evals'][e]
            ok = all(r['status'] == 'completed' and r.get('mw') for r in rs)
            rates.append(dict(step_eighths=e, q=rate(t, e), role=rs[0]['spec']['role'],
                              seeds=[dict(rng=r['spec']['rng'], status=r['status'], gates=r.get('gates'),
                                          **({k: r['mw'][k] for k in SCORES} if r.get('mw') else {})) for r in rs],
                              **(summary(rs, D) if ok else {})))
        entry = dict(v_d_kms=t['v_d'], sigma_over_m_cm2g=t['sm'], bath_gravity=t['gravity'], R_b_kpc=t['R_b'],
                     best_sampled=dict(q=t['q_s'], rmse_38=t['best_rmse']), root_verification=t['root'],
                     skipped_steps=sorted(set(t['skipped'])), rates=rates,
                     q_sel=(rate(t, t['sel']) if t['sel'] is not None else None), sel_step_eighths=t['sel'],
                     selected=(summary(t['evals'][t['sel']], D) if t['sel'] is not None else None))
        if t['name'] in hi:
            rs = hi[t['name']]
            ok = [r for r in rs if r['status'] == 'completed' and r.get('mw')]
            h = dict(completed=len(ok), of=len(rs), robust_pass=False,
                     runs=[dict({k: r.get(k) for k in KEEP}, rng=r['spec']['rng']) for r in rs])
            if len(ok) == len(rs):
                s = summary(ok, D)
                n_all = sum(bool(g['profile_gate_passed']) for g in s['seed_gates'])
                h.update(s, seeds_passing_all_three=n_all, robust_pass=bool(s['gates']['profile_gate_passed'] and n_all >= 3))
            entry['high_resolution'] = h
        out['targets'][t['name']] = entry
    out['F3_high_resolution'] = dict({k: v for k, v in f3.items() if k not in ('spec', 'traceback')}, n=n3)
    # R7: each control's seed means against the canonical-size seed means at the same rate, within two standard errors
    # of the difference or 1 km/s (RMSE) and 0.05 (slope), whichever is larger
    nc = None
    if controls:
        base_runs = ref['evals'][ref['sel']]
        base = summary(base_runs, D)
        nc = dict(target=ref['name'], q=rate(ref, ref['sel']), canonical=dict(mean=base['mean'], sd=base['sd']), variants={})
        for name, rs in controls.items():
            ok = [r for r in rs if r['status'] == 'completed' and r.get('mw')]
            if len(ok) < len(rs):
                nc['variants'][name] = dict(completed=len(ok), of=len(rs), statuses=[r['status'] for r in rs])
                continue
            s = summary(ok, D)
            agrees = {}
            for k, floor in (('rmse_38', 1.), ('slope_8_20', .05)):
                se = float(np.sqrt((s['sd'][k] or 0.)**2/len(ok) + (base['sd'][k] or 0.)**2/len(base_runs)))
                agrees[k] = bool(abs(s['mean'][k] - base['mean'][k]) <= max(2*se, floor))
            nc['variants'][name] = dict(mean=s['mean'], sd=s['sd'], agrees=agrees,
                                        runs=[dict({k: r.get(k) for k in KEEP}, rng=r['spec']['rng']) for r in rs])
            log(f"numerical control {name}: RMSE {s['mean']['rmse_38']:.2f} against {base['mean']['rmse_38']:.2f}, slope "
                f"{s['mean']['slope_8_20']:.2f} against {base['mean']['slope_8_20']:.2f}, agrees {agrees}")
    out['numerical_controls'] = nc
    out['best'] = dict(name=best_t['name'], tier=tier(best_s['gates']), changed=changed)
    out['universality'] = uni
    out['runtime_seconds'] = time.time() - F.T0
    text = json.dumps(out, indent=1, default=float)
    if F.SMOKE:
        od = evidence_io.output_dir(args, 'companion-source-revision-smoke')
        (od/'revision-smoke.json').write_text(text, encoding='utf-8')
        log(f'smoke output {od}')
        return 0
    return evidence_io.finish(args, 'companion-source-revision', text, HERE/'revision-results.json',
                              ignore={'/runtime_seconds'}, rules=((r'.*/(wall_s|runtime_s)$', None, None),))


if __name__ == '__main__':
    sys.exit(main())

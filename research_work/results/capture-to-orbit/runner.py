"""Frozen comparison runner for the capture-to-orbit feasibility test (RB-1).

    python runner.py [--output-dir DIR] [--canonical] [--limit N] [--workers N]

Results go to a fresh directory under research_work/generated/ by default;
--canonical additionally rewrites results.json and the figure beside this
file. Nothing is fitted or selected: every variant was declared in protocol.md.
"""
import argparse
import hashlib
import json
import os
import sys
import time
import uuid
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import inputs as I              # noqa: E402
import kinematics as kin        # noqa: E402
import ledger as lg             # noqa: E402
import supported_profile as sp  # noqa: E402

G = sp.G
R_EVAL = np.geomspace(1e-3, 1e5, 481)            # model-defined profile radii (kpc)
LR_EVAL = np.log(R_EVAL)
BASE = dict(n_sites=48, n_speed=24, n_angle=12, n_theta=96)
FINE = dict(n_sites=96, n_speed=48, n_angle=24, n_theta=192)
PLUNGE_RD = .05
VARIANTS = ['capture_primary', 'capture_baryon_potential', 'capture_self_consistent']
MODELS = ['baryons', 'MOND_raw', 'reference', 'mond_guided'] + VARIANTS


def fraction_at(frac, r):
    """Monotone interpolation of an enclosed-mass fraction from the model grid."""
    return PchipInterpolator(LR_EVAL, frac)(np.clip(np.log(r), LR_EVAL[0], LR_EVAL[-1]))


def production_sites(rec, field, n_sites):
    """Log-spaced sites covering 99.98% of the absorption rate rho_b J."""
    lr = np.log(I.GRID)
    dm = rec['disk_dm_dlnr'] + rec['bulge_dm_dlnr']
    W = dm*field['J']

    def cum(f):
        return cumulative_trapezoid(f, lr, initial=0)

    cW = cum(W)
    lo, hi = np.interp([1e-4*cW[-1], (1 - 1e-4)*cW[-1]], cW, lr)
    edges = np.linspace(lo, hi, n_sites + 1)

    def bins(c):
        return np.diff(np.interp(edges, lr, c))

    mass, disk, absorb, zw = bins(cum(dm)), bins(cum(rec['disk_dm_dlnr'])), bins(cW), bins(cum(W*field['zeta']))
    return dict(r=np.exp((edges[1:] + edges[:-1])/2), mass=mass, W=absorb, J=absorb/mass, zeta=zw/absorb,
                disk_fraction=np.clip(disk/mass, 0, 1), covered=float(absorb.sum()/cW[-1]))


def population(pot, sites, q, r_eval, bins=None):
    inj = sp.uniform_ball_injections(pot, sites['r'], sites['W'], q['n_speed'], q['n_angle'])
    prof = sp.population_profile(pot, inj['r'], inj['E'], inj['L'], inj['weight'], r_eval, q['n_theta'],
                                 bin_edges=bins)
    return prof, inj, prof['mass']/inj['weight'].sum()


def self_consistent(rec, sites, M_ret, frac, q, tol=1e-4, damping=.5, maxit=40):
    """Stationary fixed point: the population's own gravity replaces the reference reservoir."""
    change, new = np.inf, frac
    for it in range(1, maxit + 1):
        pot = sp.SphericalPotential(I.GRID, rec['mass'] + M_ret*fraction_at(frac, I.GRID))
        _, _, new = population(pot, sites, q, R_EVAL)
        change = float(np.max(np.abs(new - frac)))
        if change < tol:
            return new, dict(iterations=it, max_fraction_change=change, converged=True)
        frac = damping*new + (1 - damping)*frac
    return new, dict(iterations=maxit, max_fraction_change=change, converged=False)


def orbit_diagnostics(prof, inj, rd):
    w = inj['weight']/inj['weight'].sum()
    per = np.nan_to_num(prof['period_gyr'], nan=0.)
    circ = prof['circularity']
    summ = sp.summarize(R_EVAL, prof['mass'], [10*rd, 100*rd, 1e5])
    beta = [dict(r_over_rd=float(np.sqrt(a*b)/rd), mass_fraction=float(m/prof['bin_mass'].sum()), beta=float(v))
            for a, b, m, v in zip(prof['bin_edges'][:-1], prof['bin_edges'][1:], prof['bin_mass'], prof['beta'])]
    return dict(r50_kpc=summ['r50'], r90_kpc=summ['r90'], outer_density_slope=summ['outer_density_slope'],
                fraction_outside_10rd=summ['fraction_outside'][str(10*rd)],
                fraction_outside_100rd=summ['fraction_outside'][str(100*rd)],
                fraction_outside_1e5kpc=summ['fraction_outside'][str(1e5)],
                plunging_fraction=float(w[prof['rp'] < PLUNGE_RD*rd].sum()),
                period_over_10gyr_fraction=float(w[per > 10].sum()), period_over_100gyr_fraction=float(w[per > 100].sum()),
                mean_circularity=float(np.sum(w*circ)), circularity_below_0p3=float(w[circ < .3].sum()),
                circularity_above_0p9=float(w[circ > .9].sum()), anisotropy=beta)


def speeds(vb2, r, M_ret, frac):
    return np.sqrt(vb2 + G*M_ret*fraction_at(frac, r)/r)


def run_system(s):
    """One galaxy or Milky Way scenario: receivers, potentials, populations, ledger."""
    rd, rec, field, r, vb2, M_ret = (s[k] for k in ('rd', 'rec', 'field', 'r', 'vb2', 'M_ret'))
    sites = production_sites(rec, field, BASE['n_sites'])
    primary = sp.SphericalPotential(I.GRID, rec['mass'] + field['mass'])
    prof, inj, frac = population(primary, sites, BASE, np.concatenate([R_EVAL, r]), rd*np.geomspace(.25, 64, 9))
    n = len(R_EVAL)
    direct, frac = frac[n:], frac[:n]
    prof = dict(prof, mass=prof['mass'][:n], r=R_EVAL)
    fracs = dict(capture_primary=frac)
    _, _, fracs['capture_baryon_potential'] = population(sp.SphericalPotential(I.GRID, rec['mass']), sites, BASE, R_EVAL)
    fracs['capture_self_consistent'], sc = self_consistent(rec, sites, M_ret, frac, BASE)
    _, _, fine = population(primary, production_sites(rec, field, FINE['n_sites']), FINE, R_EVAL)
    _, _, subset = population(primary, sites, BASE, r[::2])
    pred = {k: speeds(vb2, r, M_ret, f) for k, f in fracs.items()}
    ref_model = np.sqrt(vb2 + G*np.interp(r, I.GRID, field['mass'])/r)
    checks = dict(
        refinement_max_speed_change_kms=float(np.max(np.abs(speeds(vb2, r, M_ret, fine) - pred['capture_primary']))),
        sampling_interpolation_max_speed_difference_kms=float(np.max(np.abs(np.sqrt(vb2 + G*M_ret*direct/r)
                                                                            - pred['capture_primary']))),
        sampling_subset_max_fraction_difference=float(np.max(np.abs(subset - direct[::2]))),
        site_coverage=sites['covered'], self_consistent=sc,
        reference_model_speed_vs_archived_kms=float(np.max(np.abs(ref_model - s['reference']))))
    book = lg.ledger(primary, sites['r'], sites['mass'], sites['J'], sites['zeta'], sites['disk_fraction'], M_ret)
    cum_rec = rec['mass']/rec['mass'][-1]
    shape = dict(receiver_r50_kpc=float(np.interp(.5, cum_rec, I.GRID)),
                 reference_r50_kpc=float(np.interp(.5, field['mass']/field['mass'][-1], I.GRID)),
                 population_fraction_inside_last_radius=float(fraction_at(frac, r[-1])),
                 reference_fraction_inside_last_radius=float(np.interp(r[-1], I.GRID, field['mass'])/field['mass'][-1]),
                 receiver_mass_Msun=float(rec['mass'][-1]), inventory_over_receivers=float(M_ret/rec['mass'][-1]))
    profiles = {k: np.round(v[::8], 10).tolist() for k, v in fracs.items()}
    profiles.update(r_kpc=R_EVAL[::8].tolist(),
                    reference=np.round(np.interp(R_EVAL[::8], I.GRID, field['mass'])/field['mass'][-1], 10).tolist(),
                    receivers=np.round(np.interp(R_EVAL[::8], I.GRID, cum_rec), 10).tolist())
    return dict(predictions={k: v.tolist() for k, v in pred.items()}, orbits=orbit_diagnostics(prof, inj, rd),
                shape=shape, ledger=book, checks=checks, profile_fraction=profiles)


def sparc_job(d):
    c = d['catalog']
    rec = I.sparc_receivers(d['rotmod'], c)
    vb2 = d['reference']**2 - G*d['reference_enclosed']/d['r']
    out = run_system(dict(rd=c['rd'], rec=rec, field=I.reference_field(c['rd'], c['L9']), r=d['r'], vb2=vb2,
                          M_ret=d['total_inventory'], reference=d['reference']))
    out.update(galaxy=d['name'], split=d['split'], capped=bool(d['capped']), rd_kpc=c['rd'], gas_fit=rec['gas_fit'],
               R_kpc=d['r'].tolist(), observed_kms=d['y'].tolist(), baryons_kms=np.sqrt(vb2).tolist(),
               reference_kms=d['reference'].tolist(), mond_guided_kms=d['mond_guided'].tolist(),
               MOND_raw_kms=d['mond_raw'].tolist())
    return out


def milky_way_job(run):
    rec = I.milky_way_receivers(run['baryons'])
    out = run_system(dict(rd=run['rd'], rec=rec, field=I.reference_field(run['rd'], run['L9']), r=run['R'],
                          vb2=run['vb']**2, M_ret=run['total_inventory'], reference=run['reference']))
    out.update(baryons=run['baryons'], Rd_kpc=run['rd'], luminosity_proxy_factor=run['lf'], R_kpc=run['R'].tolist(),
               observed_kms=run['y'].tolist(), baryons_kms=run['vb'].tolist(), reference_kms=run['reference'].tolist(),
               mond_guided_kms=run['mond_guided'].tolist(), MOND_raw_kms=run['mond_raw'].tolist(),
               split=run['split'].tolist(), stellar_mass_check=[rec['stellar_mass'], rec['stellar_mass_nominal']])
    return out


def model_speeds(row, m):
    return np.array(row['predictions'][m] if m in VARIANTS else row[m + '_kms'])


def sparc_scores(rows):
    out = {}
    for m in MODELS:
        out[m] = {}
        for s in ['train', 'validation', 'test']:
            rr = [r for r in rows if r['split'] == s]
            if rr:
                y = [np.array(r['observed_kms']) for r in rr]
                e = [np.mean((model_speeds(r, m) - v)**2) for r, v in zip(rr, y)]
                lg10 = [np.mean(np.log10(model_speeds(r, m)/v)**2) for r, v in zip(rr, y)]
                out[m][s] = dict(n=len(rr), RMSE_kms=float(np.sqrt(np.mean(e))), log_RMS=float(np.sqrt(np.mean(lg10))))
    return out


def radial_bins(rows):
    out = []
    for m in MODELS:
        for i, label in enumerate(['inner', 'middle', 'outer']):
            errs = []
            for r in rows:
                k = np.digitize(np.array(r['R_kpc'])/r['rd_kpc'], [1, 3]) == i
                if k.any():
                    errs.append(np.mean((model_speeds(r, m) - np.array(r['observed_kms']))[k]))
            out.append(dict(model=m, bin=label, mean_equal_galaxy_error_kms=float(np.mean(errs))))
    return out


def mw_scores(run):
    out = {}
    y, split = np.array(run['observed_kms']), np.array(run['split'])
    for m in MODELS:
        e = model_speeds(run, m) - y
        out[m] = {k: dict(RMSE_kms=float(np.sqrt(np.mean(e[mask]**2))), bias_kms=float(np.mean(e[mask])))
                  for k, mask in [('all', np.ones(len(y), bool)), ('inner', split == 'inner'), ('outer', split == 'outer')]}
    return out


def distribution(values):
    v = np.asarray(values, float)
    return dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
                min=float(v.min()), max=float(v.max()))


def plot(result, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = {r['galaxy']: r for r in result['sparc_rows']}
    picks = [n for n in ['NGC2403', 'DDO154', 'UGC02953'] if n in rows]
    fid = [r for r in result['milky_way_runs'] if r['Rd_kpc'] == 2.6 and r['luminosity_proxy_factor'] == 1]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5), layout='constrained')
    style = dict(baryons=('ordinary matter', 'k:'), reference=('original exact-third', 'C0-'),
                 mond_guided=('MOND-guided mixture', 'C2--'), capture_primary=('capture-to-orbit (primary)', 'C3-'),
                 capture_self_consistent=('capture-to-orbit (self-consistent)', 'C1-.'))
    for a, row in zip(list(ax[0]) + [ax[1, 0], ax[1, 1]], [rows[n] for n in picks] + fid):
        R = row['R_kpc']
        a.plot(R, row['observed_kms'], 'ko', ms=3, label='observed')
        for m, (label, fmt) in style.items():
            a.plot(R, model_speeds(row, m), fmt, lw=1.3, label=label)
        a.set(xlabel='radius (kpc)', ylabel='circular speed (km/s)',
              title=row.get('galaxy') or f"Milky Way baseline {row['baryons']} (fiducial)")
        a.grid(alpha=.2)
    ax[0, 0].legend(fontsize=7)
    a = ax[1, 2]
    books = [r['ledger'] for r in result['sparc_rows']]
    a.hist(np.log10([b['receiver_angular_momentum_loss_ideal'] for b in books]), bins=25, alpha=.7,
           label='ideal efficiency, gray drag')
    a.hist(np.log10([b['receiver_angular_momentum_loss_declared'] for b in books]), bins=25, alpha=.7,
           label='declared spectrum and drag')
    a.axvline(0, color='k', lw=1)
    a.set(xlabel='log10(receiver angular momentum lost / held)', ylabel='SPARC galaxies',
          title='Momentum the receivers must supply')
    a.legend(fontsize=7)
    fig.suptitle('Receiver-assisted threshold production (RB-1): frozen predictions, fixed retained inventory')
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--canonical', action='store_true')
    ap.add_argument('--limit', type=int)
    ap.add_argument('--workers', type=int, default=max(1, min(8, (os.cpu_count() or 2) - 2)))
    args = ap.parse_args()
    out_dir = (args.output_dir or I.ROOT/'research_work/generated'/(time.strftime('%Y%m%d-%H%M%S') + '-capture-to-orbit-'
                                                                    + uuid.uuid4().hex[:6])).resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        ap.error('Output directory must be new or empty.')
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    galaxies = I.sparc_galaxies()
    runs = I.milky_way_runs()
    if args.limit:
        galaxies, runs = galaxies[:args.limit], [r for r in runs if r['rd'] == 2.6 and r['lf'] == 1]
    with ProcessPoolExecutor(args.workers) as pool:
        rows = list(pool.map(sparc_job, galaxies))
        print(f'SPARC done {time.time() - t0:.0f}s', flush=True)
        mw = list(pool.map(milky_way_job, runs))
    print(f'Milky Way done {time.time() - t0:.0f}s', flush=True)
    for run in mw:
        run['scores'] = mw_scores(run)
    fid = [r for r in mw if r['Rd_kpc'] == 2.6 and r['luminosity_proxy_factor'] == 1]
    books = [r['ledger'] for r in rows]
    everything = rows + mw
    result = dict(
        scope='Conditional fixed-potential feasibility of one declared local interaction (RB-1): receiver-assisted '
              's-wave threshold production of one fixed-mass species on ordinary baryons, smooth companion spectrum. '
              'Steady phase-mixed spherical orbit average; fixed retained inventory equal to the reference total. '
              'No parameter fitted or selected; exposed data. Not formation, self-consistent evolution or stability.',
        protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))},
        input_sha256=I.input_hashes(), quadrature=dict(base=BASE, fine=FINE), plunge_radius_over_rd=PLUNGE_RD,
        reference_spectrum=dict(rate=kin.REFERENCE_RATE, power=kin.REFERENCE_POWER, chi=kin.REFERENCE_CHI,
                                drag_coefficient=kin.drag_coefficient(kin.REFERENCE_CHI)),
        sparc=dict(n_galaxies=len(rows), n_radii=sum(len(r['R_kpc']) for r in rows), scores=sparc_scores(rows),
                   radial=radial_bins(rows), capped_scores=sparc_scores([r for r in rows if r['capped']]),
                   capped_galaxies=[r['galaxy'] for r in rows if r['capped']]),
        milky_way=dict(fiducial=[dict(baryons=r['baryons'], scores=r['scores'], ledger=r['ledger'], orbits=r['orbits'],
                                      shape=r['shape'], checks=r['checks']) for r in fid],
                       improved_vs_reference={v: int(sum(r['scores'][v]['all']['RMSE_kms']
                                                         < r['scores']['reference']['all']['RMSE_kms'] for r in mw))
                                              for v in VARIANTS}, n_runs=len(mw)),
        ledger_summary={k: distribution([b[k] for b in books]) for k in
                        ['supply_multiplier', 'receiver_angular_momentum_loss_ideal', 'receiver_angular_momentum_loss_declared',
                         'receiver_kinetic_drain_declared', 'push_history_gyr_median', 'bound_kinetic_at_injection']},
        ideal_loss_above_one=int(sum(b['receiver_angular_momentum_loss_ideal'] > 1 for b in books)),
        orbit_summary={k: distribution([r['orbits'][k] for r in rows]) for k in
                       ['plunging_fraction', 'period_over_10gyr_fraction', 'mean_circularity', 'outer_density_slope',
                        'fraction_outside_10rd', 'fraction_outside_100rd']},
        checks=dict(max_refinement_speed_change_kms=max(r['checks']['refinement_max_speed_change_kms'] for r in everything),
                    max_sampling_interpolation_kms=max(r['checks']['sampling_interpolation_max_speed_difference_kms']
                                                       for r in everything),
                    max_sampling_subset_fraction_difference=max(r['checks']['sampling_subset_max_fraction_difference']
                                                                for r in everything),
                    max_reference_model_vs_archived_kms=max(r['checks']['reference_model_speed_vs_archived_kms']
                                                            for r in everything),
                    min_site_coverage=min(r['checks']['site_coverage'] for r in everything),
                    self_consistent_all_converged=all(r['checks']['self_consistent']['converged'] for r in everything)),
        runtime_seconds=time.time() - t0, sparc_rows=rows, milky_way_runs=mw)
    text = json.dumps(result, indent=1) + '\n'
    (out_dir/'results.json').write_text(text, encoding='utf-8', newline='\n')
    plot(result, out_dir/'capture-to-orbit.png')
    if args.canonical:
        (HERE/'results.json').write_text(text, encoding='utf-8', newline='\n')
        plot(result, HERE/'capture-to-orbit.png')
    print(json.dumps(dict(sparc=result['sparc']['scores'], checks=result['checks'], ledger=result['ledger_summary'],
                          mw=result['milky_way']['improved_vs_reference']), indent=1))
    print(f'Results: {out_dir}')


if __name__ == '__main__':
    main()

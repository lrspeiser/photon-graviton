"""RB-1 consistency revision: one incident field per control, anisotropic drag, closed ledger.

    python revision.py [--output-dir DIR] [--canonical] [--limit N] [--workers N]

Re-evaluates the archived RB-1 candidate on its own sites, receivers and fixed
primary potentials, using the unified incident field of protocol-addendum.md.
Each spectrum control (threshold-cut S1, extends-below S2) is reported on its
own. Population shapes are reused only if the per-site bound-production factor
varies by less than 1% within a control.
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import incident as inc          # noqa: E402
import inputs as I              # noqa: E402
import kinematics as kin        # noqa: E402
import runner as R              # noqa: E402
import supported_profile as sp  # noqa: E402

C, GYR = kin.C_KMS, sp.GYR
CHI, POWER = kin.REFERENCE_CHI, kin.REFERENCE_POWER
CONTROLS = ['threshold_cut', 'extends_below']
PUSH_LIMIT, SHAPE_TOLERANCE = .1, .01


def evaluate(rec, field, M_ret, rd):
    """Closed ledger for both spectrum controls on one system."""
    sites = R.production_sites(rec, field, R.BASE['n_sites'])
    pot = sp.SphericalPotential(I.GRID, rec['mass'] + field['mass'])
    r, m, fd = sites['r'], sites['mass'], sites['disk_fraction']
    a = rd*I.CP['scale_to_disk']
    T = I.CP['k0_per_kpc']*a
    vc, ve = pot.vcirc(r), pot.vesc(r)
    beta, vesc = vc/C, ve/C
    mom = [inc.moments(lambda mu, x=x: inc.reference_attenuation(x, T, mu)) for x in r/a]
    uA = np.array([q['u'] for q in mom])
    p_perp = np.array([q['p_perp'] for q in mom])
    zeta = np.array([q['zeta'] for q in mom])
    kappa_disk = 1 + (1 + CHI)*p_perp                 # tangential motion, this site's pressure tensor
    kappa_iso = 1 + (1 + CHI)/3                       # direction-averaged (bulge) receivers
    kappa = fd*kappa_disk + (1 - fd)*kappa_iso
    power = m*POWER*uA                                # absorbed power per site, per unit spectral density
    Jdisk = float(np.sum(m*fd*r*vc))
    lever = np.sum(power*fd*kappa_disk*r*vc)          # receiver angular-momentum loss rate x c^2
    gray_lever = np.sum(power*fd*(1 + p_perp)*r*vc)
    g = vc*vc/r
    out = dict(sites=len(r), p_perp_range=[float(p_perp.min()), float(p_perp.max())],
               kappa_disk_range=[float(kappa_disk.min()), float(kappa_disk.max())],
               receiver_angular_momentum=Jdisk, inventory_over_receivers=float(M_ret/rec['mass'][-1]),
               receiver_angular_momentum_loss_ideal_gray=float(M_ret*gray_lever/power.sum()/Jdisk),
               receiver_angular_momentum_loss_ideal=float(M_ret*lever/power.sum()/Jdisk))
    for spec in CONTROLS:
        B = np.array([inc.bound_rate(inc.separable(spec, lambda mu, x=x: inc.reference_attenuation(x, T, mu)),
                                     [b, 0., 0.], v, n_u=48, n_mu=32, n_phi=32) for x, b, v in zip(r/a, beta, vesc)])
        bound = m*B                                   # bound products per site, same normalization
        factor = B/(uA*vesc**3/3)                     # relative to the archived uniform-ball weight
        wts = bound/bound.sum()
        keep = wts > 1e-3
        mean_factor = float(np.sum(wts*factor))
        rest = bound.sum()
        kinetic = float(np.sum(bound*.3*vesc**2)/rest)
        multiplier = float(power.sum()/(rest*(1 + kinetic)))
        receiver = -float(np.sum(power*kappa*beta**2)/(rest*(1 + kinetic)))
        T_min = 1/PUSH_LIMIT*np.abs(zeta)*uA*multiplier*M_ret*C/(g*np.sum(m*uA))*GYR
        order = np.argsort(r)
        cw = np.cumsum(power[order])/power.sum()
        core = order[(cw > .005) & (cw < .995)]
        out[spec] = dict(
            supply_multiplier=multiplier, bound_efficiency=1/multiplier,
            bound_factor_mean=mean_factor,
            bound_factor_max_deviation=float(np.max(np.abs(factor[keep]/mean_factor - 1))),
            ledger=dict(absorbed=multiplier, bound_rest_plus_kinetic=1 + kinetic, receiver_energy_change=receiver,
                        escaping=multiplier - (1 + kinetic) - receiver,
                        recoil_bound='<= q/2 of absorbed energy, q=E\'/(M_R c^2)'),
            # Escaping is the remainder, so this only catches arithmetic slips; the independent
            # exact-kinematics channel closure is checks.ledger_closure.
            closure_relative=float((multiplier - ((1 + kinetic) + receiver + (multiplier - (1 + kinetic) - receiver)))
                                   /multiplier),
            receiver_angular_momentum_loss=float(M_ret*lever/(rest*(1 + kinetic))/Jdisk),
            push_history_gyr_median=float(np.median(T_min[core])), push_history_gyr_max=float(np.max(T_min[core])))
    return out


def sparc_job(d):
    c = d['catalog']
    rec = I.sparc_receivers(d['rotmod'], c)
    out = evaluate(rec, I.reference_field(c['rd'], c['L9']), d['total_inventory'], c['rd'])
    out.update(galaxy=d['name'], split=d['split'], capped=bool(d['capped']))
    return out


def milky_way_job(run):
    out = evaluate(I.milky_way_receivers(run['baryons']), I.reference_field(run['rd'], run['L9']),
                   run['total_inventory'], run['rd'])
    out.update(baryons=run['baryons'], Rd_kpc=run['rd'], luminosity_proxy_factor=run['lf'])
    return out


def distribution(v):
    v = np.asarray(v, float)
    return dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)),
                min=float(v.min()), max=float(v.max()))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--canonical', action='store_true')
    ap.add_argument('--limit', type=int)
    ap.add_argument('--workers', type=int, default=max(1, min(8, (os.cpu_count() or 2) - 2)))
    args = ap.parse_args()
    out_dir = (args.output_dir or I.ROOT/'research_work/generated'/(time.strftime('%Y%m%d-%H%M%S') + '-rb1-revision-'
                                                                    + uuid.uuid4().hex[:6])).resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        ap.error('Output directory must be new or empty.')
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    galaxies, runs = I.sparc_galaxies(), I.milky_way_runs()
    if args.limit:
        galaxies, runs = galaxies[:args.limit], [r for r in runs if r['rd'] == 2.6 and r['lf'] == 1]
    with ProcessPoolExecutor(args.workers) as pool:
        rows = list(pool.map(sparc_job, galaxies))
        mw = list(pool.map(milky_way_job, runs))
    archived = {r['galaxy']: r['ledger'] for r in json.loads((HERE/'results.json').read_text())['sparc_rows']}
    summary = {}
    for spec in CONTROLS:
        s = dict(supply_multiplier=distribution([r[spec]['supply_multiplier'] for r in rows]),
                 multiplier_over_archived=distribution([r[spec]['supply_multiplier']/archived[r['galaxy']]['supply_multiplier']
                                                        for r in rows]),
                 receiver_energy_change=distribution([r[spec]['ledger']['receiver_energy_change'] for r in rows]),
                 receiver_angular_momentum_loss=distribution([r[spec]['receiver_angular_momentum_loss'] for r in rows]),
                 bound_factor_mean=distribution([r[spec]['bound_factor_mean'] for r in rows]),
                 bound_factor_max_deviation=float(max(r[spec]['bound_factor_max_deviation'] for r in rows + mw)),
                 push_history_gyr_median=distribution([r[spec]['push_history_gyr_median'] for r in rows]),
                 max_relative_closure=float(max(abs(r[spec]['closure_relative']) for r in rows + mw)))
        s['population_shape_reused'] = s['bound_factor_max_deviation'] < SHAPE_TOLERANCE
        summary[spec] = s
    loss = [r['receiver_angular_momentum_loss_ideal'] for r in rows]
    result = dict(
        scope='Consistency revision of archived RB-1: unified incident field per spectrum control, per-site '
              'anisotropic drag from the attenuated pressure tensor, closed energy ledger. Same sites, receivers, '
              'fixed primary potentials and retained inventory; nothing fitted.',
        protocol_addendum_sha256=hashlib.sha256((HERE/'protocol-addendum.md').read_bytes()).hexdigest(),
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))},
        spectrum_constants=inc.spectrum_constants(), summary=summary,
        receiver_angular_momentum_loss_ideal=distribution(loss),
        receiver_angular_momentum_loss_ideal_gray=distribution([r['receiver_angular_momentum_loss_ideal_gray'] for r in rows]),
        ideal_loss_above_one=int(sum(v > 1 for v in loss)), n_galaxies=len(rows),
        archived_ideal_loss=distribution([archived[r['galaxy']]['receiver_angular_momentum_loss_ideal'] for r in rows]),
        p_perp_range=[min(r['p_perp_range'][0] for r in rows + mw), max(r['p_perp_range'][1] for r in rows + mw)],
        milky_way=mw, sparc_rows=rows, runtime_seconds=time.time() - t0)
    text = json.dumps(result, indent=1) + '\n'
    (out_dir/'revision-results.json').write_text(text, encoding='utf-8', newline='\n')
    if args.canonical:
        (HERE/'revision-results.json').write_text(text, encoding='utf-8', newline='\n')
    print(json.dumps(dict(summary=summary, ideal=result['receiver_angular_momentum_loss_ideal'],
                          ideal_gray=result['receiver_angular_momentum_loss_ideal_gray'],
                          above_one=result['ideal_loss_above_one'], p_perp=result['p_perp_range'],
                          runtime=result['runtime_seconds']), indent=1))
    print(f'Results: {out_dir}')


if __name__ == '__main__':
    main()

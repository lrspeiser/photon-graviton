"""Self-illumination pilot driver: controls, rigid ring, thin disk and galaxy-sample floor.

    python pilot.py [--output-dir DIR] [--canonical] [--workers N]

Results go to a fresh directory under research_work/generated/ by default;
--canonical also writes pilot-results.json here. Protocol: protocol.md.
Exits nonzero if any control fails.
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
import selfillum as si          # noqa: E402  (also puts capture-to-orbit on the path)
import incident as inc          # noqa: E402
import inputs as I              # noqa: E402
import kinematics as kin        # noqa: E402
import runner as R              # noqa: E402
import supported_profile as sp  # noqa: E402

C = kin.C_KMS
SPECS = ['extends_below', 'threshold_cut']


def receiver_report(pos, vel, w, x_r, beta_r, v_esc, spec):
    n, D, fw = si.beams(pos, vel, w, x_r)
    ch = si.channels(n, D, fw, beta_r, spec)
    b = np.asarray(beta_r, float)
    kappa = -(ch['dR'][1:]@(b/np.linalg.norm(b)))/(np.linalg.norm(b)*ch['P_abs'])
    eff = (si.bound_rate(n, D, fw, beta_r, v_esc, spec)/ch['P_abs'])/((v_esc**3/3)/kin.REFERENCE_POWER)
    lz = lambda p: float(np.cross(np.asarray(x_r, float), p)[2]/ch['P_abs'])
    return dict(kappa=float(kappa), bound_efficiency_factor=float(eff),
                receiver_energy_change_per_absorbed=float(ch['dR'][0]/ch['P_abs']),
                receiver_Lz_change_per_absorbed=lz(ch['dR'][1:]), products_Lz_per_absorbed=lz(ch['X'][1:]),
                incident_Lz_per_absorbed=lz(ch['K'][1:]), residual=ch['residual'])


def emitter_Lz_per_emitted(pos, vel, w, spec='extends_below'):
    P, dp = si.emitter_loss(vel, w, spec)
    return float(np.sum(np.cross(pos, dp)[:, 2])/P.sum())


def controls():
    b, out, ok = 1e-3, {}, True
    em = si.emission_check()
    ok &= all(abs(v - 1) < 1e-10 for v in em.values())
    out['emission'] = em
    for spec in SPECS:
        for wr in [0., .5, 1.]:
            pos, vel, w = si.shell(velocity=(wr*b, 0, 0))
            r = receiver_report(pos, vel, w, [0, 0, 0], [b, 0, 0], 3*b, spec)
            want = si.steady_flow_kappa(wr)
            ok &= abs(r['kappa'] - want) < (1e-4 if spec == 'extends_below' else 5e-3) and r['residual'] < 1e-8
            out[f'steady_shell/{spec}/w={wr}'] = dict(kappa=r['kappa'], first_order=want,
                                                       bound_efficiency_factor=r['bound_efficiency_factor'])
    out['passed'] = bool(ok)
    return out


def ring_runs():
    out, ok = {}, True
    for beta in [.01, 200/C]:
        pos, vel, w = si.ring(1., beta, n=4096)
        eL = emitter_Lz_per_emitted(pos, vel, w)
        for label, rr in [('interior_0.5R', .5), ('near_ring_0.95R', .95), ('exterior_1.5R', 1.5)]:
            for spec in SPECS:
                r = receiver_report(pos, vel, w, [rr, 0, 0], [0, rr*beta, 0], 3*beta, spec)
                r.update(emitter_Lz_per_emitted=eL)
                if rr < 1:                          # rigid co-rotation inside the ring: no net torque
                    ok &= abs(r['kappa']) < 1e-6
                out[f'beta={beta:.4g}/{label}/{spec}'] = r
    out['interior_zero_torque_check'] = bool(ok)
    return out


def disk_runs():
    out, beta = {}, 200/C
    for h in [.05, .1, .2]:
        pos, vel, w = si.disk(1., beta, 12., h, nR=64, nth=192, nz=6)
        eL = emitter_Lz_per_emitted(pos, vel, w)
        for Rr in ([.5, 1., 2., 4., 8.] if h == .1 else [1., 4.]):
            for spec in (SPECS if h == .1 else ['extends_below']):
                r = receiver_report(pos, vel, w, [Rr, 0, 0], [0, beta, 0], 2.5*beta, spec)
                r.update(emitter_Lz_per_emitted=eL, emitter_over_receiver_specific_Lz=eL/(Rr*beta))
                out[f'h={h}/R={Rr}/{spec}'] = r
    return out


def floor(rec, field, M_ret, rd):
    """Conservation floor: combined baryonic loss when every companion is absorbed internally."""
    sites = R.production_sites(rec, field, R.BASE['n_sites'])
    pot = sp.SphericalPotential(I.GRID, rec['mass'] + field['mass'])
    r, m, fd = sites['r'], sites['mass'], sites['disk_fraction']
    a = rd*I.CP['scale_to_disk']
    T = I.CP['k0_per_kpc']*a
    vc = pot.vcirc(r)
    mom = [inc.moments(lambda mu, x=x: inc.reference_attenuation(x, T, mu)) for x in r/a]
    uA = np.array([q['u'] for q in mom])
    p_perp = np.array([q['p_perp'] for q in mom])
    power = m*kin.REFERENCE_POWER*uA
    J = float(np.sum(m*fd*r*vc))
    carried = float(np.sum(power*fd*r*vc)/power.sum())
    external = float(np.sum(power*fd*(1 + (1 + kin.REFERENCE_CHI)*p_perp)*r*vc)/power.sum())
    return dict(floor_ideal=M_ret*carried/J, external_ideal=M_ret*external/J)


def sparc_floor(d):
    c = d['catalog']
    out = floor(I.sparc_receivers(d['rotmod'], c), I.reference_field(c['rd'], c['L9']), d['total_inventory'], c['rd'])
    out['galaxy'] = d['name']
    return out


def mw_floor(run):
    out = floor(I.milky_way_receivers(run['baryons']), I.reference_field(run['rd'], run['L9']),
                run['total_inventory'], run['rd'])
    out.update(baryons=run['baryons'], Rd_kpc=run['rd'], luminosity_proxy_factor=run['lf'])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir', type=Path)
    ap.add_argument('--canonical', action='store_true')
    ap.add_argument('--workers', type=int, default=max(1, min(8, (os.cpu_count() or 2) - 2)))
    args = ap.parse_args()
    out_dir = (args.output_dir or I.ROOT/'research_work/generated'/(time.strftime('%Y%m%d-%H%M%S') + '-self-illumination-'
                                                                    + uuid.uuid4().hex[:6])).resolve()
    if out_dir.exists() and any(out_dir.iterdir()):
        ap.error('Output directory must be new or empty.')
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    result = dict(controls=controls(), ring=ring_runs(), disk=disk_runs())
    with ProcessPoolExecutor(args.workers) as pool:
        rows = list(pool.map(sparc_floor, I.sparc_galaxies()))
        mw = list(pool.map(mw_floor, [r for r in I.milky_way_runs() if r['rd'] == 2.6 and r['lf'] == 1]))
    f, e = np.array([r['floor_ideal'] for r in rows]), np.array([r['external_ideal'] for r in rows])
    result['transfer'] = dict(
        floor_ideal=dict(median=float(np.median(f)), p10=float(np.percentile(f, 10)), p90=float(np.percentile(f, 90)),
                         min=float(f.min()), max=float(f.max())),
        floor_above_one=int((f > 1).sum()), n_galaxies=len(rows), external_ideal_median=float(np.median(e)),
        floor_over_external_median=float(np.median(f/e)), milky_way=mw, rows=rows)
    result.update(
        scope='Self-illumination pilot: RB-1 reaction under companion fields derived from declared rotating emitters '
              '(steady optically thin transport); conservation floor transferred to the galaxy sample. Nothing fitted.',
        protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.glob('*.py'))},
        all_controls_passed=bool(result['controls']['passed'] and result['ring']['interior_zero_torque_check']),
        runtime_seconds=time.time() - t0)
    text = json.dumps(result, indent=1) + '\n'
    (out_dir/'pilot-results.json').write_text(text, encoding='utf-8', newline='\n')
    if args.canonical:
        (HERE/'pilot-results.json').write_text(text, encoding='utf-8', newline='\n')
    print(json.dumps(dict(controls_passed=result['all_controls_passed'],
                          transfer={k: v for k, v in result['transfer'].items() if k != 'rows'}), indent=1))
    print(f'Results: {out_dir}')
    return 0 if result['all_controls_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

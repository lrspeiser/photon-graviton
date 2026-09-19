"""Suite job for CL-2 stage 1 (protocol-cl2.md). Reruns gates G1, G2, G3, G5 and G7 in full, G4 and G6 at reduced
size (two clusters, three galaxies, one lens), G8 on the clusters-plus-galaxies solve, and a regression anchor
against cl2-results.json. Exits non-zero if any gate fails, any control is not rejected, or an anchored number
moves. About two minutes on one core."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl2 as C  # noqa: E402
import cl2_response as CR  # noqa: E402
import cl2_sources as CS  # noqa: E402

W = CR.WIDTHS


def main():
    extract = CS.load_xcop()
    clusters = CS.cluster_block(extract, W)
    train = CS.sparc_block(W, 'train')
    g = dict(G1_analytic_members=C.gate_g1(), G2_single_gaussian_limit=C.gate_g2(), G3_mixture_routes=C.gate_g3(),
             G5_galaxy_anchor=C.gate_g5(train), G7_xcop_ingestion=C.gate_g7(extract, clusters))
    # G4 at reduced size
    keep = W >= CR.CLUSTER_MIN_WIDTH
    frac = CS.stellar_fraction_profile(extract)
    worst = 0.
    for name in ('A1795', 'A2319'):
        a = CS.build_cluster(name, extract['clusters'][name], W, n_grid=1500, frac_profile=frac)['ops']
        b = CS.build_cluster(name, extract['clusters'][name], W, n_grid=3000, frac_profile=frac)['ops']
        worst = max(worst, float(np.max(np.abs(a[:, keep] - b[:, keep])/np.max(np.abs(b[:, keep]), axis=0))))
    worst_g = 0.
    keep_g = W >= CR.GALAXY_MIN_WIDTH
    gals = [x for x in CS.I.sparc_galaxies() if x['split'] == 'train'][:3]
    for gal in gals:
        R = CS._galaxy_arrays(gal)[0]
        a, b = CS.galaxy_operator(gal, W, R, 3000, exclude_narrow=False), CS.galaxy_operator(gal, W, R, 6000, exclude_narrow=False)
        worst_g = max(worst_g, float(np.max((np.max(np.abs(a - b), axis=0)/np.max(np.abs(b), axis=0))[keep_g])))
    g['G4_grid_convergence'] = dict(cluster_worst=worst, galaxy_worst=worst_g, reduced='two clusters, three galaxies', passed=bool(worst < 1e-4 and worst_g < 1e-3))
    # G6 on one lens
    L = CS.LensSystem(CS.LENSES[0], widths=W[:1])
    g['G6_lens_anchor'] = C.gate_g6([L])
    # G8 on the clusters + galaxies solve
    sol = CR.solve([clusters, train])
    g['G8_solve'] = dict(kkt=sol['kkt'], passed=bool(sol['kkt']['max_abs_gradient_active'] < 1e-10 and sol['kkt']['min_gradient_inactive'] >= -1e-10))
    for k, v in g.items():
        C.log(f"{k}: {'pass' if v['passed'] else 'FAIL'}" + (f", control {'rejected' if v['control']['rejected'] else 'NOT REJECTED'}" if 'control' in v else ''))
    verified = all(v['passed'] for v in g.values()) and all(v['control']['rejected'] for v in g.values() if 'control' in v)
    # the anchor
    sol_c, sol_g = CR.solve([clusters]), CR.solve([train])
    sg = CS.galaxy_scores(train, train.A@sol_g['amplitudes'])
    now = dict(G1_flat=g['G1_analytic_members']['log_flat_worst'], G1_inv=g['G1_analytic_members']['inverse_width_worst'],
               G2_dlam=g['G2_single_gaussian_limit']['relative_differences'][0], G3_worst=g['G3_mixture_routes']['worst'],
               G5_baryons=g['G5_galaxy_anchor']['measured']['baryons'], G7_floor=g['G7_xcop_ingestion']['release_nfw_floor_chi2_per_point'],
               E1_clusters_alone_chi2=sol_c['chi2'][0], E1_galaxies_alone_rmse=sg['equal_galaxy_rmse_km_s'], E1_galaxies_alone_slope=sg['mass_speed_slope_model'])
    arch = json.loads((HERE/'cl2-results.json').read_text(encoding='utf-8'))['checks_short_run']
    same = all(abs(now[k] - arch[k]) <= 1e-9*max(1., abs(arch[k])) for k in now)
    out = dict(gates={k: dict(passed=v['passed'], control_rejected=v.get('control', {}).get('rejected')) for k, v in g.items()},
               numerical_verification_passed=bool(verified), anchor=now, anchor_matches_archive=bool(same), passed=bool(verified and same))
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())

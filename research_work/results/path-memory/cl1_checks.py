"""Suite job for CL-1 (protocol-cl1.md). It reruns the nine numerical gates with their negative controls --
G9 at reduced size, one lens at one width for its lens half -- and a regression anchor against
cl1-results.json. Exits non-zero if any gate fails, any control is not rejected, or an anchored number
moves. About a minute on one core."""
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl1 as C  # noqa: E402


def main():
    inputs = json.loads(C.INPUTS.read_text(encoding='utf-8'))
    kubo = json.loads(C.KUBO.read_text(encoding='utf-8'))
    lenses = C.load_lenses()
    g = C.gates(inputs, kubo, lenses)
    # G9 at reduced size
    base = C.build_clusters(inputs)
    rows, best = C.cluster_scan(base, C.CLUSTER_WIDTHS)
    rows_d, _ = C.cluster_scan(C.build_clusters(inputs, n_grid=6000), C.CLUSTER_WIDTHS)
    conv_c = max(max(abs(a['chi2']/b['chi2'] - 1), abs(a['Lambda']/b['Lambda'] - 1)) for a, b in zip(rows, rows_d))
    L = lenses[0]
    Ld = C.Lens(L.name, C.LZ.load(), n_proj=16000)
    lam = L.lambda_for_observed(L.masses['Chabrier'], 10.)
    conv_l = abs(lam/Ld.lambda_for_observed(Ld.masses['Chabrier'], 10.) - 1)
    g['G9_convergence'] = dict(cluster_scan_grid_doubling=conv_c, lens_lambda_grid_doubling=conv_l, reduced='one lens, w = 10 kpc',
                               tolerances=dict(clusters=1e-5, lenses=1e-4), passed=bool(conv_c < 1e-5 and conv_l < 1e-4))
    C.log(f"G9 (reduced): {'pass' if g['G9_convergence']['passed'] else 'FAIL'}")
    verified = all(v['passed'] for v in g.values()) and all(v['control']['rejected'] for v in g.values() if 'control' in v)
    # the regression anchor
    E3 = C.e3(kubo, best)
    now = dict(G1_worst=g['G1_shell_kernel']['worst_relative_difference'],
               G2_worst=max(g['G2_derivatives']['worst_gradient'], g['G2_derivatives']['worst_laplacian']),
               G3_w_over_ell=g['G3_link_to_rut1']['mature_ring_control']['w_over_ell'], G4_worst=g['G4_no_net_equivalent_mass']['worst'],
               G5_worst=g['G5_two_lensing_routes']['worst_relative_difference'], G6_worst_2d=g['G6_lens_anchor']['worst_2d'],
               G7_worst_mass=g['G7_xcop_ingestion']['worst_mass_deviation'],
               E1_best_w=best['w'], E1_best_Lambda=best['Lambda'], E1_best_chi2=best['chi2'],
               E1_newtonian_chi2=C.cluster_references(base)['newtonian_baryons']['chi2'],
               E2_J0037_Chabrier_Lambda_at_10kpc=lam, E2_width_index_10kpc=12,
               E3_low_bracket_L0_chi2=E3['brackets'][0]['L0']['shape_chi2'])
    arch = json.loads((HERE/'cl1-results.json').read_text(encoding='utf-8'))['checks_short_run']
    same = all(abs(now[k] - arch[k]) <= 1e-9*max(1., abs(arch[k])) for k in now)
    out = dict(gates={k: dict(passed=v['passed'], control_rejected=v.get('control', {}).get('rejected')) for k, v in g.items()},
               numerical_verification_passed=bool(verified), anchor=now, anchor_matches_archive=bool(same),
               passed=bool(verified and same))
    print(json.dumps(out, indent=1, default=float))
    return 0 if out['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())

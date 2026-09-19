"""CL-2 stage 1: one written-field response tested jointly on galaxies, lenses and clusters (protocol-cl2.md).

    python cl2.py [--output-dir DIR] [--canonical]

Runs the eight numerical gates with their controls, then E1 (the universal spectrum per block and jointly, with
the declared decision rule), E2 (the lens tests), E3 (transfer) and E4 (the cluster diagnostics); writes
cl2-results.json and cl2-geometry.json and compares the results with the archive (the first run creates it).
About twelve minutes on one core; cl2_checks.py is the suite job.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky, solve_triangular

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(RESULTS/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE))
import steady_field as SF  # noqa: E402
import cl1 as C1  # noqa: E402
import cl2_response as CR  # noqa: E402
import cl2_sources as CS  # noqa: E402
import cl2_xcop_acquire as ACQ  # noqa: E402

G = SF.G
W = CR.WIDTHS
T0 = time.time()
CR2_ARCHIVE = RESULTS/'supported-reservoir/cr2-results.json'
DESCRIBED = dict(cluster_floor_factor=2., galaxy_rmse_factor=1.1, galaxy_slope_tolerance=.05, lens_einstein_tolerance=.03, lens_nfw_factor=1.5)


def log(msg):
    print(f'[{time.time() - T0:6.1f}s] {msg}', flush=True)


def scale_rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b))/np.max(np.abs(b)))


# ================================================================ gates
def gate_g1():
    src = SF.SphericalSource(np.array([1e-9, 2e-9]), np.array([1., 1.]))
    ws = np.geomspace(1e-3, 1e3, 601)
    dlnw = float(np.log(ws[1]/ws[0]))
    r = np.array([.03, .3, 3.])
    flat = CR.point_source_sum(src, r, ws, np.full(len(ws), dlnw))
    e_flat = float(np.max(np.abs(flat/CR.log_flat_force(1., r, ws[0], ws[-1]) - 1)))
    inv = CR.point_source_sum(src, r, ws, dlnw/ws)
    e_inv = float(np.max(np.abs(inv/CR.inverse_width_force(1., r) - 1)))
    ctrl = float(np.max(np.abs(inv/CR.inverse_width_force(1., r, np.sqrt(np.pi)) - 1)))
    return dict(log_flat_worst=e_flat, inverse_width_worst=e_inv, tolerance=1e-6, passed=bool(e_flat < 1e-6 and e_inv < 1e-6),
                control=dict(name='the 1/w member against sqrt(pi) M/r^2 in place of sqrt(pi/2)', worst=ctrl, rejected=bool(ctrl > 1e-6)))


def gate_g2():
    inputs = json.loads(C1.INPUTS.read_text(encoding='utf-8'))
    arch = json.loads((HERE/'cl1-results.json').read_text(encoding='utf-8'))['E1_clusters']['best']
    cls = C1.build_clusters(inputs)
    w = arch['w']
    A = np.vstack([np.array([c['source'].g_mem(c['radii'], w)]).T for c in cls])
    y = np.concatenate([c['g_req'] - c['g_N'] for c in cls])
    s = np.concatenate([c['sigma'] for c in cls])
    one = CR.solve([CR.Block('cl1', A, y, s)], weights=[1.])
    A2 = np.hstack([A, np.vstack([np.array([c['source'].g_mem(c['radii'], 100.)]).T for c in cls])])
    two = CR.solve([CR.Block('cl1', A2, y, s)], weights=[1.])
    d_lam, d_chi = abs(one['amplitudes'][0]/arch['Lambda'] - 1), abs(one['chi2'][0]/arch['chi2'] - 1)
    return dict(cl1_archived=dict(w=w, Lambda=arch['Lambda'], chi2=arch['chi2']), spectrum_one_width=dict(Lambda=float(one['amplitudes'][0]), chi2=one['chi2'][0]),
                relative_differences=[d_lam, d_chi], tolerance=1e-9, passed=bool(d_lam < 1e-9 and d_chi < 1e-9),
                control=dict(name='a second width (100 kpc) with free amplitude', chi2=two['chi2'][0], rejected=bool(abs(two['chi2'][0]/arch['chi2'] - 1) > 1e-9)))


def gate_g3():
    widths, amps = (1., 3., 10.), (1., .3, .1)
    r3 = np.geomspace(1e-4, 300., 8000)
    src = SF.SphericalSource(r3, SF.plummer_mass(1e12, 5., r3))
    R2 = np.geomspace(1e-3, 300., 6000)
    proj = SF.ProjectedSource(R2, SF.plummer_projected_mass(1e12, 5., R2))
    worst, ctrl = 0., 1.
    rows = []
    for b in (1., 4., 10., 25.):
        a3 = SF.deflection_from_g(lambda x: sum(a*src.g_mem([x], w)[0] for w, a in zip(widths, amps)), b)
        a2 = SF.deflection_from_projected_mass(sum(a*proj.projected_equivalent_mass([b], w)[0] for w, a in zip(widths, amps)), b)
        ac = SF.deflection_from_projected_mass(sum(a*proj.projected_equivalent_mass([b], w, kernel_dR=SF.shell_kernel_dr)[0] for w, a in zip(widths, amps)), b)
        rows.append(dict(b=b, three_dimensional=a3, two_dimensional=a2, relative_difference=abs(a3/a2 - 1), control=abs(ac/a3 - 1)))
        worst, ctrl = max(worst, abs(a3/a2 - 1)), min(ctrl, abs(ac/a3 - 1))
    return dict(rows=rows, worst=worst, tolerance=1e-6, passed=bool(worst < 1e-6),
                control=dict(name='the two-dimensional route with the shell kernel', smallest=ctrl, rejected=bool(ctrl > 1e-6)))


def gate_g4(extract, gal_block):
    keep = W >= CR.CLUSTER_MIN_WIDTH
    frac = CS.stellar_fraction_profile(extract)
    worst, excluded = 0., 0.
    rows = []
    for name, c in extract['clusters'].items():
        a = CS.build_cluster(name, c, W, n_grid=1500, frac_profile=frac)['ops']
        b = CS.build_cluster(name, c, W, n_grid=3000, frac_profile=frac)['ops']
        d = float(np.max(np.abs(a[:, keep] - b[:, keep])/np.max(np.abs(b[:, keep]), axis=0)))
        # the excluded widths, computed here only to report their non-convergence
        rows.append(dict(cluster=name, retained_widths=d))
        worst = max(worst, d)
    keep_g = W >= CR.GALAXY_MIN_WIDTH
    worst_g, worst_excluded = 0., 0.
    for p, gal in zip(gal_block.meta['galaxies'], [g for g in CS.I.sparc_galaxies() if g['split'] == gal_block.meta['split']]):
        a = CS.galaxy_operator(gal, W, p['R'], n_grid=3000, exclude_narrow=False)
        b = CS.galaxy_operator(gal, W, p['R'], n_grid=6000, exclude_narrow=False)
        d = np.max(np.abs(a - b), axis=0)/np.max(np.abs(b), axis=0)
        worst_g, worst_excluded = max(worst_g, float(np.max(d[keep_g]))), max(worst_excluded, float(np.max(d[~keep_g])))
    return dict(cluster_rows=rows, cluster_worst=worst, galaxy_worst_retained=worst_g, galaxy_worst_excluded=worst_excluded,
                tolerances=dict(clusters=1e-4, galaxies=1e-3), excluded_cluster_widths_kpc=W[~keep].tolist(), excluded_galaxy_widths_kpc=W[~keep_g].tolist(),
                declared_gate_G4='FAILED as declared on the excluded galaxy width (amendment 2); G4b applies to the retained columns',
                passed=bool(worst < 1e-4 and worst_g < 1e-3))


def gate_g5(train):
    refs = CS.galaxy_references(train)
    got = dict(baryons=refs['baryons']['equal_galaxy_rmse_km_s'], pm1_fitted=refs['pm1_law_fitted_a_star']['equal_galaxy_rmse_km_s'],
               simple_mond=refs['simple_mond_repo_a_star']['equal_galaxy_rmse_km_s'])
    arch = dict(baryons=52.57, pm1_fitted=20.21, simple_mond=19.89)
    ctrl = CS.galaxy_references(CS.sparc_block(W[:2], 'train', n_grid=300, upsilon_disk=1.))['baryons']['equal_galaxy_rmse_km_s']
    return dict(measured=got, archived=arch, tolerance=.02, passed=bool(all(abs(got[k] - arch[k]) <= .02 for k in got)),
                control=dict(name='Upsilon_disk = 1 must move the baryon score by more than 5 km/s', baryons_rmse=ctrl, rejected=bool(abs(ctrl - got['baryons']) > 5)))


def gate_g6(lenses):
    arch = json.loads(CR2_ARCHIVE.read_text(encoding='utf-8'))['geometries']['G1_flat_FLRW']['benchmarks']
    rows, worst, ctrl = [], 0., 1.
    for L in lenses:
        so = L.stars_only()
        a = arch[L.name]['stars_only']
        d = max(abs(so['chi2']/a['chi2'] - 1), abs(so['beta'] - a['beta']))
        rows.append(dict(lens=L.name, chi2=so['chi2'], archived_chi2=a['chi2'], beta=so['beta'], archived_beta=a['beta'], difference=d))
        worst = max(worst, d)
    Lp = CS.LensSystem(lenses[0].name, psf_scale=2., widths=W[:1])
    ctrl = abs(Lp.stars_only()['chi2']/rows[0]['chi2'] - 1)
    return dict(rows=rows, worst=worst, tolerance=1e-6, passed=bool(worst < 1e-6),
                control=dict(name='the PSF doubled must move chi2 by more than 1%', relative_change=ctrl, rejected=bool(ctrl > .01)))


def gate_g7(extract, block):
    prov_ok = extract['provenance']['sha256'] == ACQ.SHA256
    checks = {c['name']: c['gas_mass_check'] for c in block.meta['clusters']}
    refs = CS.cluster_references(block)
    floor = refs['release_nfw_floor']['chi2_per_point']
    frac = CS.stellar_fraction_profile(extract)
    ctrl_cls = [CS.build_cluster(n, c, W, frac_profile=frac, mu=1.2) for n, c in extract['clusters'].items()]
    ctrl = sum(CS.fit_boundary(c, c['P_nfw'])[0] for c in ctrl_cls)/sum(len(c['rp']) for c in ctrl_cls)
    return dict(archive_hash_pinned=prov_ok, gas_mass_integral_over_release=checks, release_nfw_floor_chi2_per_point=floor,
                tolerances=dict(gas_mass=.05, floor=10.), passed=bool(prov_ok and all(abs(v - 1) < .05 for v in checks.values()) and floor < 10),
                control=dict(name='mu = 1.2 must raise the floor above 10 times its value', chi2_per_point=ctrl, rejected=bool(ctrl > 10*floor)))


# ================================================================ lenses: rows, fits and predictions
def lens_rows(lenses, imf, betas):
    """Einstein block (one row per lens) and whitened kinematics block at the given betas."""
    Ae, ye, se, Ak, yk = [], [], [], [], []
    for L, beta in zip(lenses, betas):
        m = L.pop[imf]/1e11
        a, y, s = L.einstein_row(m)
        Ae.append(a); ye.append(y); se.append(s)
        A, yy = L.kinematic_rows(m, beta)
        Lc = cholesky(np.diag(2*L.y)@L.cov@np.diag(2*L.y), lower=True)
        Ak.append(solve_triangular(Lc, A, lower=True)); yk.append(solve_triangular(Lc, yy, lower=True))
    E = CR.Block('lens_einstein', np.array(Ae), np.array(ye), np.array(se))
    K = CR.Block('lens_kinematics', np.vstack(Ak), np.concatenate(yk), np.ones(sum(len(v) for v in yk)))
    return E, K


def lens_evaluate(lenses, imf, amplitudes):
    """At fixed amplitudes: Einstein residuals, and the kinematics chi2 (in V, release covariance) with beta refitted."""
    out = []
    for L in lenses:
        m = L.pop[imf]/1e11
        a, y, s = L.einstein_row(m)
        chi, beta = L.fit_beta(m, amplitudes)
        out.append(dict(lens=L.name, imf=imf, einstein_relative_residual=float((a@amplitudes - y)/L.need), kinematics_chi2=chi, beta=beta, n_bins=len(L.y)))
    return out


def alternating_solve(blocks_fixed, lenses, imf, betas0, use_einstein=True, use_kinematics=True, iters=6):
    """Nonnegative solve over fixed blocks plus the lens blocks, alternating with per-lens beta refits."""
    betas = list(betas0)
    best = None
    for it in range(iters):
        E, K = lens_rows(lenses, imf, betas)
        blocks = list(blocks_fixed) + ([E] if use_einstein else []) + ([K] if use_kinematics else [])
        sol = CR.solve(blocks)
        ev = lens_evaluate(lenses, imf, sol['amplitudes'])
        betas_new = [e['beta'] for e in ev]
        objective = sum(c/b.N for c, b in zip(sol['chi2'], blocks))
        if best is not None and abs(objective - best[0]) < 1e-6*max(1., abs(best[0])):
            best = (objective, sol, ev, blocks); break
        best = (objective, sol, ev, blocks)
        betas = betas_new
    return best


def described(block_name, sol_chi2, N, refs=None, scores=None, lens_ev=None, nfw_total=None):
    if block_name == 'clusters':
        return bool(sol_chi2/N <= DESCRIBED['cluster_floor_factor']*refs['release_nfw_floor']['chi2_per_point'])
    if block_name.startswith('sparc'):
        return bool(scores['equal_galaxy_rmse_km_s'] <= DESCRIBED['galaxy_rmse_factor']*refs['simple_mond_repo_a_star']['equal_galaxy_rmse_km_s']
                    and abs(scores['mass_speed_slope_model'] - scores['mass_speed_slope_observed']) <= DESCRIBED['galaxy_slope_tolerance'])
    if block_name == 'lenses':
        return bool(all(abs(e['einstein_relative_residual']) <= DESCRIBED['lens_einstein_tolerance'] for e in lens_ev)
                    and sum(e['kinematics_chi2'] for e in lens_ev) <= DESCRIBED['lens_nfw_factor']*nfw_total)
    return None


# ================================================================ driver
def main():
    args = evidence_io.parse(__doc__)
    extract = CS.load_xcop()
    log('building blocks')
    clusters = CS.cluster_block(extract, W)
    train = CS.sparc_block(W, 'train')
    valid, test = CS.sparc_block(W, 'validation'), CS.sparc_block(W, 'test')
    mw = CS.milky_way_blocks(W)
    lenses = [CS.LensSystem(n, widths=W) for n in CS.LENSES]
    try:
        held = CS.LensSystem(CS.HELD_OUT_LENS, widths=W)
    except KeyError:
        held = None                       # amendment 1: no published light-profile components; withdrawn
    log('blocks built')
    cref, gref = CS.cluster_references(clusters), CS.galaxy_references(train)
    cr2 = json.loads(CR2_ARCHIVE.read_text(encoding='utf-8'))['geometries']['G1_flat_FLRW']['benchmarks']
    nfw_total = sum(cr2[n]['nfw']['chi2'] for n in CS.LENSES)
    betas0 = [cr2[n]['stars_only']['beta'] for n in CS.LENSES]
    gates = dict(G1_analytic_members=gate_g1(), G2_single_gaussian_limit=gate_g2(), G3_mixture_routes=gate_g3(),
                 G4_grid_convergence=gate_g4(extract, train), G5_galaxy_anchor=gate_g5(train), G6_lens_anchor=gate_g6(lenses),
                 G7_xcop_ingestion=gate_g7(extract, clusters))
    for k, v in gates.items():
        log(f"{k}: {'pass' if v['passed'] else 'FAIL'}" + (f", control {'rejected' if v['control']['rejected'] else 'NOT REJECTED'}" if 'control' in v else ''))

    # ---------------- E1: the universal spectrum, per block and jointly
    def gal_eval(block, amps):
        return CS.galaxy_scores(block, block.A@amps)
    E1 = dict(widths_kpc=W.tolist(), excluded_cluster_widths_kpc=W[W < CR.CLUSTER_MIN_WIDTH].tolist(), excluded_galaxy_widths_kpc=W[W < CR.GALAXY_MIN_WIDTH].tolist(), references=dict(clusters=cref, galaxies_train=gref),
              decision_rule=DESCRIBED)
    sol_c = CR.solve([clusters])
    E1['clusters_alone'] = dict(chi2=sol_c['chi2'][0], chi2_per_point=sol_c['chi2'][0]/clusters.N, spectrum=CR.spectrum_summary(W, sol_c['amplitudes'], G),
                                described=described('clusters', sol_c['chi2'][0], clusters.N, refs=cref),
                                per_cluster_chi2_per_point={c['name']: float(np.sum(((c['ops']@sol_c['amplitudes'] + nu*c['thermal'] + c['P_N']*c['thermal'] - c['Pobs'])/c['eP'])**2))/len(c['rp'])
                                                            for c, nu in zip(clusters.meta['clusters'], sol_c['nuisances'][0])})
    sol_g = CR.solve([train])
    sg = gal_eval(train, sol_g['amplitudes'])
    E1['galaxies_alone'] = dict(train=sg, validation=gal_eval(valid, sol_g['amplitudes']), test=gal_eval(test, sol_g['amplitudes']),
                                spectrum=CR.spectrum_summary(W, sol_g['amplitudes'], G), described=described('sparc', sg['chi2'], train.N, refs=gref, scores=sg))
    lens_alone = alternating_solve([], lenses, 'Chabrier', betas0)
    E1['lenses_alone'] = dict(objective=lens_alone[0], evaluation=lens_alone[2], spectrum=CR.spectrum_summary(W, lens_alone[1]['amplitudes'], G),
                              described=described('lenses', None, None, lens_ev=lens_alone[2], nfw_total=nfw_total), cr2_free_nfw_total_chi2=nfw_total)
    joint = alternating_solve([clusters, train], lenses, 'Chabrier', betas0)
    obj, sol_j, ev_j, blocks_j = joint
    sj = gal_eval(train, sol_j['amplitudes'])
    E1['joint'] = dict(objective=obj, chi2_per_block={b.name: c for b, c in zip(blocks_j, sol_j['chi2'])}, N_per_block={b.name: b.N for b in blocks_j},
                       clusters=dict(chi2_per_point=sol_j['chi2'][0]/clusters.N, described=described('clusters', sol_j['chi2'][0], clusters.N, refs=cref)),
                       galaxies=dict(train=sj, validation=gal_eval(valid, sol_j['amplitudes']), test=gal_eval(test, sol_j['amplitudes']),
                                     described=described('sparc', sj['chi2'], train.N, refs=gref, scores=sj)),
                       lenses=dict(evaluation=ev_j, described=described('lenses', None, None, lens_ev=ev_j, nfw_total=nfw_total)),
                       spectrum=CR.spectrum_summary(W, sol_j['amplitudes'], G), kkt=sol_j['kkt'])
    alone_ok = E1['clusters_alone']['described'] and E1['galaxies_alone']['described'] and E1['lenses_alone']['described']
    joint_ok = E1['joint']['clusters']['described'] and E1['joint']['galaxies']['described'] and E1['joint']['lenses']['described']
    E1['outcome'] = ('a: a shared spectrum works provisionally' if joint_ok else 'b: only separate spectra work' if alone_ok else 'c: even flexible positive spectra fail')
    E1['outcome_detail'] = dict(clusters_alone=E1['clusters_alone']['described'], galaxies_alone=E1['galaxies_alone']['described'],
                                lenses_alone=E1['lenses_alone']['described'], joint=[E1['joint'][k]['described'] for k in ('clusters', 'galaxies', 'lenses')])
    log(f"E1 outcome: {E1['outcome']}")
    # width-grid and cluster sensitivities on the clusters+galaxies joint solve
    sens = {}
    for label, grid in CR.WIDTH_GRIDS.items():
        if label == 'primary':
            continue
        cb, gb = CS.cluster_block(extract, grid), CS.sparc_block(grid, 'train')
        s = CR.solve([cb, gb])
        sens[label] = dict(clusters_chi2_per_point=s['chi2'][0]/cb.N, galaxies=CS.galaxy_scores(gb, gb.A@s['amplitudes']), active_widths=[float(grid[i]) for i in s['widths_active']])
    for label, kw in (('non_thermal_pressure', dict(alpha_nt=True)), ('mu_0.59', dict(mu=.59)), ('mu_0.61', dict(mu=.61)), ('stars_none', dict(stars='none'))):
        cb = CS.cluster_block(extract, W, **kw)
        s = CR.solve([cb])
        sens[label] = dict(clusters_alone_chi2_per_point=s['chi2'][0]/cb.N, floor=CS.cluster_references(cb)['release_nfw_floor']['chi2_per_point'],
                           spectrum=CR.spectrum_summary(W, s['amplitudes'], G))
    E1['sensitivities'] = sens
    log('E1 sensitivities done')

    # ---------------- E2: the lens tests
    E2 = dict(cr2_benchmarks={n: cr2[n] for n in CS.LENSES}, held_out=CS.HELD_OUT_LENS if held is not None else 'withdrawn (amendment 1): no published light-profile components')
    for imf in ('Chabrier', 'Salpeter'):
        rows = {}
        lo = alternating_solve([], lenses, imf, betas0, use_einstein=True, use_kinematics=False)
        ko = alternating_solve([], lenses, imf, betas0, use_einstein=False, use_kinematics=True)
        jo = alternating_solve([], lenses, imf, betas0)
        for label, res in (('lensing_only', lo), ('kinematics_only', ko), ('joint', jo)):
            amps = res[1]['amplitudes']
            rows[label] = dict(spectrum=CR.spectrum_summary(W, amps, G), evaluation=res[2], held_out=lens_evaluate([held], imf, amps)[0] if held is not None else 'not available',
                               Lambda_at_4p64_kpc=float(amps[10]), kinematics_chi2_total=float(sum(e['kinematics_chi2'] for e in res[2])),
                               einstein_max_abs_residual=float(max(abs(e['einstein_relative_residual']) for e in res[2])))
        E2[imf] = rows
    E2['held_out_stars_only'] = held.stars_only() if held is not None else 'not available'
    log('E2 done')

    # ---------------- E3: transfer without retuning
    def cluster_fixed(amps):
        per = {c['name']: CS.fit_boundary(c, c['ops']@amps + c['P_N'])[0] for c in clusters.meta['clusters']}
        return dict(chi2_per_point=float(sum(per.values()))/clusters.N, per_cluster={k: v/len(c['rp']) for (k, v), c in zip(per.items(), clusters.meta['clusters'])})
    def mw_eval(amps):
        return {b.name: dict(rmse_km_s=float(np.sqrt(np.mean((np.sqrt(np.maximum(b.meta['R']*(b.meta['gN'] + b.A@amps), 0.)) - b.meta['y'])**2))),
                             baryons_rmse=float(np.sqrt(np.mean((np.sqrt(b.meta['R']*b.meta['gN']) - b.meta['y'])**2)))) for b in mw if b is not None}
    E3 = {}
    E3['calibrated_on_galaxies'] = dict(clusters=cluster_fixed(sol_g['amplitudes']), lenses=lens_evaluate(lenses, 'Chabrier', sol_g['amplitudes']),
                                        held_out=lens_evaluate([held], 'Chabrier', sol_g['amplitudes'])[0] if held is not None else 'not available', milky_way=mw_eval(sol_g['amplitudes']))
    E3['calibrated_on_clusters'] = dict(galaxies=dict(train=gal_eval(train, sol_c['amplitudes']), validation=gal_eval(valid, sol_c['amplitudes']), test=gal_eval(test, sol_c['amplitudes'])),
                                        lenses=lens_evaluate(lenses, 'Chabrier', sol_c['amplitudes']), held_out=lens_evaluate([held], 'Chabrier', sol_c['amplitudes'])[0] if held is not None else 'not available',
                                        milky_way=mw_eval(sol_c['amplitudes']))
    cg = CR.solve([clusters, train])
    E3['calibrated_on_galaxies_and_clusters'] = dict(spectrum=CR.spectrum_summary(W, cg['amplitudes'], G), clusters_chi2_per_point=cg['chi2'][0]/clusters.N,
                                                     galaxies=dict(train=gal_eval(train, cg['amplitudes']), validation=gal_eval(valid, cg['amplitudes']), test=gal_eval(test, cg['amplitudes'])),
                                                     lenses=lens_evaluate(lenses, 'Chabrier', cg['amplitudes']), held_out=lens_evaluate([held], 'Chabrier', cg['amplitudes'])[0] if held is not None else 'not available',
                                                     milky_way=mw_eval(cg['amplitudes']))
    log('E3 done')

    # ---------------- E4: cluster diagnostics
    E4 = dict(per_cluster=[], truncation=dict(), cl1_regression_case=None)
    for c in clusters.meta['clusters']:
        E4['per_cluster'].append(dict(cluster=c['name'], stars=c['star_source'], points=len(c['rp']), r_out_kpc=c['r_out'],
                                      chi2_per_point=dict(release_nfw=cref['release_nfw_floor']['per_cluster_chi2_per_point'][c['name']],
                                                          baryons=cref['newtonian_baryons']['per_cluster_chi2_per_point'][c['name']],
                                                          pm1_law=cref['pm1_law_repo_a_star']['per_cluster_chi2_per_point'][c['name']],
                                                          spectrum_alone=E1['clusters_alone']['per_cluster_chi2_per_point'][c['name']]),
                                      shell_decomposition=CS.shell_decomposition(c, W, sol_c['amplitudes'])))
    for x in (1., 1.5, 2., 3.):
        cb = CS.cluster_block(extract, W, truncate_gas_at=x)
        E4['truncation'][f'gas_truncated_at_{x:g}_R500'] = dict(chi2_per_point_at_fixed_amplitudes=float(sum(CS.fit_boundary(c, c['ops']@sol_c['amplitudes'] + c['P_N'])[0] for c in cb.meta['clusters']))/cb.N,
                                                                  refit=CR.solve([cb])['chi2'][0]/cb.N)
    inputs = json.loads(C1.INPUTS.read_text(encoding='utf-8'))
    cl1 = {c['name']: c for c in C1.build_clusters(inputs)}
    chi_new, chi_n, chi_p = 0., 0., 0.
    for c in clusters.meta['clusters']:
        k = cl1[c['name']]
        gN = G*np.interp(k['radii'], c['r'], c['M_b'])/k['radii']**2
        gm = c['source'].g_mem(k['radii'], 1.)*0.
        for w, a in zip(W, sol_c['amplitudes']):
            if a > 0 and w >= CR.CLUSTER_MIN_WIDTH:
                gm = gm + a*c['source'].g_mem(k['radii'], w)
        chi_new += float(np.sum(((gN + gm - k['g_req'])/k['sigma'])**2))
        chi_n += float(np.sum(((gN - k['g_req'])/k['sigma'])**2))
        chi_p += float(np.sum(((gN + np.sqrt(CS.A_STAR_REPO*gN) - k['g_req'])/k['sigma'])**2))
    E4['cl1_regression_case'] = dict(note='CL-1 five-point NFW comparison with the release baryons and the clusters-alone spectrum, at fixed amplitudes',
                                     spectrum_chi2=chi_new, newtonian_chi2=chi_n, pm1_law_chi2=chi_p, cl1_archived_best_single_width_chi2=json.loads((HERE/'cl1-results.json').read_text())['E1_clusters']['best']['chi2'])
    log('E4 done')

    geometry = {L.name: {sc: CS.geometry(L.OBS[L.name]['zFG'], L.OBS[L.name]['zBG'], sc) for sc in ('G1_flat_FLRW', 'PF1_static_euclidean', 'G2_coscaling_coasting')}
                for L in lenses + ([held] if held is not None else [])}
    for name, rec in geometry.items():
        rec['conditional_geometry_PF1_archived'] = dict(Dl_kpc=CS.LensSystem.COND[name]['conditional_Dl_Mpc']*1000, Dls_over_Ds=CS.LensSystem.COND[name]['conditional_Dls_over_Ds'])
        rec['PF1_matches_archived'] = bool(abs(rec['PF1_static_euclidean']['Dl_kpc']/rec['conditional_geometry_PF1_archived']['Dl_kpc'] - 1) < 1e-6
                                           and abs(rec['PF1_static_euclidean']['Dls_over_Ds'] - rec['conditional_geometry_PF1_archived']['Dls_over_Ds']) < 1e-6)
    gates['G8_solve'] = dict(kkt=sol_j['kkt'], tolerance=1e-10, passed=bool(sol_j['kkt']['max_abs_gradient_active'] < 1e-10 and sol_j['kkt']['min_gradient_inactive'] >= -1e-10))
    verified = all(v['passed'] for v in gates.values()) and all(v['control']['rejected'] for v in gates.values() if 'control' in v)
    result = dict(experiment='CL-2 stage 1: one written-field response tested jointly on galaxies, lenses and clusters', protocol='protocol-cl2.md',
                  statuses=dict(numerical_verification='passed' if verified else 'FAILED', scientific_outcome=E1['outcome']),
                  gates=gates, numerical_verification_passed=bool(verified), E1_spectrum=E1, E2_lenses=E2, E3_transfer=E3, E4_cluster_diagnostics=E4,
                  primary_geometry='G1_flat_FLRW', input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                                                                   (HERE/'protocol-cl2.md', CS.XCOP, HERE/'cl2_response.py', HERE/'cl2_sources.py')},
                  checks_short_run=dict(G1_flat=gates['G1_analytic_members']['log_flat_worst'], G1_inv=gates['G1_analytic_members']['inverse_width_worst'],
                                        G2_dlam=gates['G2_single_gaussian_limit']['relative_differences'][0], G3_worst=gates['G3_mixture_routes']['worst'],
                                        G5_baryons=gates['G5_galaxy_anchor']['measured']['baryons'], G7_floor=gates['G7_xcop_ingestion']['release_nfw_floor_chi2_per_point'],
                                        E1_clusters_alone_chi2=E1['clusters_alone']['chi2'], E1_galaxies_alone_rmse=E1['galaxies_alone']['train']['equal_galaxy_rmse_km_s'],
                                        E1_galaxies_alone_slope=E1['galaxies_alone']['train']['mass_speed_slope_model'], E1_joint_objective=E1['joint']['objective']),
                  runtime_seconds=round(time.time() - T0, 1))
    (HERE/'cl2-geometry.json').write_text(json.dumps(geometry, indent=1) + '\n', encoding='utf-8', newline='\n')
    text = json.dumps(result, indent=1, default=float) + '\n'
    canonical = HERE/'cl2-results.json'
    if not canonical.exists() or args.canonical:
        canonical.write_text(text, encoding='utf-8', newline='\n')
        log('first run: archive created')
    status = evidence_io.finish(args, 'path-memory-cl2', text, canonical, ignore={'/runtime_seconds'})
    print(json.dumps(dict(numerical_verification=result['statuses']['numerical_verification'], outcome=E1['outcome'], detail=E1['outcome_detail'],
                          clusters_alone=E1['clusters_alone']['chi2_per_point'], floor=cref['release_nfw_floor']['chi2_per_point'],
                          galaxies_alone=E1['galaxies_alone']['train'], joint=E1['joint']['chi2_per_block']), indent=1, default=float))
    return status if verified else 1


if __name__ == '__main__':
    raise SystemExit(main())

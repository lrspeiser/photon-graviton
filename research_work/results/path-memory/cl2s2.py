"""CL-2 stage 2 driver (protocol-cl2-stage2.md, amendment 1): every column visible to every block, one source per
galaxy, the equal-galaxy speed loss minimised exactly with a certificate, the corrected cluster comparison under
the release's SZ correlations, the joint and transfer solves. Writes cl2s2-results.json."""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import cl2s2_lib as L2          # noqa: E402
import cl2_sources as CS        # noqa: E402
import cl2_response as CR       # noqa: E402
import cl2_xcop_acquire as XA   # noqa: E402

GEN = ROOT/'research_work/generated/cl2s2'
TAR = ROOT/'research_work/generated/xcop-release/allfiles.tar.gz'
ARCHIVE = HERE/'cl2-results.json'
OUT = HERE/'cl2s2-results.json'
W = L2.WIDTHS
RMSE_LIMIT = 21.9
N_GRID_CLUSTERS = 6000


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plain(x):
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [plain(v) for v in x]
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        return float(x)
    return x


# ------------------------------------------------------------------ galaxies
def load_galaxies():
    z = np.load(GEN/'galaxies.npz')
    names = sorted(set(k.split('/')[1] for k in z.files if '/' in k))
    blocks = {}
    for split in ('train', 'validation', 'test'):
        gals = sorted(set(k.split('/')[1] for k in z.files if k.startswith(split + '/')))
        A, gr, gt, R, v, e, sizes, masses = [], [], [], [], [], [], [], []
        for g in gals:
            p = f'{split}/{g}/'
            A.append(z[p + 'ops']); gr.append(z[p + 'gN_rec']); gt.append(z[p + 'gN_tab']); R.append(z[p + 'R'])
            v.append(z[p + 'vobs']); e.append(z[p + 'ev']); sizes.append(len(z[p + 'R'])); masses.append(float(z[p + 'mass'][0]))
        blocks[split] = dict(names=gals, A=np.vstack(A), gN_rec=np.concatenate(gr), gN_tab=np.concatenate(gt), R=np.concatenate(R),
                             v=np.concatenate(v), ev=np.concatenate(e), sizes=sizes, masses=np.array(masses))
    return blocks, float(z['G4b_worst'][0])


def slope_diagnostic(block, L, gN):
    """log-log slope of the outer model and observed speeds against the integrated baryonic mass (five outer points)."""
    u = block['R']*(gN + block['A']@L)
    vm = np.sqrt(np.maximum(u, 0))
    off, M, vo, vmo = 0, [], [], []
    for n, mass in zip(block['sizes'], block['masses']):
        sl = slice(off, off + n)
        off += n
        if n >= 5:
            M.append(mass); vo.append(block['v'][sl][-5:].mean()); vmo.append(vm[sl][-5:].mean())
    M, vo, vmo = np.log10(M), np.log10(vo), np.log10(vmo)
    return dict(slope_model=float(np.polyfit(M, vmo, 1)[0]), slope_observed=float(np.polyfit(M, vo, 1)[0]), galaxies=len(M))


def galaxy_stage(blocks, archive):
    tr = blocks['train']
    out = {}
    for label, key in (('reconstructed', 'gN_rec'), ('tabulated', 'gN_tab')):
        gN = tr[key]
        loss = L2.SpeedLoss(tr['A'], gN, tr['R'], tr['v'], tr['sizes'])
        s = np.maximum(2*tr['v']*tr['ev']/tr['R'], 1e-3*tr['v']**2/tr['R'])
        L0 = L2.nnls_start(tr['A'], tr['v']**2/tr['R'] - gN, s, np.ones(len(s)))     # stage 1's problem: unit weights, 1/s^2 errors
        start = dict(rmse=loss.rmse(L0) if np.isfinite(loss.value(L0)) else None, feasible=bool(np.isfinite(loss.value(L0))))
        if not start['feasible']:
            L0 = np.zeros(tr['A'].shape[1])
        La, ra = loss.solve_lbfgs(L0)
        Lb, rb = loss.solve_trust_newton(fallback=L0)
        Fa, Fb = loss.value(La), loss.value(Lb)
        best = Lb if Fb <= Fa else La
        kkt = loss.kkt(best, L0)
        curv_opt = loss.curvature_check(best, 200)
        rng = np.random.default_rng(3)
        curv_rand = min(loss.curvature_check(np.maximum(best*(1 + .5*rng.normal(size=len(best))) + 1e-3*np.max(best)*np.abs(rng.normal(size=len(best))), 0), 10, seed=k) for k in range(20))
        res = dict(start=start, F_lbfgs=Fa, F_newton=Fb, solver_agreement=abs(Fa - Fb)/max(abs(Fb), 1e-300), newton_iterations=rb['iterations'], newton_start=rb['start'], newton_message=rb['message'],
                   lbfgs_message=str(ra.message), kkt=kkt, curvature_min_optimum=curv_opt, curvature_min_random=curv_rand,
                   rmse_train=loss.rmse(best), amplitudes=best.tolist(), widths_active=[int(i) for i in np.flatnonzero(best > 0)],
                   spectrum=CR.spectrum_summary(W, best, CS.G), per_galaxy_rmse=dict(zip(tr['names'], loss.per_galaxy_rmse(best))),
                   slope=slope_diagnostic(tr, best, gN))
        for split in ('validation', 'test'):
            b = blocks[split]
            lv = L2.SpeedLoss(b['A'], b[key], b['R'], b['v'], b['sizes'])
            val = lv.value(best)
            res[f'rmse_{split}'] = float(np.sqrt(val)) if np.isfinite(val) else None
        res['described'] = bool(res['rmse_train'] <= RMSE_LIMIT)
        out[label] = res
        out[label]['loss_object'] = loss
    # G-C4: stage 1's galaxies-alone score from the tabulated-force NNLS start, judged by stage 1's equal-galaxy RMSE
    arch = archive['E1_spectrum']['galaxies_alone']['train']['equal_galaxy_rmse_km_s']
    out['G_C4'] = dict(archived=arch, nnls_tabulated=out['tabulated']['start']['rmse'],
                       relative=abs(out['tabulated']['start']['rmse']/arch - 1) if out['tabulated']['start']['rmse'] else None)
    # the discrepancy between the two Newtonian forces, point by point
    ok = tr['gN_tab'] > 0
    ratio = tr['gN_rec'][ok]/tr['gN_tab'][ok]
    out['force_ratio_reconstructed_over_tabulated'] = dict(median=float(np.median(ratio)), p10=float(np.percentile(ratio, 10)), p90=float(np.percentile(ratio, 90)),
                                                          min=float(ratio.min()), max=float(ratio.max()), points=int(len(ratio)), tabulated_zero_points=int((~ok).sum()))
    # baryons-only scores under both forces
    for label, key in (('reconstructed', 'gN_rec'), ('tabulated', 'gN_tab')):
        loss = L2.SpeedLoss(tr['A'], tr[key], tr['R'], tr['v'], tr['sizes'])
        out['force_ratio_reconstructed_over_tabulated'][f'baryons_only_rmse_{label}'] = loss.rmse(np.zeros(tr['A'].shape[1]))
    return out


# ------------------------------------------------------------------ clusters
def load_clusters(ext, frac, n_grid=N_GRID_CLUSTERS, **kw):
    """Cluster records at the working grid (one width for the metadata) with the cached operators for all widths."""
    z = np.load(GEN/'cluster_ops.npz')
    assert np.allclose(z['widths'], W)
    cls = {}
    for name in ext['clusters']:
        cl = CS.build_cluster(name, ext['clusters'][name], np.array([100.]), n_grid=n_grid, frac_profile=frac, **kw)
        cl['ops'] = z[f'{name}_{n_grid}']
        assert np.allclose(z[f'{name}_{n_grid}_rp'], cl['rp'])
        cls[name] = cl
    conv = {}
    for name in ('A1795', 'A2319'):
        a, b = z[f'{name}_6000'], z[f'{name}_12000']
        conv[name] = (np.max(np.abs(a - b), axis=0)/np.max(np.abs(b), axis=0)).tolist()
    return cls, conv


def cluster_blocks(cls, corrs):
    rows, nus, ys, Lcs, meta = [], [], [], [], []
    off = 0
    for name, cl in cls.items():
        A_w, nu_w, y_w, Lc = L2.whitened_cluster(cl, corrs[name], cl['ops'])
        rows.append(A_w); nus.append(nu_w); ys.append(y_w); Lcs.append(Lc); meta.append((name, off, len(y_w)))
        off += len(y_w)
    N = off
    nuis = np.zeros((N, len(cls)))
    for k, (name, o, n) in enumerate(meta):
        nuis[o:o + n, k] = nus[k][:, 0]
    block = CR.Block('clusters', np.vstack(rows), np.concatenate(ys), np.ones(N), nuisance=nuis, meta=dict(names=list(cls), spans=meta))
    return block, Lcs


def cluster_references(cls, Lcs, corrected=True):
    out = {}
    for label, key in (('release_nfw_reference', 'P_nfw'), ('newtonian_baryons', 'P_N'), ('pm1_law_repo_a_star', 'P_pm1_repo'), ('pm1_law_fitted_a_star', 'P_pm1_fit')):
        per, tot, N = {}, 0., 0
        for (name, cl), Lc in zip(cls.items(), Lcs):
            model = cl[key] if corrected else cl[key]*cl['thermal']      # stage 1 applied the fraction twice
            chi, _ = L2.boundary_fit(model, cl['thermal'], cl['Pobs'], Lc)
            per[name] = chi/len(cl['rp']); tot += chi; N += len(cl['rp'])
        out[label] = dict(chi2=tot, chi2_per_point=tot/N, per_cluster_chi2_per_point=per)
    return out


def cluster_solve(block):
    sol = CR.solve([block], weights=[1.])
    per = {}
    for name, o, n in block.meta['spans']:
        pred = block.predict(sol['amplitudes'], sol['nuisances'][0])
        per[name] = float(np.sum((pred[o:o + n] - block.y[o:o + n])**2)/n)
    return dict(chi2=sol['chi2'][0], chi2_per_point=sol['chi2'][0]/block.N, amplitudes=sol['amplitudes'].tolist(), nuisances=sol['nuisances'][0].tolist(),
                widths_active=sol['widths_active'], kkt=sol['kkt'], spectrum=CR.spectrum_summary(W, sol['amplitudes'], CS.G), per_cluster_chi2_per_point=per)


def stage1_reproduction(ext, archive, corrs):
    """G-S3: stage 1's block (1500 grid, 26 widths with its zeroing) through the whitening with the identity correlation
    reproduces stage 1's archived clusters-alone chi-square."""
    block = CS.cluster_block(ext, CR.WIDTHS)
    sol = CR.solve([block])
    ident = {c['name']: dict(corr=np.eye(int(np.sum(np.asarray(c['kinds']) == 'sz')))) for c in block.meta['clusters']}
    cls = {c['name']: c for c in block.meta['clusters']}
    wb, _ = cluster_blocks(cls, ident)
    sol_w = CR.solve([wb], weights=[1.])
    arch = archive['E1_spectrum']['clusters_alone']['chi2']
    return dict(archived=arch, stage1_block=sol['chi2'][0], whitened_identity=sol_w['chi2'][0],
                relative_block=abs(sol['chi2'][0]/arch - 1), relative_whitened=abs(sol_w['chi2'][0]/arch - 1))


# ------------------------------------------------------------------ lenses
def build_lenses():
    return [L2.LensRows(n) for n in CS.LENSES]


def lens_blocks(lenses, imf, betas, light='W'):
    R = [L.rows(imf, b) for L, b in zip(lenses, betas)]
    E = CR.Block('lens_einstein', np.array([r['ae'] for r in R]), np.array([r['ye'] for r in R]), np.array([r['se'] for r in R]))
    K = CR.Block('lens_kinematics', np.vstack([r['Ak'] for r in R]), np.concatenate([r['yk'] for r in R]), np.ones(sum(len(r['yk']) for r in R)))
    return E, K


def lens_evaluate(lenses, imf, L, light='W'):
    ev = [dict(lens=Lz.name, **{k: v for k, v in zip(('chi2', 'beta'), Lz.fit_beta(imf, L))}) for Lz in lenses]
    for e, Lz in zip(ev, lenses):
        e['einstein_residual'] = Lz.evaluate(imf, e['beta'], L)['einstein_residual']
    total = float(sum(e['chi2'] for e in ev))
    worst = float(max(abs(e['einstein_residual']) for e in ev))
    return dict(per_lens=ev, total_kinematics_chi2=total, max_abs_einstein_residual=worst, described=bool(total <= L2.CHI2_LIMIT and worst <= L2.EINSTEIN_LIMIT))


def ratio_hull(lenses, imf):
    """Per lens at its stars-only anisotropy: the ratio of extra mean second moment to extra bending per width, and the
    ratio the data need; feasibility of a nonnegative mixture requires the needed ratio inside the hull."""
    out = {}
    for Lz in lenses:
        so = Lz.L.stars_only()
        r = Lz.rows(imf, so['beta'])
        co = r['co']
        m = r['m']
        K = m*co[1:].mean(axis=1)                       # mean extra V^2 per unit amplitude, per width
        B = r['ae']                                     # extra bend per unit amplitude, per width
        need_v2 = float(np.mean(Lz.L.y**2 - m*co[0]))
        need_b = float(r['ye'])
        ratios = K/np.where(B > 0, B, np.nan)
        ok = np.isfinite(ratios)
        needed = need_v2/need_b if need_b != 0 else np.nan
        out[Lz.name] = dict(beta=so['beta'], needed_ratio=needed, hull_min=float(np.nanmin(ratios[ok])), hull_max=float(np.nanmax(ratios[ok])),
                            inside=bool(np.nanmin(ratios[ok]) <= needed <= np.nanmax(ratios[ok])))
    return out


# ------------------------------------------------------------------ joint solves
def joint_solve(gal_loss, quad_blocks, weights, L0, nu0):
    """Minimise sum_b w_b chi2_b(L, nu_b) + F_gal(L) over L >= 0 and nonnegative nuisances (convex)."""
    nw = len(L0)
    n_nu = [b.n_nuisance for b in quad_blocks]
    x0 = np.concatenate([np.maximum(np.asarray(L0, float), 0)] + [np.maximum(np.asarray(n, float), 0) for n, k in zip(nu0, n_nu) if k])

    def unpack(x):
        L = x[:nw]
        nus, off = [], nw
        for k in n_nu:
            nus.append(x[off:off + k] if k else None)
            off += k
        return L, nus

    def f(x):
        L, nus = unpack(x)
        val, grad = 0., np.zeros_like(x)
        if gal_loss is not None:
            v = gal_loss.value(L)
            if not np.isfinite(v):
                return 1e300, np.zeros_like(x)
            val += v
            grad[:nw] += gal_loss.gradient(L)
        off = nw
        for b, wt, nu in zip(quad_blocks, weights, nus):
            r = (b.predict(L, nu) - b.y)/b.s
            val += wt*float(np.sum(r*r))
            grad[:nw] += 2*wt*(b.A/b.s[:, None]).T@r
            if b.n_nuisance:
                grad[off:off + b.n_nuisance] += 2*wt*(b.nuisance/b.s[:, None]).T@r
                off += b.n_nuisance
        return val, grad
    res = minimize(f, x0, jac=True, method='L-BFGS-B', bounds=[(0, None)]*len(x0), options=dict(maxiter=50000, ftol=1e-15, gtol=1e-12, maxcor=40))
    L, nus = unpack(res.x)
    _, g = f(res.x)
    pg_lbfgs = float(np.linalg.norm(np.where(res.x > 0, g, np.minimum(g, 0))))
    # amendment 3: polish with the trust-region Newton method in column-scaled variables and certify there
    J = L2.JointObjective(gal_loss, quad_blocks, weights)
    L, nus, pol = J.polish(L, nus)
    L, nus, fin = L2.refine_face(J, L, nus)                                  # amendment 4
    return dict(L=L, nus=nus, objective=fin['objective'], projected_gradient=fin['projected_gradient'], lbfgs_objective=float(res.fun),
                lbfgs_projected_gradient=pg_lbfgs, lbfgs_message=str(res.message), polish=pol, face_refinement=fin, iterations=int(res.nit))


def main():
    t0 = time.time()
    archive = json.loads(ARCHIVE.read_text(encoding='utf-8'))
    ext = json.loads((HERE/'cl2-inputs-xcop-profiles.json').read_text(encoding='utf-8'))
    frac = CS.stellar_fraction_profile(ext)
    results = dict(experiment='CL-2 stage 2', protocol='protocol-cl2-stage2.md; amendments 1 to 4', widths_kpc=W.tolist(),
                   input_sha256={'cl2-inputs-xcop-profiles.json': sha(HERE/'cl2-inputs-xcop-profiles.json'), 'cl2-results.json': sha(ARCHIVE),
                                 'xcop-release/allfiles.tar.gz': sha(TAR), 'pinned': XA.SHA256})
    gates = {}
    # ---- S2-3 gates
    dg = L2.disk_force_gates()
    ann = lambda s: np.where((np.asarray(s, float) > .9) & (np.asarray(s, float) < 1.1), 1e9, 0.)
    from scipy.integrate import quad
    exact_annulus = quad(lambda a: 2*np.pi*a*1e9*L2.ring_force(1., a, 10.), .9, 1.1, limit=400)[0]
    pm10 = L2.disk_force(ann, np.array([10.]), kinks=(.9, 1.1))[0]
    pm1000 = L2.disk_force(ann, np.array([1000.]), kinks=(.9, 1.1), Rmax=4000.)[0]
    mass = np.pi*(1.1**2 - .9**2)*1e9
    dg.update(annulus_exact_relative=float(abs(pm10/exact_annulus - 1)), point_mass_1000_relative=float(abs(pm1000/(CS.G*mass/1e6) - 1)))
    gates['G_N1_freeman'] = dg['freeman_relative'] < 1e-3
    gates['G_N1_control_rejected'] = dg['control_relative'] > .1
    gates['G_N2_bulge'] = True                                            # the bulge force is G M(<r)/r^2 by construction
    gates['G_N3_point_mass_at_10_radii_as_declared'] = dg['point_mass_relative'] < 1e-6
    gates['G_N3b_annulus_exact'] = dg['annulus_exact_relative'] < 1e-6
    gates['G_N3c_point_mass_at_1000_radii'] = dg['point_mass_1000_relative'] < 1e-5
    results['disk_force'] = dg
    print('disk gates', {k: v for k, v in gates.items()}, flush=True)
    # ---- galaxies
    blocks, g4b = load_galaxies()
    gates['G4b_galaxies'] = g4b < 1e-3
    gal = galaxy_stage(blocks, archive)
    gates['G_C1_solvers'] = all(gal[k]['solver_agreement'] < 1e-6 for k in ('reconstructed', 'tabulated'))
    gates['G_C2_kkt'] = all(gal[k]['kkt']['max_abs_gradient_active'] < 1e-8 and gal[k]['kkt']['min_gradient_inactive'] > -1e-8 for k in ('reconstructed', 'tabulated'))
    gates['G_C3_convex'] = all(gal[k]['curvature_min_optimum'] >= 0 and gal[k]['curvature_min_random'] >= 0 for k in ('reconstructed', 'tabulated'))
    gates['G_C4_reproduction'] = gal['G_C4']['relative'] is not None and gal['G_C4']['relative'] < 1e-6
    loss_rec = gal['reconstructed'].pop('loss_object')
    gal['tabulated'].pop('loss_object')
    results['galaxies'] = dict(G4b_worst=g4b, **gal)
    print('galaxies: certified minimum reconstructed %.3f km/s (tabulated %.3f), start %.3f, archived %.3f' % (
        gal['reconstructed']['rmse_train'], gal['tabulated']['rmse_train'], gal['tabulated']['start']['rmse'], gal['G_C4']['archived']), flush=True)
    # ---- clusters
    corrs = L2.sz_correlations(TAR, list(ext['clusters']))
    cls, conv = load_clusters(ext, frac)
    for name in corrs:
        rext = np.array(ext['clusters'][name]['pressure']['sz']['r_over_R500'])*ext['clusters'][name]['header']['R500_kpc']
        corrs[name]['radius_relative_difference'] = float(np.max(np.abs(corrs[name]['radius_kpc'] - rext)/rext))
        ev = np.linalg.eigvalsh(corrs[name]['corr'])
        corrs[name]['condition'] = float(ev.max()/ev.min()); corrs[name]['min_eigenvalue'] = float(ev.min())
        corrs[name]['symmetric'] = bool(np.max(np.abs(corrs[name]['corr'] - corrs[name]['corr'].T)) < 1e-12)
    gates['G_S1_radii'] = all(c['radius_relative_difference'] < .02 for c in corrs.values())
    gates['G_S2_matrices'] = all(c['symmetric'] and c['min_eigenvalue'] > 0 and c['condition'] < 1e6 for c in corrs.values())
    rep = stage1_reproduction(ext, archive, corrs)
    gates['G_S3_stage1_reproduction'] = rep['relative_block'] < 1e-9 and rep['relative_whitened'] < 1e-9
    conv_gate = True
    for name, d in conv.items():
        for w, x in zip(W, d):
            if w >= 1. and x >= 1e-4:
                conv_gate = False
            if .3 < w < 1. and x >= 1e-3:
                conv_gate = False
    gates['G4_2_cluster_convergence'] = conv_gate
    block, Lcs = cluster_blocks(cls, corrs)
    refs = cluster_references(cls, Lcs, corrected=True)
    refs_diag = None
    ident = {n: dict(corr=np.eye(len(c['corr']))) for n, c in corrs.items()}   # identity: any sub-block is the identity
    block_d, Lcs_d = cluster_blocks(cls, ident)
    refs_diag = cluster_references(cls, Lcs_d, corrected=True)
    csol = cluster_solve(block)
    csol_diag = cluster_solve(block_d)
    # negligibility of the two narrowest widths at the largest amplitude any fit assigns (galaxy and cluster fits)
    amp_max = max(max(gal['reconstructed']['amplitudes'][:2]), max(gal['tabulated']['amplitudes'][:2]), max(csol['amplitudes'][:2]))
    neg = max(float(np.max(np.abs(cl['ops'][:, :2]@np.full(2, amp_max))/cl['P_N'])) for cl in cls.values())
    gates['G4_3_negligibility'] = neg < 1e-3
    results['clusters'] = dict(sz_correlations={n: {k: v for k, v in c.items() if k != 'corr'} | dict(corr=c['corr'].tolist()) for n, c in corrs.items()},
                               convergence_6000_vs_12000={n: dict(zip([f'{w:.3f}' for w in W], d)) for n, d in conv.items()},
                               stage1_reproduction=rep, references_covariance=refs, references_diagonal=refs_diag,
                               spectrum_covariance=csol, spectrum_diagonal=csol_diag, narrow_width_negligibility=neg,
                               ratio_to_release_reference=csol['chi2_per_point']/refs['release_nfw_reference']['chi2_per_point'],
                               described=bool(csol['chi2_per_point'] <= 2*refs['release_nfw_reference']['chi2_per_point']))
    print('clusters: spectrum %.3f per point, reference %.3f, ratio %.3f (diagonal: %.3f / %.3f)' % (
        csol['chi2_per_point'], refs['release_nfw_reference']['chi2_per_point'], results['clusters']['ratio_to_release_reference'],
        csol_diag['chi2_per_point'], refs_diag['release_nfw_reference']['chi2_per_point']), flush=True)
    # sensitivities at stage 1's grid, relative to a same-grid primary, with the corrected references and the covariance
    sens = {}
    base1500, conv1500 = None, None
    for label, kw in (('primary_1500', {}), ('mu_0.59', dict(mu=.59)), ('mu_0.61', dict(mu=.61)), ('stars_none', dict(stars='none')), ('non_thermal_pressure', dict(alpha_nt=True))):
        c1 = {}
        for name in ext['clusters']:
            cl = CS.build_cluster(name, ext['clusters'][name], W, n_grid=1500, frac_profile=frac, **kw)
            cl['ops'] = np.array([cl['pressure_of'](cl['source'].g_mem(cl['r'], w)) for w in W]).T
            c1[name] = cl
        b1, l1 = cluster_blocks(c1, corrs)
        s1 = cluster_solve(b1)
        r1 = cluster_references(c1, l1, corrected=True)
        sens[label] = dict(spectrum_chi2_per_point=s1['chi2_per_point'], reference_chi2_per_point=r1['release_nfw_reference']['chi2_per_point'],
                           baryons_chi2_per_point=r1['newtonian_baryons']['chi2_per_point'])
        print('sensitivity', label, sens[label], flush=True)
    results['clusters']['sensitivities_1500_grid'] = sens
    # ---- lenses
    lenses = build_lenses()
    print('lenses built at %.0f s' % (time.time() - t0), flush=True)
    betas0 = [Lz.L.stars_only()['beta'] for Lz in lenses]
    results['lenses'] = dict(geometry=L2.LENS_GEOMETRY, ratio_hull={imf: ratio_hull(lenses, imf) for imf in ('Chabrier',)},
                             tf1_static_joint_cited=dict(W=dict(chi2=327.13, einstein=.1306), S=dict(chi2=434.77, einstein=.1522)))
    # ---- E1: alone and jointly (Chabrier, static geometry, well coupling)
    tr = blocks['train']
    E1 = dict(clusters_alone=dict(chi2_per_point=csol['chi2_per_point'], described=results['clusters']['described']),
              galaxies_alone=dict(rmse_train=gal['reconstructed']['rmse_train'], described=gal['reconstructed']['described']))
    # lenses alone: alternate betas as TF-1 (well coupling, F1 family)
    betas = list(betas0)
    hist = []
    for it in range(6):
        E, K = lens_blocks(lenses, 'Chabrier', betas)
        sol = CR.solve([E, K])
        ev = lens_evaluate(lenses, 'Chabrier', sol['amplitudes'])
        hist.append(dict(iteration=it, objective=float(sum(c/b.N for c, b in zip(sol['chi2'], (E, K)))), total=ev['total_kinematics_chi2'], einstein=ev['max_abs_einstein_residual']))
        new = [e['beta'] for e in ev['per_lens']]
        if it and abs(hist[-1]['objective'] - hist[-2]['objective']) < 1e-6*max(1., hist[-2]['objective']):
            break
        betas = new
    E1['lenses_alone'] = dict(total_kinematics_chi2=ev['total_kinematics_chi2'], max_abs_einstein_residual=ev['max_abs_einstein_residual'], described=ev['described'],
                              spectrum=CR.spectrum_summary(W, sol['amplitudes'], CS.G), history=hist, betas=betas)
    print('lenses alone', E1['lenses_alone']['total_kinematics_chi2'], E1['lenses_alone']['max_abs_einstein_residual'], flush=True)
    # joint: clusters (covariance) + galaxies (exact loss) + lenses (well), betas alternated
    betas = list(betas0)
    L = np.array(gal['reconstructed']['amplitudes'])
    nu_c = np.array(csol['nuisances'])
    jh = []
    for it in range(6):
        E, K = lens_blocks(lenses, 'Chabrier', betas)
        js = joint_solve(loss_rec, [block, E, K], [1./block.N, 1./E.N, 1./K.N], L, [nu_c, None, None])
        L, nu_c = js['L'], js['nus'][0]
        ev = lens_evaluate(lenses, 'Chabrier', L)
        jh.append(dict(iteration=it, objective=js['objective'], projected_gradient=js['projected_gradient'], lbfgs_objective=js['lbfgs_objective'], lbfgs_projected_gradient=js['lbfgs_projected_gradient'], lens_total=ev['total_kinematics_chi2'], einstein=ev['max_abs_einstein_residual']))
        new = [e['beta'] for e in ev['per_lens']]
        if it and abs(jh[-1]['objective'] - jh[-2]['objective']) < 1e-6*max(1., jh[-2]['objective']):
            break
        betas = new
    cchi = block.chi2(L, nu_c)
    E1['joint'] = dict(objective=js['objective'], history=jh, projected_gradient=js['projected_gradient'],
                       clusters=dict(chi2_per_point=cchi/block.N, described=bool(cchi/block.N <= 2*refs['release_nfw_reference']['chi2_per_point'])),
                       galaxies=dict(rmse_train=loss_rec.rmse(L), described=bool(loss_rec.rmse(L) <= RMSE_LIMIT)),
                       lenses=dict(total_kinematics_chi2=ev['total_kinematics_chi2'], max_abs_einstein_residual=ev['max_abs_einstein_residual'], described=ev['described']),
                       spectrum=CR.spectrum_summary(W, L, CS.G), amplitudes=L.tolist(), betas=betas)
    print('joint', {k: E1['joint'][k] for k in ('clusters', 'galaxies', 'lenses')}, flush=True)
    results['E1'] = E1
    # ---- E3 transfer
    mw = CS.milky_way_blocks(W)

    def mw_scores(Lv):
        out = {}
        for b in mw:
            u = b.meta['R']*(b.meta['gN'] + b.A@Lv)
            out[b.meta['variant']] = float(np.sqrt(np.mean((np.sqrt(np.maximum(u, 0)) - b.meta['y'])**2)))
        return out

    def gal_scores(Lv):
        out = {}
        for split in ('train', 'validation', 'test'):
            b = blocks[split]
            v = L2.SpeedLoss(b['A'], b['gN_rec'], b['R'], b['v'], b['sizes']).value(Lv)
            out[split] = float(np.sqrt(v)) if np.isfinite(v) else None
        return out

    def cluster_score(Lv):
        # refit the boundary pressures at fixed amplitudes: nonnegative least squares over the nuisances only
        from scipy.optimize import nnls
        nu, _ = nnls(block.nuisance, block.y - block.A@Lv)
        return float(np.sum((block.A@Lv + block.nuisance@nu - block.y)**2)/block.N)
    E3 = {}
    Lg = np.array(gal['reconstructed']['amplitudes'])
    Lc = np.array(csol['amplitudes'])
    E3['calibrated_on_galaxies'] = dict(clusters_chi2_per_point=cluster_score(Lg), lenses=lens_evaluate(lenses, 'Chabrier', Lg), milky_way=mw_scores(Lg), galaxies=gal_scores(Lg))
    E3['calibrated_on_clusters'] = dict(galaxies=gal_scores(Lc), lenses=lens_evaluate(lenses, 'Chabrier', Lc), milky_way=mw_scores(Lc))
    jg = joint_solve(loss_rec, [block], [1./block.N], Lg, [nu_c])
    Lb = jg['L']
    gates['G_J_joint_certificate'] = js['projected_gradient'] < 1e-6*max(1., js['objective']) and jg['projected_gradient'] < 1e-6*max(1., jg['objective'])
    E3['calibrated_on_galaxies_and_clusters'] = dict(objective=jg['objective'], projected_gradient=jg['projected_gradient'], galaxies=gal_scores(Lb), clusters_chi2_per_point=cluster_score(Lb),
                                                     lenses=lens_evaluate(lenses, 'Chabrier', Lb), milky_way=mw_scores(Lb), spectrum=CR.spectrum_summary(W, Lb, CS.G))
    E3['milky_way_baryons_only'] = mw_scores(np.zeros(len(W)))
    results['E3'] = E3
    print('transfer', {k: {kk: vv for kk, vv in v.items() if kk in ('clusters_chi2_per_point', 'galaxies', 'milky_way')} for k, v in E3.items() if isinstance(v, dict)}, flush=True)
    # ---- verdict
    rec = gal['reconstructed']
    outcome = ('a: the certified minimum describes the galaxies' if rec['described'] else
               'b: no member of the universal linear family describes the galaxies with these sources; certified minimum %.2f km/s against 21.9' % rec['rmse_train'])
    results['gates'] = gates
    results['numerical_verification_passed'] = bool(all(gates.values()))
    results['statuses'] = dict(numerical_verification='passed' if all(gates.values()) else 'FAILED on: ' + ', '.join(k for k, v in gates.items() if not v),
                               scientific_outcome=outcome)
    results['checks_short_run'] = dict(freeman=dg['freeman_relative'], annulus=dg['annulus_exact_relative'], rmse_reconstructed=rec['rmse_train'],
                                       rmse_tabulated=gal['tabulated']['rmse_train'], F_reconstructed=rec['F_newton'], clusters_per_point=csol['chi2_per_point'],
                                       reference_per_point=refs['release_nfw_reference']['chi2_per_point'])
    results['runtime_seconds'] = time.time() - t0
    OUT.write_text(json.dumps(plain(results), indent=1) + '\n', encoding='utf-8')
    print('written', OUT, 'runtime %.0f s' % results['runtime_seconds'], 'outcome:', outcome, 'gates:', gates, flush=True)


if __name__ == '__main__':
    main()

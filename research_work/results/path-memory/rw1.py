"""RW-1 driver (protocol-rw1.md, amendment 1): the whirlpool beyond the rut, backed into regime by regime. Reads the
caches of rw1_build.py, stage 2's caches and archive, the X-COP extract and correlations, CR-2's lens interface and the
Coma bins; writes rw1-results.json."""
import hashlib
import json
import os
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl2_sources as CS      # noqa: E402
import cl2s2_lib as L2        # noqa: E402
import cl2s2 as D             # noqa: E402
import rw1_lib as RW          # noqa: E402
import steady_field as SF     # noqa: E402
import cl1 as C1              # noqa: E402
import cluster_lensing as CLZ  # noqa: E402

ROOT = HERE.parents[2]
GEN = ROOT/'research_work/generated/routes'
OUT = HERE/'rw1-results.json'
V_ACCEPT = 21.9
G = SF.G
REGIMES = ('inner', 'outer', 'edge', 'whole')
BS_KPC = (100., 300., 1000.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plain(x):
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return plain(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def rmse_of(gals, lam_by_gal, k, source='rec', regime='whole'):
    """Equal-galaxy speed error over the galaxies that carry the regime, with per-galaxy strengths on column k."""
    vals = []
    for g, lam in zip(gals, lam_by_gal):
        m = g['masks'][regime]
        if m.sum() < 2:
            continue
        u = g['R'][m]*(g['gN_' + source][m] + lam*g['cols'][m, k])
        if np.any(u <= 0):
            return None
        vals.append(np.mean((np.sqrt(u) - g['v'][m])**2))
    return float(np.sqrt(np.mean(vals)))


def kernel_index(kernels, p, w):
    for i, (pp, ww) in enumerate(kernels):
        if abs(pp - p) < 1e-12 and abs(ww - w) < 1e-12:
            return i
    return None


def main():
    t0 = time.time()
    arch2 = json.loads((HERE/'cl2s2-results.json').read_text(encoding='utf-8'))
    arch1 = json.loads((HERE/'cl2-results.json').read_text(encoding='utf-8'))
    archn = json.loads((HERE/'nl1-results.json').read_text(encoding='utf-8'))
    gates, out = {}, dict(experiment='RW-1', protocol='protocol-rw1.md; amendment 1', kernels_galaxies=RW.KERNELS_GAL, kernels_clusters=RW.KERNELS_CL,
                         input_sha256={f: sha(HERE/f) for f in ('cl2s2-results.json', 'cl2-results.json', 'nl1-results.json', 'cl2-inputs-xcop-profiles.json')})
    for f in ('rw1_galaxies.npz', 'rw1_clusters.npz', 'rw1_milky_way.npz'):
        if (GEN/f).exists() or not os.environ.get('RW1_STOP_AFTER_GALAXIES'):
            out['input_sha256'][f] = sha(GEN/f)
    # ---------------------------------------------------------------- K1-K6
    kg = RW.kernel_gates()
    out['kernel_gates'] = kg
    for k, v in kg['passes'].items():
        gates[k] = bool(v)
    print('kernel gates', kg['passes'], flush=True)
    # ---------------------------------------------------------------- galaxies
    blocks, _ = D.load_galaxies()
    z = np.load(GEN/'rw1_galaxies.npz')
    assert np.allclose(z['kernels'], np.array(RW.KERNELS_GAL))
    KG = RW.KERNELS_GAL
    nk = len(KG)
    gals_by_split = {}
    for split in ('train', 'validation', 'test'):
        b = blocks[split]
        gals, off = [], 0
        for i, name in enumerate(b['names']):
            n = b['sizes'][i]
            key = f'{split}/{name}/'
            R = z[key + 'R']
            assert np.allclose(R, b['R'][off:off + n]), name
            geo = z[key + 'geometry']
            gals.append(dict(name=name, R=R, v=b['v'][off:off + n], gN_rec=b['gN_rec'][off:off + n], gN_tab=b['gN_tab'][off:off + n],
                             cols=z[key + 'cols'], point=z[key + 'point'], R_half=float(geo[0]), M_grid=float(geo[1]), rd=float(geo[2]), M=float(geo[3]),
                             masks=RW.regime_masks(R, float(geo[0]))))
            off += n
        gals_by_split[split] = gals
    train = gals_by_split['train']
    # K7
    k7 = {}
    for name in ('NGC2403', 'NGC3198', 'DDO154'):
        key = next(k for k in z.files if k.endswith(f'/{name}/nodes256'))[:-len('nodes256')]      # the galaxy's own split
        c, c256, c6000 = z[key + 'cols'], z[key + 'nodes256'], z[key + 'grid6000']
        rel = lambda a, b: float(np.max(np.max(np.abs(a - b), axis=0)/np.max(np.abs(b), axis=0)))
        k7[name] = dict(nodes=rel(c, c256), grid=rel(c, c6000))
    out['K7'] = k7
    gates['K7_nodes'] = all(v['nodes'] < 1e-8 for v in k7.values())
    gates['K7_grid'] = all(v['grid'] < 1e-3 for v in k7.values())
    # R1: lambda = 0 with the tabulated force against stage 1's archived baryons reference
    r1 = float(np.sqrt(np.mean([np.mean((np.sqrt(g['R']*g['gN_tab']) - g['v'])**2) for g in train])))
    ref1 = arch1['E1_spectrum']['references']['galaxies_train']['baryons']['equal_galaxy_rmse_km_s']
    out['R1'] = dict(baryons_tabulated=r1, archived=ref1, relative=abs(r1/ref1 - 1))
    gates['R1_galaxies'] = out['R1']['relative'] < 1e-9
    print('R1 %.6f vs %.6f' % (r1, ref1), flush=True)
    # F1: individual fits per galaxy, regime, source, kernel
    excluded = {r: [g['name'] for g in train if g['masks'][r].sum() < 2] for r in REGIMES}
    out['excluded_from_regime'] = excluded
    F1 = {}
    worst_agree, worst_dF = 0., 0.
    for source in ('rec', 'tab'):
        for regime in REGIMES:
            Fm, Lm = np.full((len(train), nk), np.nan), np.full((len(train), nk), np.nan)
            for i, g in enumerate(train):
                m = g['masks'][regime]
                if m.sum() < 2:
                    continue
                wts = np.full(int(m.sum()), 1./m.sum())
                for k in range(nk):
                    f = RW.fit_lambda(g['cols'][m, k], g['gN_' + source][m], g['R'][m], g['v'][m], wts)
                    Fm[i, k], Lm[i, k] = f['F'], f['lam']
                    worst_agree = max(worst_agree, f['agree'])
                    if f['how'] == 'interior':
                        worst_dF = max(worst_dF, abs(f['dF'])/abs(f['dF0']))
            F1[(source, regime)] = (Fm, Lm)
            print('F1 %s %s done at %.0f s' % (source, regime, time.time() - t0), flush=True)
    # N2: the point-mass whirlpool, whole regime, reconstructed source
    Fp, Lp = np.full((len(train), nk), np.nan), np.full((len(train), nk), np.nan)
    for i, g in enumerate(train):
        wts = np.full(len(g['R']), 1./len(g['R']))
        for k in range(nk):
            f = RW.fit_lambda(g['point'][:, k], g['gN_rec'], g['R'], g['v'], wts)
            Fp[i, k], Lp[i, k] = f['F'], f['lam']
            worst_agree = max(worst_agree, f['agree'])
    out['C1_worst_agreement'] = worst_agree
    out['C1_worst_derivative_ratio'] = worst_dF
    # F2 and F3 per regime and source; S1; S2
    galaxies = {}
    for source in ('rec', 'tab'):
        for regime in REGIMES:
            Fm, Lm = F1[(source, regime)]
            carry = ~np.isnan(Fm[:, 0])
            mean_F = np.nanmean(Fm, axis=0)
            k2 = int(np.argmin(mean_F))
            best_k = np.array([int(np.nanargmin(Fm[i])) if carry[i] else -1 for i in range(len(train))])
            # F3 universal
            idx = np.flatnonzero(carry)
            W = np.concatenate([train[i]['cols'][train[i]['masks'][regime]] for i in idx])
            gN = np.concatenate([train[i]['gN_' + source][train[i]['masks'][regime]] for i in idx])
            Rr = np.concatenate([train[i]['R'][train[i]['masks'][regime]] for i in idx])
            vv = np.concatenate([train[i]['v'][train[i]['masks'][regime]] for i in idx])
            wts = np.concatenate([np.full(int(train[i]['masks'][regime].sum()), 1./(train[i]['masks'][regime].sum()*len(idx))) for i in idx])
            F3 = [RW.fit_lambda(W[:, k], gN, Rr, vv, wts) for k in range(nk)]
            f3 = np.array([f['F'] for f in F3])
            k3 = int(np.argmin(f3))
            worst_agree = max(worst_agree, max(f['agree'] for f in F3))
            # C2 separability at k2
            lam2 = Lm[:, k2]
            assembled = rmse_of(train, np.nan_to_num(lam2), k2, source, regime)
            c2 = abs(assembled**2 - mean_F[k2])/mean_F[k2]
            # N1: the Newtonian kernel
            p2 = [k for k, (p, w) in enumerate(KG) if p == 2.]
            k2_newton = p2[int(np.argmin(mean_F[p2]))]
            k3_newton = p2[int(np.argmin(f3[p2]))]
            # N3
            n3 = float(np.sqrt(np.mean([np.mean((np.sqrt(train[i]['R'][train[i]['masks'][regime]]*train[i]['gN_' + source][train[i]['masks'][regime]]) - train[i]['v'][train[i]['masks'][regime]])**2) for i in idx])))
            # S1, S2
            hist_p = {str(p): int(sum(1 for i in idx if KG[best_k[i]][0] == p)) for p in RW.P_GRID}
            hist_w = {str(w): int(sum(1 for i in idx if KG[best_k[i]][1] == w)) for w in RW.W_GAL}
            M = np.array([train[i]['M'] for i in idx])
            Rh = np.array([train[i]['R_half'] for i in idx])
            lam_h = lam2[idx]*(Rh/RW.R_STAR)**(2 - KG[k2][0])
            entry = dict(n_galaxies=int(carry.sum()), individual=dict(best_kernel=[(KG[best_k[i]] if carry[i] else None) for i in range(len(train))],
                                                                     best_rmse=[float(np.sqrt(Fm[i, best_k[i]])) if carry[i] else None for i in range(len(train))],
                                                                     lambda_at_best=[float(Lm[i, best_k[i]]) if carry[i] else None for i in range(len(train))],
                                                                     hist_p=hist_p, hist_w=hist_w, equal_galaxy_rmse_at_own_best=float(np.sqrt(np.mean([Fm[i, best_k[i]] for i in idx])))),
                         F2=dict(kernel=KG[k2], rmse=float(np.sqrt(mean_F[k2])), lambda_by_galaxy=[float(x) if not np.isnan(x) else None for x in lam2],
                                 rmse_by_kernel=[float(np.sqrt(x)) for x in mean_F], C2_relative=float(c2), matched=bool(np.sqrt(mean_F[k2]) <= V_ACCEPT),
                                 regression=RW.regress_loglog(M, lam2[idx]), regression_h=RW.regress_loglog(M, lam_h),
                                 newton_control=dict(kernel=KG[k2_newton], rmse=float(np.sqrt(mean_F[k2_newton])))),
                         F3=dict(kernel=KG[k3], lam=F3[k3]['lam'], rmse=float(np.sqrt(f3[k3])), rmse_by_kernel=[float(np.sqrt(x)) for x in f3],
                                 lambda_by_kernel=[f['lam'] for f in F3], how=F3[k3]['how'], matched=bool(np.sqrt(f3[k3]) <= V_ACCEPT),
                                 newton_control=dict(kernel=KG[k3_newton], lam=F3[k3_newton]['lam'], rmse=float(np.sqrt(f3[k3_newton])))),
                         baryons_alone_rmse=n3)
            galaxies[f'{source}_{regime}'] = entry
            print('%s %s: F2 %s %.2f km/s, F3 %s lam %.4g %.2f km/s, baryons %.2f, C2 %.1e' % (source, regime, KG[k2], np.sqrt(mean_F[k2]), KG[k3], F3[k3]['lam'], np.sqrt(f3[k3]), n3, c2), flush=True)
    meanFp = np.nanmean(Fp, axis=0)
    kp = int(np.argmin(meanFp))
    galaxies['N2_point_source_whole_rec'] = dict(kernel=KG[kp], rmse=float(np.sqrt(meanFp[kp])), rmse_by_kernel=[float(np.sqrt(x)) for x in meanFp],
                                                 matter_sourced_rmse=galaxies['rec_whole']['F2']['rmse'], difference_km_s=float(np.sqrt(meanFp[kp]) - galaxies['rec_whole']['F2']['rmse']))
    gates['C1_solvers'] = worst_agree < 1e-9
    gates['C1_derivative'] = worst_dF < 1e-8
    gates['C2_separability'] = all(galaxies[k]['F2']['C2_relative'] < 1e-12 for k in galaxies if 'F2' in galaxies[k])
    out['galaxies'] = galaxies
    out['galaxy_table'] = [dict(name=g['name'], M=g['M'], R_half=g['R_half'], rd=g['rd'], n=len(g['R']), n_inner=int(g['masks']['inner'].sum())) for g in train]
    if os.environ.get('RW1_STOP_AFTER_GALAXIES'):                      # a plumbing test of the galaxy section only; never the archived run
        Path(os.environ['RW1_STOP_AFTER_GALAXIES']).write_text(json.dumps(plain(out), indent=1), encoding='utf-8')
        print('stopped after the galaxies (test mode) at %.0f s' % (time.time() - t0), flush=True)
        return
    # ---------------------------------------------------------------- clusters
    ext = CS.load_xcop()
    frac = CS.stellar_fraction_profile(ext)
    names = list(ext['clusters'])
    corrs = L2.sz_correlations(D.TAR, names)
    cls, _ = D.load_clusters(ext, frac)
    zc = np.load(GEN/'rw1_clusters.npz')
    assert np.allclose(zc['kernels'], np.array(RW.KERNELS_CL))
    KC = RW.KERNELS_CL
    nkc = len(KC)
    rows, Lcs = {}, []
    for name in names:
        cl = cls[name]
        assert np.allclose(zc[name + '/rp'], cl['rp'])
        A_w, nu_w, y_w, Lc = L2.whitened_cluster(cl, corrs[name], zc[name + '/ops'])
        rows[name] = (A_w, nu_w, y_w)
        Lcs.append(Lc)
    refs = D.cluster_references(cls, Lcs)
    N_c = sum(len(cls[n]['rp']) for n in names)
    ref_release = refs['release_nfw_reference']['chi2_per_point']
    out['R2'] = dict(release_reference=ref_release, newtonian_baryons=refs['newtonian_baryons']['chi2_per_point'],
                     archived_release=arch2['clusters']['references_covariance']['release_nfw_reference']['chi2_per_point'],
                     archived_baryons=arch2['clusters']['references_covariance']['newtonian_baryons']['chi2_per_point'])
    out['R2']['relative'] = max(abs(out['R2']['release_reference']/out['R2']['archived_release'] - 1), abs(out['R2']['newtonian_baryons']/out['R2']['archived_baryons'] - 1))
    gates['R2_clusters'] = out['R2']['relative'] < 1e-9
    k8 = {}
    for name in ('A1795', 'A2255'):
        a, f, g12 = zc[name + '/K8_interp'], zc[name + '/K8_full'], zc[name + '/K8_grid12000']
        rel = lambda x, y: float(np.max(np.max(np.abs(x - y), axis=0)/np.max(np.abs(y), axis=0)))
        k8[name] = dict(interp=rel(a, f), grid=rel(a, g12))
    out['K8'] = k8
    gates['K8_interp'] = all(v['interp'] < 1e-5 for v in k8.values())
    gates['K8_grid'] = all(v['grid'] < 1e-3 for v in k8.values())
    chi = np.zeros((len(names), nkc))
    lamc = np.zeros((len(names), nkc))
    kkt_worst = 0.
    for i, name in enumerate(names):
        A_w, nu_w, y_w = rows[name]
        for k in range(nkc):
            f = RW.cluster_fit(A_w[:, k], nu_w, y_w)
            chi[i, k], lamc[i, k] = f['chi2'], f['lam']
            kkt_worst = max(kkt_worst, f['kkt'])
    tot = chi.sum(axis=0)/N_c
    kc2 = int(np.argmin(tot))
    uni = [RW.cluster_fit_universal([rows[n][0][:, k] for n in names], [rows[n][1] for n in names], [rows[n][2] for n in names]) for k in range(nkc)]
    kkt_worst = max(kkt_worst, max(u['kkt'] for u in uni))
    tot3 = np.array([u['chi2'] for u in uni])/N_c
    kc3 = int(np.argmin(tot3))
    p2c = [k for k, (p, w) in enumerate(KC) if p == 2.]
    cl_geo = {}
    for name in names:
        cl = cls[name]
        Mb = cl['M_b']
        cl_geo[name] = dict(M_b=float(Mb[-1]), R_half=float(np.interp(.5*Mb[-1], Mb, cl['r'])), R500=float(cl['R500']), r_out=float(cl['r_out']), n=len(cl['rp']))
    Mc = np.array([cl_geo[n]['M_b'] for n in names])
    Rhc = np.array([cl_geo[n]['R_half'] for n in names])
    best_kc = [int(np.argmin(chi[i])) for i in range(len(names))]
    clusters = dict(per_point_reference=ref_release, N=N_c, geometry=cl_geo,
                    individual={n: dict(best_kernel=KC[best_kc[i]], chi2_per_point=float(chi[i, best_kc[i]]/len(cls[n]['rp'])), lam=float(lamc[i, best_kc[i]])) for i, n in enumerate(names)},
                    F2=dict(kernel=KC[kc2], chi2_per_point=float(tot[kc2]), per_cluster_chi2_per_point={n: float(chi[i, kc2]/len(cls[n]['rp'])) for i, n in enumerate(names)},
                            lambda_by_cluster={n: float(lamc[i, kc2]) for i, n in enumerate(names)}, chi2_per_point_by_kernel=tot.tolist(),
                            described=bool(tot[kc2] <= 2*ref_release), regression=RW.regress_loglog(Mc, lamc[:, kc2]),
                            regression_h=RW.regress_loglog(Mc, lamc[:, kc2]*(Rhc/RW.R_STAR)**(2 - KC[kc2][0])),
                            newton_control=dict(kernel=KC[p2c[int(np.argmin(tot[p2c]))]], chi2_per_point=float(tot[p2c].min()))),
                    F3=dict(kernel=KC[kc3], lam=uni[kc3]['lam'], chi2_per_point=float(tot3[kc3]), chi2_per_point_by_kernel=tot3.tolist(), lambda_by_kernel=[u['lam'] for u in uni],
                            per_cluster_chi2_per_point={n: float(uni[kc3]['per_cluster_chi2'][i]/len(cls[n]['rp'])) for i, n in enumerate(names)}, described=bool(tot3[kc3] <= 2*ref_release),
                            newton_control=dict(kernel=KC[p2c[int(np.argmin(tot3[p2c]))]], lam=uni[p2c[int(np.argmin(tot3[p2c]))]]['lam'], chi2_per_point=float(tot3[p2c].min()))),
                    baryons_alone_per_point=refs['newtonian_baryons']['chi2_per_point'], hist_p={str(p): int(sum(1 for i in range(len(names)) if KC[best_kc[i]][0] == p)) for p in RW.P_GRID})
    out['C1_worst_kkt'] = kkt_worst
    gates['C1_nnls'] = kkt_worst < 1e-10
    print('clusters: F2 %s %.3f per point, F3 %s lam %.4g %.3f, reference %.3f' % (KC[kc2], tot[kc2], KC[kc3], uni[kc3]['lam'], tot3[kc3], ref_release), flush=True)
    # ---------------------------------------------------------------- the universal candidates
    gw = galaxies['rec_whole']
    U = {}
    U['U1'] = dict(kernel=tuple(gw['F3']['kernel']), lam=gw['F3']['lam'], law='universal strength')
    reg = gw['F2']['regression']
    U['U2'] = dict(kernel=tuple(gw['F2']['kernel']), slope=reg['slope'], intercept=reg['intercept'], law='lambda = 10^(intercept) M^slope from the G-whole line')
    regc = clusters['F2']['regression']
    U['U3'] = dict(kernel=tuple(clusters['F2']['kernel']), slope=regc['slope'], intercept=regc['intercept'], law='lambda = 10^(intercept) M^slope from the cluster line')

    def lam_of(u, M):
        return u['lam'] if 'lam' in u else (10**(u['intercept'] + u['slope']*np.log10(M)) if u['slope'] is not None else 0.)

    def gal_columns(g, split, kernel):
        k = kernel_index(KG, *kernel)
        if k is not None:
            return g['cols'][:, k]
        if 'extra' not in g:
            g['extra'] = {}
        if kernel not in g['extra']:
            rec = next(x for x in CS.I.sparc_galaxies() if x['name'] == g['name'])
            g['extra'][kernel] = RW.galaxy_whirl_columns(rec, g['R'], [kernel])[0][:, 0]
        return g['extra'][kernel]

    def cl_columns(name, kernel):
        k = kernel_index(KC, *kernel)
        if k is not None:
            return zc[name + '/ops'][:, k], zc[name + '/g'][:, k]
        ops, gs = RW.cluster_whirl_columns(cls[name], [kernel])
        return ops[:, 0], gs[:, 0]

    def score_galaxies(u, split, regime='whole'):
        gals = gals_by_split[split]
        vals = []
        for g in gals:
            m = g['masks'][regime]
            if m.sum() < 2:
                continue
            lam = lam_of(u, g['M'])
            uu = g['R'][m]*(g['gN_rec'][m] + lam*gal_columns(g, split, u['kernel'])[m])
            if np.any(uu <= 0):
                return None
            vals.append(np.mean((np.sqrt(uu) - g['v'][m])**2))
        return float(np.sqrt(np.mean(vals)))

    def score_clusters(u):
        tot, per = 0., {}
        for name in names:
            cl = cls[name]
            col, _ = cl_columns(name, u['kernel'])
            lam = lam_of(u, cl_geo[name]['M_b'])
            A_w, nu_w, y_w, _ = L2.whitened_cluster(cl, corrs[name], col[:, None])
            nu, res = nnls(nu_w, y_w - lam*A_w[:, 0])
            per[name] = float(res*res/len(cl['rp']))
            tot += res*res
        return dict(chi2_per_point=float(tot/N_c), per_cluster=per, described=bool(tot/N_c <= 2*ref_release))
    for key, u in U.items():
        u['galaxies'] = {split: score_galaxies(u, split) for split in ('train', 'validation', 'test')}
        u['galaxy_regimes'] = {r: score_galaxies(u, 'train', r) for r in REGIMES}
        u['clusters'] = score_clusters(u)
        print(key, u['kernel'], u['galaxies'], u['clusters']['chi2_per_point'], flush=True)
    # S3: the clusters at the galaxies' F2 kernel with their own strengths, against the galaxy line extrapolated to their masses
    kg_cl = kernel_index(KC, *gw['F2']['kernel'])
    s3 = {}
    for i, name in enumerate(names):
        if kg_cl is not None:
            f = dict(lam=float(lamc[i, kg_cl]), chi2=float(chi[i, kg_cl]))
        else:
            col, _ = cl_columns(name, tuple(gw['F2']['kernel']))
            A1, nu1, y1, _ = L2.whitened_cluster(cls[name], corrs[name], col[:, None])
            ff = RW.cluster_fit(A1[:, 0], nu1, y1)
            f = dict(lam=ff['lam'], chi2=ff['chi2'])
        line = float(lam_of(U['U2'], cl_geo[name]['M_b']))
        s3[name] = dict(lam_at_galaxy_kernel=f['lam'], chi2_per_point_at_galaxy_kernel=f['chi2']/len(cls[name]['rp']), galaxy_line_lambda=line,
                        ratio_at_galaxy_kernel=(f['lam']/line if line > 0 else None), lam_at_cluster_kernel=float(lamc[i, kc2]),
                        ratio_cluster_kernel_to_line=(float(lamc[i, kc2])/line if line > 0 else None))
    med = lambda key: (float(np.median([v[key] for v in s3.values() if v[key] is not None])) if any(v[key] is not None for v in s3.values()) else None)
    out['S3'] = dict(per_cluster=s3, galaxy_kernel=gw['F2']['kernel'], cluster_kernel=KC[kc2],
                     chi2_per_point_at_galaxy_kernel=float(sum(v['chi2_per_point_at_galaxy_kernel']*len(cls[n]['rp']) for n, v in s3.items())/N_c),
                     median_ratio_at_galaxy_kernel=med('ratio_at_galaxy_kernel'), median_ratio_cluster_kernel=med('ratio_cluster_kernel_to_line'))
    print('S3', out['S3']['chi2_per_point_at_galaxy_kernel'], out['S3']['median_ratio_at_galaxy_kernel'], out['S3']['median_ratio_cluster_kernel'], flush=True)
    # ---------------------------------------------------------------- the Milky Way
    zm = np.load(GEN/'rw1_milky_way.npz')
    mw = {}
    r3 = {}
    for variant in ('I', 'II'):
        R, y, gN, cols = zm[variant + '/R'], zm[variant + '/y'], zm[variant + '/gN'], zm[variant + '/cols']
        geo = zm[variant + '/geometry']
        r3[variant] = float(np.sqrt(np.mean((np.sqrt(R*gN) - y)**2)))
        mw[variant] = dict(baryons_alone=r3[variant], archived_nl1=archn['E3']['milky_way_baryons_only'][variant], R_half=float(geo[0]), M=float(geo[1]))
        for key, u in U.items():
            k = kernel_index(KG, *u['kernel'])
            col = cols[:, k] if k is not None else RW.milky_way_columns(variant, R, [u['kernel']])[0][:, 0]
            lam = lam_of(u, float(geo[1]))
            uu = R*(gN + lam*col)
            mw[variant][key] = float(np.sqrt(np.mean((np.sqrt(np.maximum(uu, 0)) - y)**2))) if np.all(uu > 0) else None
    out['R3'] = {v: abs(r3[v]/mw[v]['archived_nl1'] - 1) for v in ('I', 'II')}
    gates['R3_milky_way'] = max(out['R3'].values()) < 1e-9
    out['milky_way'] = mw
    print('Milky Way', mw, flush=True)
    # ---------------------------------------------------------------- the lenses
    lens_kernels = list(KG)
    for u in U.values():
        if kernel_index(lens_kernels, *u['kernel']) is None:
            lens_kernels.append(tuple(u['kernel']))
    lens_cache = GEN/'rw1_lenses.npz'
    key_k = hashlib.sha256(np.array(lens_kernels).tobytes()).hexdigest()[:16]
    cached = dict(np.load(lens_cache)) if lens_cache.exists() else {}
    lenses = []
    for n in CS.LENSES:
        kb, kd = f'{n}/{key_k}/basis', f'{n}/{key_k}/bend'
        if kb in cached:
            lenses.append(RW.LensWhirl(n, lens_kernels, cached=(cached[kb], cached[kd])))
        else:
            L = RW.LensWhirl(n, lens_kernels)
            cached[kb], cached[kd] = L.basis, L.basis_bend
            lenses.append(L)
    np.savez_compressed(lens_cache, **cached)
    print('lenses built at %.0f s' % (time.time() - t0), flush=True)
    # K9
    sub = [(p, 1.) for p in RW.P_GRID]
    Lf = RW.LensWhirl('J0037-0942', sub, full=True)
    Li = RW.LensWhirl('J0037-0942', sub)
    out['K9'] = float(np.max(np.max(np.abs(Li.basis - Lf.basis), axis=1)/np.max(np.abs(Lf.basis), axis=1)))
    gates['K9_lens_interp'] = out['K9'] < 1e-5
    # R5: stage 1's stars-only benchmarks under its geometry
    r5 = {}
    for n in CS.LENSES:
        so = CS.LensSystem(n, arch1['primary_geometry']).stars_only()
        a = arch1['E2_lenses']['cr2_benchmarks'][n]['stars_only']
        r5[n] = dict(chi2=so['chi2'], archived=a['chi2'], relative=abs(so['chi2']/a['chi2'] - 1), beta=so['beta'], archived_beta=a['beta'])
    out['R5'] = r5
    gates['R5_lenses'] = max(v['relative'] for v in r5.values()) < 1e-9
    lens_out = dict(per_lens={}, at_galaxy_F2_kernel={}, transfers={})
    kg2 = kernel_index(lens_kernels, *gw['F2']['kernel'])
    for L in lenses:
        name = L.name
        chis, lams, betas, eins = [], [], [], []
        for k in range(len(KG)):
            lam = L.einstein_lambda('Chabrier', k)
            amps = np.zeros(len(lens_kernels))
            amps[k] = lam
            chi2, beta = L.fit_beta('Chabrier', amps)
            ev = L.evaluate('Chabrier', beta, amps)
            chis.append(chi2); lams.append(lam); betas.append(beta); eins.append(ev['einstein_residual'])
        kb = int(np.argmin(chis))
        lens_out['per_lens'][name] = dict(best_kernel=KG[kb], chi2=float(chis[kb]), lam=float(lams[kb]), beta=float(betas[kb]), einstein_residual=float(eins[kb]),
                                          chi2_by_kernel=[float(c) for c in chis], lambda_by_kernel=[float(l) for l in lams], stellar_mass=float(L.L.pop['Chabrier']),
                                          stars_only_lambda_zero=bool(lams[kg2] == 0.))
        lens_out['at_galaxy_F2_kernel'][name] = dict(lam=float(lams[kg2]), chi2=float(chis[kg2]), beta=float(betas[kg2]), einstein_residual=float(eins[kg2]))
        print('lens %s: best %s chi2 %.1f lam %.3g; at galaxy kernel chi2 %.1f lam %.3g' % (name, KG[kb], chis[kb], lams[kb], chis[kg2], lams[kg2]), flush=True)
    tot_gal_kernel = float(sum(v['chi2'] for v in lens_out['at_galaxy_F2_kernel'].values()))
    lens_out['at_galaxy_F2_kernel_total_chi2'] = tot_gal_kernel
    lens_out['at_galaxy_F2_kernel_matched'] = bool(tot_gal_kernel <= L2.CHI2_LIMIT)
    lens_out['lambda_regression_at_galaxy_kernel'] = RW.regress_loglog([v['stellar_mass'] for v in lens_out['per_lens'].values()], [v['lam'] for v in lens_out['at_galaxy_F2_kernel'].values()], n_boot=200)
    for key, u in U.items():
        k = kernel_index(lens_kernels, *u['kernel'])
        per, tot, worst = [], 0., 0.
        for L in lenses:
            amps = np.zeros(len(lens_kernels))
            amps[k] = lam_of(u, L.L.pop['Chabrier'])
            chi2, beta = L.fit_beta('Chabrier', amps)
            ev = L.evaluate('Chabrier', beta, amps)
            per.append(dict(lens=L.name, chi2=float(chi2), beta=float(beta), einstein_residual=float(ev['einstein_residual']), lam=float(amps[k])))
            tot += chi2
            worst = max(worst, abs(ev['einstein_residual']))
        u['lenses'] = dict(per_lens=per, total_kinematics_chi2=float(tot), max_abs_einstein_residual=float(worst), described=bool(tot <= L2.CHI2_LIMIT and worst <= L2.EINSTEIN_LIMIT))
        print(key, 'lenses', tot, worst, flush=True)
    out['lenses'] = lens_out
    # ---------------------------------------------------------------- the cluster lensing curves
    curves = {}
    for i, name in enumerate(names):
        cl = cls[name]
        r = cl['r']
        h = ext['clusters'][name]['hydro_mass']
        M_nfw = CS._loglog_interp(r, h['RADIUS'], h['M_NFW'], left_power=2.)
        g_ref = G*M_nfw/r**2
        g_N = G*cl['M_b']/r**2
        a_ref = RW.deflection_curve(r, g_ref, 2., BS_KPC)
        a_N = RW.deflection_curve(r, g_N, 2., BS_KPC)
        r_ext = np.geomspace(r[0], 100*r[-1], 3000)
        cases = {'F2': (KC[kc2], lamc[i, kc2])}
        for key, u in U.items():
            cases[key] = (tuple(u['kernel']), lam_of(u, cl_geo[name]['M_b']))
        rec = dict(alpha_reference_rad=a_ref.tolist(), alpha_baryons_over_reference=(a_N/a_ref).tolist())
        for key, (kernel, lam) in cases.items():
            gw_ext = RW.shell_apply(r_ext, cl['source'].nodes, cl['source'].dM, kernel[0], kernel[1])
            a_w = RW.deflection_curve(r_ext, lam*gw_ext, kernel[0], BS_KPC)
            rec[key] = dict(kernel=kernel, lam=float(lam), alpha_over_reference=((a_N + a_w)/a_ref).tolist())
        curves[name] = rec
    out['cluster_lensing_curves'] = dict(impact_parameters_kpc=list(BS_KPC), per_cluster=curves,
                                         median_over_clusters={key: [float(np.median([curves[n][key]['alpha_over_reference'][j] for n in names])) for j in range(3)] for key in cases})
    print('cluster lensing curves done at %.0f s' % (time.time() - t0), flush=True)
    # ---------------------------------------------------------------- Coma
    kubo = json.loads(C1.KUBO.read_text(encoding='utf-8'))
    Rk = np.array([x['published_radius_h_inverse_Mpc'] for x in kubo['rows']])*1000/C1.H_COMA
    yk = np.array([x['shear_t'] for x in kubo['rows']])
    ek = np.array([x['plotted_sigma_t'] for x in kubo['rows']])

    def shape_chi2(ds):
        amp = max(float(np.sum(ds*yk/ek**2)/np.sum((ds/ek)**2)), 0.)
        return float(np.sum(((amp*ds - yk)/ek)**2)), amp
    coma = dict(brackets=[])
    cfref = json.loads((HERE/'cluster-lensing-forward/results.json').read_text(encoding='utf-8'))
    for j, (ne0, ms) in enumerate(C1.COMA_BRACKET):
        model = C1.coma_model(ne0, ms)
        r, rho = model['r'], model['rho']
        from scipy.integrate import cumulative_trapezoid
        mass = cumulative_trapezoid(4*np.pi*r*r*rho, r, initial=0) + 4*np.pi/3*r[0]**3*rho[0]
        Mp_b, Sig_b = CLZ.BaryonShells(r, mass).project(Rk)
        ds_b = Mp_b/(np.pi*Rk*Rk) - Sig_b
        chi_b, amp_b = shape_chi2(ds_b)
        archived = cfref['coma_transfer']['brackets'][j]['fits']['baryons']['shape_chi2']
        rec = dict(ne0=ne0, mstar=ms, M_b=model['M_b'], baryons_shape_chi2=chi_b, archived_baryons=archived, relative=abs(chi_b/archived - 1), candidates={})
        src = SF.SphericalSource(r, mass)
        for key, u in U.items():
            kernel = tuple(u['kernel'])
            lam = lam_of(u, model['M_b'])
            vals = {}
            for r_max in (1e5, 1e6):
                r_ext = np.geomspace(r[0], r_max, 4000)
                gwl = lam*RW.shell_apply(r_ext, src.nodes, src.dM, kernel[0], kernel[1])
                Rg, Sig_f, Mp_f = RW.project_effective(r_ext, gwl, r_max)
                ds_f = np.interp(Rk, Rg, Mp_f)/(np.pi*Rk*Rk) - np.interp(Rk, Rg, Sig_f)
                vals[str(int(r_max))] = shape_chi2(ds_b + ds_f)[0]
            rec['candidates'][key] = dict(kernel=kernel, lam=float(lam), shape_chi2_field_to_1e6_kpc=vals['1000000'], shape_chi2_field_to_1e5_kpc=vals['100000'])
        coma['brackets'].append(rec)
    out['R4'] = max(b['relative'] for b in coma['brackets'])
    gates['R4_coma'] = out['R4'] < 1e-6
    out['coma'] = coma
    print('Coma', [(b['baryons_shape_chi2'], {k: v['shape_chi2_field_to_1e6_kpc'] for k, v in b['candidates'].items()}) for b in coma['brackets']], flush=True)
    # ---------------------------------------------------------------- decision
    for key, u in U.items():
        u['passes'] = dict(galaxies=bool(u['galaxies']['train'] is not None and u['galaxies']['train'] <= V_ACCEPT), clusters=u['clusters']['described'], lenses=u['lenses']['described'])
        u['passes_all'] = all(u['passes'].values())
    regimes_matched = {r: galaxies[f'rec_{r}']['F2']['matched'] for r in REGIMES}
    regimes_matched['clusters'] = clusters['F2']['described']
    regimes_matched['lenses'] = lens_out['at_galaxy_F2_kernel_matched']
    if any(u['passes_all'] for u in U.values()):
        decision = 'a'
        detail = 'a: ' + ', '.join(k for k, u in U.items() if u['passes_all']) + ' passes the galaxies, the clusters and the lenses at once'
    elif all(regimes_matched.values()):
        decision = 'b'
        detail = 'b: every regime is matched by its own formula, no universal candidate passes all three; ' + '; '.join(
            '%s fails %s' % (k, ', '.join(n for n, p in u['passes'].items() if not p)) for k, u in U.items())
    else:
        decision = 'c'
        detail = 'c: not matched even by its own formula: ' + ', '.join(r for r, m in regimes_matched.items() if not m)
    out['universal_candidates'] = U
    out['clusters'] = clusters
    out['regimes_matched'] = regimes_matched
    out['gates'] = gates
    out['decision'] = decision
    out['numerical_verification_passed'] = bool(all(gates.values()))
    out['statuses'] = dict(numerical_verification='passed' if all(gates.values()) else 'FAILED on: ' + ', '.join(k for k, v in gates.items() if not v), scientific_outcome=detail)
    out['checks_short_run'] = dict(galaxy_F2_whole_kernel=gw['F2']['kernel'], galaxy_F2_whole_rmse=gw['F2']['rmse'], galaxy_F3_whole=[gw['F3']['kernel'], gw['F3']['lam'], gw['F3']['rmse']],
                                   lambda_NGC2403=gw['F2']['lambda_by_galaxy'][[g['name'] for g in train].index('NGC2403')],
                                   lambda_DDO064=gw['F2']['lambda_by_galaxy'][[g['name'] for g in train].index('DDO064')],
                                   lambda_NGC2841=gw['F2']['lambda_by_galaxy'][[g['name'] for g in train].index('NGC2841')],
                                   cluster_F2=[clusters['F2']['kernel'], clusters['F2']['chi2_per_point']],
                                   rmse_F2_whole_by_name={n: float(np.sqrt(F1[('rec', 'whole')][0][[g['name'] for g in train].index(n), kernel_index(KG, *gw['F2']['kernel'])])) for n in ('NGC2403', 'DDO064', 'NGC2841')})
    out['runtime_seconds'] = time.time() - t0
    OUT.write_text(json.dumps(plain(out), indent=1), encoding='utf-8')
    print('written', OUT, 'gates', gates, 'decision', detail, 'runtime %.0f s' % out['runtime_seconds'], flush=True)


if __name__ == '__main__':
    main()

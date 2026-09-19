"""NL-1 driver (protocol-nl1.md, amendment 1): the root spectrum, linear in the roots of its amplitudes, on the galaxies,
clusters and lenses, jointly and in transfer. Writes nl1-results.json."""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import cl2s2 as D              # noqa: E402
import cl2s2_lib as L2         # noqa: E402
import routes_lib as RL        # noqa: E402
import cl2_response as CR      # noqa: E402
import cl2_sources as CS       # noqa: E402

GEN = ROOT/'research_work/generated/routes'
OUT = HERE/'nl1-results.json'
W = RL.WIDTHS
V_ACCEPT, CLUSTER_REF = RL.V_ACCEPT, RL.CLUSTER_REF


def galaxy_root_blocks(blocks):
    info = {}
    for split, b in blocks.items():
        Ar, inf = RL.root_basis(b['A'])
        b['A_root'] = {key: np.hstack([np.sqrt(np.maximum(b[key], 0))[:, None], Ar]) for key in ('gN_rec', 'gN_tab')}
        info[split] = inf
    return info


def certified(A, gN, block, weights=None):
    loss = L2.SpeedLoss(A, gN, block['R'], block['v'], block['sizes'])
    if weights is not None:
        loss.w = weights
    s = np.maximum(2*block['v']*block['ev']/block['R'], 1e-3*block['v']**2/block['R'])
    L0 = L2.nnls_start(A, block['v']**2/block['R'] - gN, s, np.ones(len(s)))
    if not np.isfinite(loss.value(L0)):
        L0 = np.zeros(A.shape[1])
    La, ra = loss.solve_lbfgs(L0)
    Lb, rb = loss.solve_trust_newton(fallback=L0)
    Fa, Fb = loss.value(La), loss.value(Lb)
    best = Lb if Fb <= Fa else La
    return loss, best, L0, dict(F_lbfgs=Fa, F_newton=Fb, solver_agreement=abs(Fa - Fb)/max(abs(Fb), 1e-300), kkt=loss.kkt(best, L0),
                                curvature_min_optimum=loss.curvature_check(best, 200), newton_start=rb['start'], newton_message=rb['message'])


def multistart(loss, best, L0, a_local, n=20, seed=11):
    """The declared twenty starts: the local law alone, the local law plus a scaled spectrum, eighteen random; each must
    reach the certified optimum within 1e-6 (amendment 1)."""
    rng = np.random.default_rng(seed)
    F = loss.value(best)
    starts = [np.r_[np.sqrt(a_local), np.zeros(len(best) - 1)], np.r_[np.sqrt(a_local), L0[1:]*.5]]
    starts += [np.abs(rng.normal(size=len(best)))*np.max(best)*rng.uniform(.1, 2.) for _ in range(n - 2)]
    finals = []
    for x0 in starts:
        if not np.isfinite(loss.value(x0)):
            x0 = np.maximum(x0*.1, 0)
            if not np.isfinite(loss.value(x0)):
                x0 = np.zeros_like(x0)
        x, _ = loss.solve_lbfgs(x0)
        x, _ = loss.solve_trust_newton(L0=x, fallback=x)
        finals.append(loss.value(x))
    finals = np.array(finals)
    return dict(best=F, finals=finals.tolist(), within_1e6=int(np.sum(np.abs(finals/F - 1) < 1e-6)), worst_relative=float(np.max(np.abs(finals/F - 1))))


def main():
    t0 = time.time()
    ext = json.loads((HERE/'cl2-inputs-xcop-profiles.json').read_text(encoding='utf-8'))
    frac = CS.stellar_fraction_profile(ext)
    arch2 = json.loads((HERE/'cl2s2-results.json').read_text(encoding='utf-8'))
    archive1 = json.loads((HERE/'cl2-results.json').read_text(encoding='utf-8'))
    blocks, _ = D.load_galaxies()
    neg = galaxy_root_blocks(blocks)
    tr = blocks['train']
    out = dict(experiment='NL-1', protocol='protocol-nl1.md; amendments 1 and 2', widths_kpc=W.tolist(), members=['local'] + W.tolist(),
               negative_entries=neg, input_sha256={'cl2s2-results.json': hashlib.sha256((HERE/'cl2s2-results.json').read_bytes()).hexdigest()})
    gates = {}
    # ---- R1: the local law at PM-1's fitted a* on the tabulated force reproduces the archived reference
    loss_tab = L2.SpeedLoss(tr['A_root']['gN_tab'], tr['gN_tab'], tr['R'], tr['v'], tr['sizes'])
    s_pm1 = np.r_[np.sqrt(CS.A_STAR_PM1), np.zeros(len(W))]
    r1 = loss_tab.rmse(s_pm1)
    ref = archive1['E1_spectrum']['references']['galaxies_train']['pm1_law_fitted_a_star']['equal_galaxy_rmse_km_s']
    gates['R1_local_law'] = abs(r1/ref - 1) < 1e-3
    out['R1'] = dict(local_law_tabulated=r1, archived=ref, relative=abs(r1/ref - 1))
    # ---- galaxies alone (certified; the declared multi-start must agree)
    gal = {}
    for label, key in (('reconstructed', 'gN_rec'), ('tabulated', 'gN_tab')):
        loss, best, L0, cert = certified(tr['A_root'][key], tr[key], tr)
        rng = np.random.default_rng(3)
        cert['curvature_min_random'] = min(loss.curvature_check(np.maximum(best*(1 + .5*rng.normal(size=len(best))) + 1e-3*np.max(best)*np.abs(rng.normal(size=len(best))), 0), 10, seed=k) for k in range(20))
        ms = multistart(loss, best, L0, CS.A_STAR_PM1)
        # the local-law-only member, a* free (one-dimensional), for the reading
        from scipy.optimize import minimize_scalar
        f_local = lambda a: loss.value(np.r_[np.sqrt(max(a, 0)), np.zeros(len(W))])
        o = minimize_scalar(f_local, bounds=(1., 1e5), method='bounded')
        r = dict(**cert, multistart=ms, rmse_train=loss.rmse(best), amplitudes_s=best.tolist(), amplitudes_a=(best**2).tolist(),
                 active_members=[out['members'][i] for i in np.flatnonzero(best > 0)],
                 local_only=dict(a_star=float(o.x), rmse=float(np.sqrt(o.fun))),
                 per_galaxy_rmse=dict(zip(tr['names'], loss.per_galaxy_rmse(best))), slope=D.slope_diagnostic(dict(tr, A=tr['A_root'][key]), best, tr[key]))
        for split in ('validation', 'test'):
            b = blocks[split]
            v = L2.SpeedLoss(b['A_root'][key], b[key], b['R'], b['v'], b['sizes']).value(best)
            r[f'rmse_{split}'] = float(np.sqrt(v)) if np.isfinite(v) else None
        r['described'] = bool(r['rmse_train'] <= V_ACCEPT)
        gal[label] = r
        if label == 'reconstructed':
            loss_rec, best_rec = loss, best
        print(label, 'root spectrum certified minimum %.3f km/s (local law alone %.3f at a* %.0f); multistart within 1e-6: %d of 20' % (r['rmse_train'], r['local_only']['rmse'], r['local_only']['a_star'], ms['within_1e6']), flush=True)
    gates['G1_kkt'] = all(gal[k]['kkt']['max_abs_gradient_active'] < 1e-8 and gal[k]['kkt']['min_gradient_inactive'] > -1e-8 for k in gal)
    gates['G2_multistart'] = all(gal[k]['multistart']['within_1e6'] == 20 for k in gal)
    gates['G3_gradient'] = True   # set below from the finite-difference check
    # G3: analytic gradient of the loss on the root basis against finite differences at 20 random feasible points
    rng = np.random.default_rng(5)
    worst = 0.
    for _ in range(20):
        x = np.maximum(best_rec*(1 + .3*rng.normal(size=len(best_rec))), 0) + 1e-4*np.max(best_rec)
        g = loss_rec.gradient(x)
        for j in rng.choice(len(x), 5, replace=False):
            e = np.zeros_like(x)
            e[j] = 1e-6*max(1., abs(x[j]))
            fd = (loss_rec.value(x + e) - loss_rec.value(x - e))/(2*e[j])
            worst = max(worst, abs(fd - g[j])/max(abs(g[j]), 1e-12))
    gates['G3_gradient'] = worst < 1e-6
    out['G3_worst'] = worst
    gates['C1_solvers'] = all(gal[k]['solver_agreement'] < 1e-6 for k in gal)
    gates['C3_convex'] = all(gal[k]['curvature_min_optimum'] >= 0 and gal[k]['curvature_min_random'] >= 0 for k in gal)
    gates['G5_negative_entries'] = neg['train']['negative_entries']/neg['train']['entries'] < 1e-3
    out['galaxies'] = gal
    # ---- clusters: the root columns on the 1500-point grid, the SZ correlations, the corrected boundary fit
    # amendment 2: the root columns on stage 2's 6000-point grid, every width included; the same build's linear columns
    # must reproduce stage 2's archived cache exactly
    z = np.load(GEN/'nl1_clusters6000.npz')
    corrs = L2.sz_correlations(D.TAR, list(ext['clusters']))
    cls = {}
    g4 = 0.
    z6 = np.load(D.GEN/'cluster_ops.npz')
    for name in ext['clusters']:
        cl = CS.build_cluster(name, ext['clusters'][name], np.array([100.]), n_grid=6000, frac_profile=frac)
        cl['ops'] = z[name + '/root']
        assert np.allclose(z[name + '/rp'], cl['rp'])
        lin, ref6 = z[name + '/linear6000'], z6[f'{name}_6000']
        g4 = max(g4, float(np.max(np.abs(lin - ref6))/np.max(np.abs(ref6))))
        cls[name] = cl
    gates['G4b_cluster_identity'] = g4 < 1e-12
    out['G4b_worst'] = g4
    out['cluster_negative_samples'] = json.loads((GEN/'nl1_cluster_negatives6000.json').read_text())
    block, Lcs = D.cluster_blocks(cls, corrs)
    csol = D.cluster_solve(block)
    refs = D.cluster_references(cls, Lcs, corrected=True)
    # R1 for the clusters: the local member at PM-1's a* reproduces stage 2's PM-1 references under the correlations
    pm = {}
    for label, a_star in (('pm1_law_fitted_a_star', CS.A_STAR_PM1), ('pm1_law_repo_a_star', CS.A_STAR_REPO)):
        s0 = np.r_[np.sqrt(a_star), np.zeros(len(W))]
        nu, _ = nnls(block.nuisance, block.y - block.A@s0)
        chi = float(np.sum((block.A@s0 + block.nuisance@nu - block.y)**2))/block.N
        pm[label] = dict(chi2_per_point=chi, stage2=arch2['clusters']['references_covariance'][label]['chi2_per_point'], relative=abs(chi/arch2['clusters']['references_covariance'][label]['chi2_per_point'] - 1))
    gates['R1_cluster_local_law'] = all(v['relative'] < 1e-9 for v in pm.values())
    out['clusters'] = dict(alone=dict(chi2_per_point=csol['chi2_per_point'], amplitudes_s=csol['amplitudes'], active_members=[out['members'][i] for i in csol['widths_active']],
                                      kkt=csol['kkt'], per_cluster_chi2_per_point=csol['per_cluster_chi2_per_point'], described=bool(csol['chi2_per_point'] <= 2*CLUSTER_REF)),
                           references=refs, pm1_reproduction=pm, reference_per_point=refs['release_nfw_reference']['chi2_per_point'])
    print('clusters alone: root spectrum %.3f per point (reference %.3f, PM-1 local %.1f)' % (csol['chi2_per_point'], refs['release_nfw_reference']['chi2_per_point'], pm['pm1_law_fitted_a_star']['chi2_per_point']), flush=True)
    # ---- lenses
    lenses = [RL.LensRoot(n) for n in CS.LENSES]
    print('lenses built at %.0f s' % (time.time() - t0), flush=True)
    betas0 = [Lz.L.stars_only()['beta'] for Lz in lenses]
    betas = list(betas0)
    hist = []
    for it in range(6):
        E, K = D.lens_blocks(lenses, 'Chabrier', betas)
        sol = CR.solve([E, K])
        ev = D.lens_evaluate(lenses, 'Chabrier', sol['amplitudes'])
        hist.append(dict(iteration=it, objective=float(sum(c/b.N for c, b in zip(sol['chi2'], (E, K)))), total=ev['total_kinematics_chi2'], einstein=ev['max_abs_einstein_residual']))
        new = [e['beta'] for e in ev['per_lens']]
        if it and abs(hist[-1]['objective'] - hist[-2]['objective']) < 1e-6*max(1., hist[-2]['objective']):
            break
        betas = new
    out['lenses'] = dict(alone=dict(total_kinematics_chi2=ev['total_kinematics_chi2'], max_abs_einstein_residual=ev['max_abs_einstein_residual'], described=ev['described'],
                                    amplitudes_s=sol['amplitudes'].tolist(), active_members=[out['members'][i] for i in sol['widths_active']], history=hist, betas=betas, per_lens=ev['per_lens']))
    print('lenses alone', ev['total_kinematics_chi2'], ev['max_abs_einstein_residual'], flush=True)
    # ---- the balanced joint solve (convex): F/21.9^2 + chi2_c/(N_c 2 ref) + chi2_K/128 + chi2_E/6
    loss_b = L2.SpeedLoss(tr['A_root']['gN_rec'], tr['gN_rec'], tr['R'], tr['v'], tr['sizes'])
    loss_b.w = loss_b.w/V_ACCEPT**2
    betas = list(betas0)
    s = np.array(best_rec)
    nu_c = np.array(csol['nuisances'])
    jh = []
    for it in range(6):
        E, K = D.lens_blocks(lenses, 'Chabrier', betas)
        js = D.joint_solve(loss_b, [block, E, K], [1./(block.N*2*CLUSTER_REF), 1./RL.KIN_LIMIT, 1./RL.N_LENSES], s, [nu_c, None, None])
        s, nu_c = js['L'], js['nus'][0]
        ev = D.lens_evaluate(lenses, 'Chabrier', s)
        terms = dict(galaxies=loss_rec.value(s)/V_ACCEPT**2, clusters=block.chi2(s, nu_c)/block.N/(2*CLUSTER_REF), kinematics=ev['total_kinematics_chi2']/RL.KIN_LIMIT,
                     einstein=float(sum((e['einstein_residual']/L2.EINSTEIN_LIMIT)**2 for e in ev['per_lens']))/RL.N_LENSES)
        jh.append(dict(iteration=it, objective=js['objective'], projected_gradient=js['projected_gradient'], terms=terms))
        new = [e['beta'] for e in ev['per_lens']]
        if it and abs(jh[-1]['objective'] - jh[-2]['objective']) < 1e-6*max(1., jh[-2]['objective']):
            break
        betas = new
    gates['GJ_joint_certificate'] = js['projected_gradient'] < 1e-6*max(1., js['objective'])
    joint = dict(objective=js['objective'], projected_gradient=js['projected_gradient'], terms=terms, history=jh, amplitudes_s=s.tolist(),
                 active_members=[out['members'][i] for i in np.flatnonzero(s > 0)], betas=betas,
                 galaxies_rmse=loss_rec.rmse(s), clusters_chi2_per_point=block.chi2(s, nu_c)/block.N, lenses=dict(total_kinematics_chi2=ev['total_kinematics_chi2'], max_abs_einstein_residual=ev['max_abs_einstein_residual']),
                 described=dict(galaxies=bool(terms['galaxies'] <= 1), clusters=bool(terms['clusters'] <= 1), lenses=bool(ev['described'])))
    out['joint'] = joint
    print('joint terms', terms, flush=True)
    # ---- transfer
    mw = CS.milky_way_blocks(W)

    def mw_scores(sv):
        o = {}
        for b in mw:
            Ar, _ = RL.root_basis(b.A)
            basis = np.hstack([np.sqrt(np.maximum(b.meta['gN'], 0))[:, None], Ar])
            u = b.meta['R']*(b.meta['gN'] + basis@sv)
            o[b.meta['variant']] = float(np.sqrt(np.mean((np.sqrt(np.maximum(u, 0)) - b.meta['y'])**2)))
        return o

    def gal_scores(sv):
        return {split: (lambda v: float(np.sqrt(v)) if np.isfinite(v) else None)(L2.SpeedLoss(blocks[split]['A_root']['gN_rec'], blocks[split]['gN_rec'], blocks[split]['R'], blocks[split]['v'], blocks[split]['sizes']).value(sv)) for split in ('train', 'validation', 'test')}

    def cluster_score(sv):
        nu, _ = nnls(block.nuisance, block.y - block.A@sv)
        return float(np.sum((block.A@sv + block.nuisance@nu - block.y)**2)/block.N)
    sg, sc = np.array(best_rec), np.array(csol['amplitudes'])
    E3 = dict(calibrated_on_galaxies=dict(clusters_chi2_per_point=cluster_score(sg), lenses={k: v for k, v in D.lens_evaluate(lenses, 'Chabrier', sg).items() if k != 'per_lens'}, milky_way=mw_scores(sg), galaxies=gal_scores(sg)),
              calibrated_on_clusters=dict(galaxies=gal_scores(sc), lenses={k: v for k, v in D.lens_evaluate(lenses, 'Chabrier', sc).items() if k != 'per_lens'}, milky_way=mw_scores(sc)),
              joint=dict(galaxies=gal_scores(s), clusters_chi2_per_point=cluster_score(s), milky_way=mw_scores(s)),
              milky_way_baryons_only=mw_scores(np.zeros(len(W) + 1)))
    out['E3'] = E3
    # ---- decision
    rec = gal['reconstructed']
    if all(joint['described'].values()):
        outcome = 'a: the joint minimiser describes galaxies, clusters and lenses at once'
    elif rec['described']:
        outcome = 'b: the galaxies are described (%.2f km/s) but no joint minimiser brings the clusters (term %.2f) and lenses (kinematics term %.2f, Einstein term %.2f) inside' % (rec['rmse_train'], terms['clusters'], terms['kinematics'], terms['einstein'])
    else:
        outcome = 'c: the root spectrum does not describe the galaxies with the consistent source (%.2f km/s; tabulated %.2f)' % (rec['rmse_train'], gal['tabulated']['rmse_train'])
    out.update(gates=gates, numerical_verification_passed=bool(all(gates.values())),
               statuses=dict(numerical_verification='passed' if all(gates.values()) else 'FAILED on: ' + ', '.join(k for k, v in gates.items() if not v), scientific_outcome=outcome),
               checks_short_run=dict(rmse_reconstructed=rec['rmse_train'], F_reconstructed=rec['F_newton'], local_law_tabulated=r1, clusters_alone=csol['chi2_per_point']),
               runtime_seconds=time.time() - t0)
    OUT.write_text(json.dumps(D.plain(out), indent=1) + '\n', encoding='utf-8')
    print('written', OUT, 'runtime %.0f s' % out['runtime_seconds'], 'gates', gates, 'outcome:', outcome, flush=True)


if __name__ == '__main__':
    main()

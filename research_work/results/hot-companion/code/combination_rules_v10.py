"""Round 10, part 2: how the companion adds up. Three rules, tested on the galaxies and the clusters.

A. The local form of the pull (exact, for any field of one frequency). An emitter locked a quarter
   cycle ahead of the field it sits in, psi = Im[Psi(x) e^{i w t}], feels on average
       <F> = (q0/2) Im(Psi* grad Psi) / |Psi| = (q0/2) |Psi| grad(theta),     Psi = |Psi| e^{i theta}:
   the local amplitude times how fast the local phase turns. For one travelling wave |grad theta| = k
   and the pull is (q0 k/2)|Psi| (part 1, first_principles_v10.py). Where companion arrives from
   several directions at once, the phase turns more slowly on average and the pull is diluted: its
   average is the net current divided by the amplitude. Checked against the exact time average for
   random superpositions, and worked out for two coherent counter-propagating waves.
B. Three ways the companion of many pieces of matter could add up, written in the law's terms
   (S_N = G int rho/d^2 is the plain total of every piece's Newtonian pull; S and g_hot are the
   heat term's plain total and net flow; every rule points along the net flow g_N + g_hot):
     1. independent waves (part 1):       extra = sqrt(a) |g_N + g_hot| / sqrt(S_N + S)
     2. one stream carrying all of it:     extra = sqrt(a (S_N + S))
     3. round 3, the adopted law:          extra = sqrt(a (|g_N| + S))
   Rule 2 is what a companion that always moves at speed u, all in one direction at each point,
   would give if the ordered matter's opposing flows did not cancel.
C. SPARC, the 149 rotation curves, with S_N computed from each galaxy's surface densities (stars from
   its 3.6-micron profile at M/L 0.5; gas by inverting V_gas into thick annuli; the bulge spherical),
   and a and g_d refitted for each rule. Two disk thicknesses.
D. X-COP, the 12 clusters, u refitted with a and g_d held (as in part 1).
E. The Milky Way: a local estimate of the circular speed under each rule, with its own fitted a and g_d.

    python code/combination_rules_v10.py --output-dir run-combination-rules-v10
"""
from __future__ import annotations
import argparse, io, json, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar, nnls
from scipy.special import ellipk

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
import law as L                                   # noqa: E402
import run as RUN                                 # noqa: E402
ROOT = RUN.ROOT                                   # the repository root (data files)

G = L.G


# --------------------------------------------------------------------------- A. the local pull
def local_pull(src_pos, Qs, probe, k, alphas):
    """(q0 = 1) (1/2) Im(Psi* grad Psi)/|Psi| at the probe, for point sources with phases alphas [nreal, nsrc]."""
    d = probe[None, :] - src_pos; r = np.linalg.norm(d, axis=1); rhat = d / r[:, None]
    Z = (Qs / (4 * np.pi * r))[None, :] * np.exp(1j * (alphas - k * r[None, :]))
    Psi = Z.sum(1)
    dPsi = (Z * (-1 / r - 1j * k)[None, :]) @ rhat                       # [nreal, 3]
    return 0.5 * np.imag(np.conj(Psi)[:, None] * dPsi) / np.abs(Psi)[:, None]


def part_a(rng):
    import first_principles_v10 as FP
    k = 1.0; rows = []
    cases = dict(one_source=(np.array([[0.0, 0, 0]]), np.array([1.0]), np.array([60.0, 0, 0])),
                 two_same_side=(np.array([[0.0, 0, 0], [10.0, 0, 0]]), np.array([1.0, 0.7]), np.array([80.0, 0, 0])),
                 two_opposite=(np.array([[-60.0, 0, 0], [60.0, 0, 0]]), np.array([1.0, 0.6]), np.zeros(3)))
    p = rng.normal(0, 40.0, (400, 3)); cases['cloud_of_400'] = (p, np.ones(400) / 20.0, np.array([30.0, 0, 0]))
    for name, (pos, Qs, probe) in cases.items():
        al = rng.uniform(0, 2 * np.pi, (500, len(Qs)))
        exact = FP.locked_force_many(pos, Qs, probe, k, np.pi / 2, al)
        form = local_pull(pos, Qs, probe, k, al)
        rows.append(dict(config=name, max_abs_difference=float(np.abs(exact - form).max()),
                         typical_force=float(np.abs(exact).mean())))
    # two coherent waves travelling in opposite directions, amplitudes 1 and b: the pull at x is
    # (q0 k/2)(1 - b^2)/|Psi(x)|; its average over a wavelength is (q0 k/2)(1 - b)(2/pi) K(4b/(1+b)^2)
    two = []
    for b in (0.0, 0.25, 0.5, 0.75, 0.9):
        x = np.linspace(0, np.pi, 20001)[:-1]
        mod = np.sqrt(1 + b * b + 2 * b * np.cos(2 * x))
        avg = np.mean((1 - b * b) / mod)
        closed = (1 - b) * (2 / np.pi) * ellipk(4 * b / (1 + b) ** 2)
        two.append(dict(counter_amplitude=b, mean_pull_over_single_wave=float(avg), closed_form=float(closed),
                        net_flux_over_single_wave=float(1 - b * b), round3_pull_over_single_wave=float(np.sqrt(1 - b * b)),
                        energy_over_single_wave=float(1 + b * b)))
    return dict(identity_checks=rows, counter_propagating=two)


# --------------------------------------------------------------------------- C. SPARC plain totals
def zgrid(h, n=64):
    """nodes (bin midpoints, never z = 0) and weights for p(z) = exp(-|z|/h)/(2h), folded onto z > 0."""
    xe = np.sinh(np.linspace(0, np.arcsinh(14.0 / 1e-5), n + 1)) * 1e-5
    x = 0.5 * (xe[1:] + xe[:-1]); w = np.exp(-xe[:-1]) - np.exp(-xe[1:])
    return x * h, w / w.sum()


def ring_phi(R, a, z):
    """potential per unit mass of a thin ring of radius a (G = 1), at radius R and height z."""
    s = (R + a) ** 2 + z ** 2
    return -(2 / np.pi) * ellipk(np.clip(4 * R * a / s, 0, 1 - 1e-15)) / np.sqrt(s)


def agrid(Rf, rmax, h):
    g1 = np.geomspace(1e-3, 3 * rmax, 500)
    d = np.sinh(np.linspace(-np.arcsinh(8 / 1e-4), np.arcsinh(8 / 1e-4), 800)) * 1e-4 * h   # never exactly Rf
    g2 = Rf + d
    return np.unique(np.r_[g1, g2[g2 > 1e-4]])


def kernels(Rf, h, rmax):
    """on a grid of source radii a: S_N/G and g_R/G per unit surface density, for a disk of exponential
    thickness h, at midplane radius Rf."""
    a = agrid(Rf, rmax, h); z, wz = zgrid(h)
    A, Z = np.meshgrid(a, z, indexing='ij')
    kS = (2 * np.pi * A / np.sqrt(((Rf - A) ** 2 + Z ** 2) * ((Rf + A) ** 2 + Z ** 2))) @ wz
    dR = 1e-5 * Rf
    kg = 2 * np.pi * a * (((ring_phi(Rf + dR, A, Z) - ring_phi(Rf - dR, A, Z)) / (2 * dR)) @ wz)
    return a, kS, kg


def cumint(a, f):
    return np.r_[0.0, np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(a))]


def disk_sigma(r, SB, Rd):
    """stellar surface density (Msun/kpc^2) at M/L 0.5 from the 3.6-micron profile, exponential beyond the data."""
    S = 0.5 * np.maximum(SB, 1e-6) * 1e6; ls = np.log(S)
    def f(a):
        a = np.asarray(a, float); out = np.exp(np.interp(a, r, ls))
        out = np.where(a > r[-1], S[-1] * np.exp(-(a - r[-1]) / Rd), out)
        return np.where(a < r[0], S[0], out) * (SB.max() > 0)
    return f


def gas_sigma(r, vg2, K):
    """V_gas^2 (with its sign) inverted into piecewise-constant annuli: non-negative least squares with a
    light smoothness penalty. K[i] = kernels at r[i]."""
    e = np.r_[0.0, 0.5 * (r[1:] + r[:-1]), r[-1] + 0.5 * (r[-1] - r[-2])]
    e = np.r_[e, e[-1] * np.array([1.15, 1.35, 1.6])]
    nb = e.size - 1
    B = np.array([np.diff(np.interp(e, a, cumint(a, kg))) * r[i] * G for i, (a, kS, kg) in enumerate(K)])
    lam = 0.03 * np.abs(B).max()
    D = np.zeros((nb - 1, nb)); D[np.arange(nb - 1), np.arange(nb - 1)] = -1; D[np.arange(nb - 1), np.arange(1, nb)] = 1
    w, _ = nnls(np.vstack([B, lam * D]), np.r_[vg2, np.zeros(nb - 1)], maxiter=5000)
    def f(a):
        a = np.asarray(a, float); idx = np.searchsorted(e, a, side='right') - 1
        return np.where((idx >= 0) & (idx < nb), w[np.clip(idx, 0, nb - 1)], 0.0)
    return f


def sparc_plain_totals(gals, h_gas, zfac):
    """S_N at every measured radius of every galaxy; also the model in-plane pulls, to check the
    reconstruction against SPARC's own V_disk and V_gas."""
    zf = zipfile.ZipFile(ROOT / 'temporal_candidate_audit/data/Rotmod_LTG.zip')
    chk_d, chk_g = [], []
    for g in gals:
        a = np.atleast_2d(np.loadtxt(io.BytesIO(zf.read(g['name'] + '_rotmod.dat')))); a = a[a[:, 0] > 0]
        r, vg, vd, vb, SBd = a[:, 0], a[:, 3], a[:, 4], a[:, 5], a[:, 6]
        ok = (vg * np.abs(vg) + 0.5 * vd ** 2 + 0.7 * vb ** 2) / r > 0
        assert ok.sum() == g['r'].size and np.allclose(r[ok], g['r'])
        Rd = g['Rd']; hz = zfac * 0.196 * Rd ** 0.633
        rmax = max(r[-1] * 1.6, 5 * Rd)
        Kg = [kernels(Rf, h_gas, rmax) for Rf in r]
        fg = gas_sigma(r, vg * np.abs(vg), Kg); fd = disk_sigma(r, SBd, Rd)
        SN, gmod = [], []
        for i in np.where(ok)[0]:
            ad, kSd, kgd = kernels(r[i], hz, rmax); sd = fd(ad)
            ag, kSg, kgg = Kg[i]; sg = fg(ag)
            SN.append(np.trapezoid(sd * kSd, ad) + np.trapezoid(sg * kSg, ag))
            gd_, gg_ = np.trapezoid(sd * kgd, ad) * G, np.trapezoid(sg * kgg, ag) * G
            gmod.append(gd_ + gg_)
            sdk, sgk = 0.5 * vd[i] ** 2 / r[i], vg[i] * abs(vg[i]) / r[i]
            if sdk > 0.05 * np.max(0.5 * vd ** 2 / r): chk_d.append(gd_ / sdk)
            if abs(sgk) > 0.2 * np.max(np.abs(vg * np.abs(vg) / r)): chk_g.append(gg_ / sgk)
        SN = G * np.array(SN)
        if g['sigb'] > 0:                        # the bulge's own plain total (spherical shells)
            SN = SN + L.scalar_sum(g['r'], g['sfine'], g['dmb'], np.ones_like(g['dmb']))
        gmod = np.array(gmod) + g['vb2'] / g['r']
        # keep S_N / g_N from the reconstruction, applied to SPARC's own g_N (and never below it)
        g['SN'] = np.maximum(SN * g['gN'] / np.maximum(gmod, 1e-12 * g['gN']), g['gN'])
    pct = lambda x: [float(v) for v in np.percentile(x, [16, 50, 84])]
    ratio = np.concatenate([g['SN'] / g['gN'] for g in gals])
    return dict(h_gas_kpc=h_gas, stellar_thickness=f'{zfac:g} x 0.196 R_d^0.633 kpc',
                disk_pull_model_over_sparc_16_50_84=pct(chk_d), gas_pull_model_over_sparc_16_50_84=pct(chk_g),
                SN_over_gN_16_50_84=pct(ratio))


def heat(g, u):
    k = L.heat_weight(g['sigb'], u) if g['sigb'] > 0 else 0.0
    S = L.scalar_sum(g['r'], g['sfine'], g['dmb'], k * np.ones_like(g['dmb'])) if g['sigb'] > 0 else 0.0
    return S, k * g['vb2'] / g['r']


def galaxy_rule(rule, g, a, u, lam):
    S, ghot = heat(g, u); rel = L.released(g['gN'], a, lam)
    if rule == 'independent':
        ex = np.sqrt(a) * (g['gN'] + ghot) / np.sqrt(g['SN'] + S)
    elif rule == 'one_stream':
        ex = np.sqrt(a * (g['SN'] + S))
    else:
        ex = np.sqrt(a * (g['gN'] + S))
    return g['gN'] + rel * ex


def sparc_fits(gals, law):
    a0, lam0, u = law['a_code'], law['lam'], law['u_kms']
    out = {}
    for rule in ('round3', 'independent', 'one_stream'):
        fn = lambda a, lam: (lambda g: galaxy_rule(rule, g, a, u, lam))
        at0 = RUN.sparc_score(gals, fn(a0, lam0))
        cost = lambda p: RUN.sparc_score(gals, fn(10 ** p[0], 10 ** p[1]))[1]
        best = min((minimize(cost, [np.log10(a0), l0], method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-8))
                    for l0 in (0.3, 0.6, 1.2)), key=lambda r: r.fun)
        a, lam = 10 ** best.x[0], 10 ** best.x[1]
        sc = {s or 'all': RUN.sparc_score(gals, fn(a, lam), s)[0] for s in (None, 'test', 'validation')}
        out[rule] = dict(rms_at_round9_constants=at0[0], a_SI=float(a * L.KMS2_PER_KPC), lam=float(lam),
                         g_d_SI=float(a * lam * L.KMS2_PER_KPC), rms_refitted=sc, msq_refitted=float(best.fun))
    return out


# --------------------------------------------------------------------------- D. X-COP
def xcop(law):
    import common as C
    import run_v3 as R3
    ctx = C.Context(tier='quick', verbose=False)
    cls = ctx.xcop(); a, lam = law['a_code'], law['lam']

    def pred(rule):
        def f(c, u):
            k = L.heat_weight(np.sqrt(c['sig2_star_hse']), u)
            gN = G * c['Mb'] / c['Rk'] ** 2
            S = G * (c['W'] @ (k * c['dms'])) / c['Rk'] ** 2
            SN = G * (c['W'] @ (c['dms'] + c['dmg'])) / c['Rk'] ** 2
            ghot = G * np.array([np.sum((k * c['dms'])[c['s'] < R]) for R in c['Rk']]) / c['Rk'] ** 2
            ex = dict(round3=np.sqrt(a * (gN + S)), one_stream=np.sqrt(a * (SN + S)),
                      independent=np.sqrt(a) * (gN + ghot) / np.sqrt(SN + S))[rule]
            return (gN + L.released(gN, a, lam) * ex) * c['Rk'] ** 2 / G
        return f

    out = {}
    for rule in ('round3', 'independent', 'one_stream'):
        p = pred(rule)
        r = minimize_scalar(lambda lu: np.mean(R3.resid(cls, lambda c: p(c, 10 ** lu)) ** 2), bounds=(0.5, 3.5), method='bounded')
        u = 10 ** r.x; res = R3.resid(cls, lambda c: p(c, u))
        out[rule] = dict(u_best=float(u), rms=float(np.sqrt(np.mean(res ** 2))),
                         mean_residual_by_radius=[float(x) for x in res.mean(0)])
    c0 = cls[0]
    out['heat_over_net_pull_by_radius'] = [float(x) for x in np.mean([
        (G * (c['W'] @ (L.heat_weight(np.sqrt(c['sig2_star_hse']), law['u_kms']) * c['dms'])) / c['Rk'] ** 2) / (G * c['Mb'] / c['Rk'] ** 2)
        for c in cls], axis=0)]
    out['radii_over_R500'] = [float(x) for x in c0['Rk'] / c0['R5']]
    return out


# --------------------------------------------------------------------------- E. Milky Way
def milky_way(fits):
    """The local estimate of part 1 (cold matter only), for each rule with its own SPARC-fitted a and g_d."""
    import first_principles_v10 as FP
    rows = FP.milky_way_factor()
    for r in rows:
        gN, SN, R = r['g_N'], r['S_N'], r['R_kpc']
        for rule, f in fits.items():
            a, gd = f['a_SI'] / L.KMS2_PER_KPC, f['g_d_SI'] / L.KMS2_PER_KPC
            ex = dict(round3=np.sqrt(a * gN), one_stream=np.sqrt(a * SN), independent=np.sqrt(a) * gN / np.sqrt(SN))[rule]
            r[f'v_{rule}_own_constants'] = float(np.sqrt(R * (gN + np.exp(-gN / gd) * ex)))
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir; out.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic(); rng = np.random.default_rng(20260925)
    from law_config import load_law
    law = load_law('round9')
    res = dict(experiment='round 10, part 2: how the companion adds up (three rules, galaxies and clusters)',
               law=dict(a_SI=law['a_SI'], lam=law['lam'], u_kms=law['u_kms']))

    A = part_a(rng); res['A_local_pull'] = A
    print('A. pull = (q0/2) Im(Psi* grad Psi)/|Psi|: largest difference from the exact time average',
          ', '.join(f"{r['config']} {r['max_abs_difference']:.1e} (force ~{r['typical_force']:.1e})" for r in A['identity_checks']))
    print('   two coherent waves in opposite directions (amplitudes 1 and b): mean pull / one wave\'s pull')
    for r in A['counter_propagating']:
        print(f"     b {r['counter_amplitude']:.2f}: {r['mean_pull_over_single_wave']:.3f} (closed form {r['closed_form']:.3f});"
              f" net flux {r['net_flux_over_single_wave']:.3f}, round-3 pull {r['round3_pull_over_single_wave']:.3f}, energy {r['energy_over_single_wave']:.2f}")

    print('C. SPARC ...', flush=True)
    gals = RUN.load_sparc()
    res['C_sparc'] = []
    for h_gas, zfac in ((0.2, 1.0), (0.5, 2.0)):
        info = sparc_plain_totals(gals, h_gas, zfac)
        info['fits'] = sparc_fits(gals, law)
        res['C_sparc'].append(info)
        print(f"   gas thickness {h_gas} kpc, stars {info['stellar_thickness']}: reconstruction of SPARC's pulls, disk "
              f"{info['disk_pull_model_over_sparc_16_50_84'][1]:.3f}, gas {info['gas_pull_model_over_sparc_16_50_84'][1]:.3f};"
              f" S_N/|g_N| median {info['SN_over_gN_16_50_84'][1]:.2f} (16-84%: {info['SN_over_gN_16_50_84'][0]:.2f}-{info['SN_over_gN_16_50_84'][2]:.2f})", flush=True)
        for rule, f in info['fits'].items():
            print(f"     {rule:12s} at the round-9 constants {f['rms_at_round9_constants']:.2f} km/s; refitted a {f['a_SI']:.3g} m/s2,"
                  f" g_d {f['g_d_SI']:.3g}: {f['rms_refitted']['all']:.2f} km/s (test {f['rms_refitted']['test']:.2f},"
                  f" validation {f['rms_refitted']['validation']:.2f})", flush=True)

    print('D. X-COP ...', flush=True)
    res['D_xcop'] = D = xcop(law)
    for rule in ('round3', 'independent', 'one_stream'):
        r = D[rule]
        print(f"   {rule:12s} best u {r['u_best']:6.1f} km/s  rms {r['rms']:.3f}  mean residual by radius " + ' '.join(f'{x:+.2f}' for x in r['mean_residual_by_radius']))
    print('   heat term / net Newtonian pull by radius: ' + ' '.join(f'{x:.1f}' for x in D['heat_over_net_pull_by_radius']))

    res['E_milky_way'] = E = milky_way(res['C_sparc'][0]['fits'])
    for r in E:
        print(f"E. Milky Way R {r['R_kpc']:5.1f} kpc: S_N/|g_N| {1 / r['ratio_gN_over_SN']:.2f}; v (local estimate, each rule with its own constants)"
              f" round 3 {r['v_round3_own_constants']:.1f}, independent {r['v_independent_own_constants']:.1f}, one stream {r['v_one_stream_own_constants']:.1f} km/s"
              f" (Newton {r['v_newton']:.1f})")

    base = res['C_sparc'][0]['fits']
    res['summary'] = {rule: dict(sparc_rms_kms=base[rule]['rms_refitted']['all'], xcop_rms=D[rule]['rms'],
                                 sun_local_kms=[r for r in E if r['R_kpc'] == 8.2][0][f'v_{rule}_own_constants'])
                      for rule in ('round3', 'independent', 'one_stream')}
    res['summary']['benchmarks'] = dict(sparc_mond_simple_kms=16.1328, note='SPARC: typical speed miss after refitting a and g_d; X-COP: rms of ln(M_obs/M_pred) after refitting u')
    res['seconds'] = time.monotonic() - t0
    (out / 'combination_rules_v10.json').write_text(json.dumps(res, indent=1, default=float) + '\n')
    print(f"wrote {out / 'combination_rules_v10.json'} ({res['seconds']:.0f} s)")


if __name__ == '__main__':
    main()

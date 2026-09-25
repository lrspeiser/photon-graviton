"""Round 19: the blocker mechanism for screening and release (supplied 25 September 2026), carried to the Solar System.

The derivation supplied by the user (research_work/results/hot-companion/screening-blockers-v19/, checked here: its
master-equation and event checks reproduce) proposes one elementary process for both of the law's factors: a
coupling element of the companion holds n temporary blockers, created at rate eta g and each removed at rate gamma; it
couples only with none. Then the equilibrium open fraction is exp(-g/g_d), g_d = gamma/eta, and, if a newly emitted
element inherits the local equilibrium population plus one launch blocker, its release in a constant field is
exp(-g/g_d)(1 - exp(-r/L)), L = u/gamma: the adopted law's two factors, from one lifetime. In a changing field the
mean load m obeys dm/dt = eta g - gamma m, and the open fraction is (1 - s) exp(-m) (s: the launch blocker survives).

What that means around a star. A star's companion is launched in the star's own field, g/g_d ~ 10^12 at the Sun's
surface, and then moves outward at u, where the field falls faster than the blockers can clear (for r << L). The load
it carries out is cleared only after about ln(load) lifetimes, so the release becomes a sharp switch far beyond L. This
script computes the Sun's release profile along its outward path for three launch preparations, each against the
adopted scalar release (1 - exp(-r/L)), and the resulting Cassini quadrupole Q2 and wide-binary boosts with the same
exact solver as the regression suite (regression/t_precision.sun_in_galaxy, here with the release profile as input):

  equilibrium    the element inherits the local equilibrium load at the surface, x = g/g_d, plus one launch blocker
                 (the supplied derivation's preparation)
  path only      it starts empty (plus one launch blocker) at the surface and gathers blockers on the way out
  capacity K     K sites per element, filled at launch; the open fraction is (1 - p)^K per element, each site emptying
                 at gamma (the supplied finite-capacity control), K = 10 and 100

With the release written as the excess over the local equilibrium, R_b = (1 - s) exp(-(m - g/g_d)), a constant field
gives back the adopted 1 - exp(-r/L) exactly (checked). The Galaxy's field at the Sun is taken from the round-19 Milky
Way fit (code/mw_joint_v19.py) and, for comparison, from the suite's inversion of 230 km/s.

    python code/blocker_release_v19.py --output run-blocker-release-v19/blocker_release_v19.json
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / 'regression'))
from t_precision import G, AU, MSUN, galactic_pull                     # noqa: E402

PC_AU = 206264.806
R_SUN = 6.957e8
Q2_MEAS = (1.6e-27, 1.8e-27)                                           # Park et al. 2026


def sun_in_galaxy_R(M, ge, S, H, a, gd, R_of_r, rmin_au=10.0, rmax_au=1e9, n_r=8000, n_mu=800, seps_au=(3000, 7000, 20000)):
    """regression/t_precision.sun_in_galaxy with the release profile R(r) given as a function (its derivative taken
    numerically): the same phantom density rho = R rho0 - R' D w_p / (4 pi G)."""
    mu, w = np.polynomial.legendre.leggauss(n_mu)
    r = np.geomspace(rmin_au, rmax_au, n_r) * AU
    gi = G * M / r ** 2
    st = np.sqrt(1 - mu ** 2)
    gp = ge * mu[None, :] - gi[:, None]; gq = ge * st[None, :]
    g = np.sqrt(gp ** 2 + gq ** 2)
    wp = gp + H * mu[None, :]; wq = gq + H * st[None, :]
    w_dot_gh = (wp * gp + wq * gq) / g
    F = np.exp(-g / gd) * np.sqrt(a * (g + S)); dF = F * (-1 / gd + 0.5 / (g + S))
    D = F / (g + H); dD = (dF * (g + H) - F) / (g + H) ** 2
    rho0 = (M / (4 * np.pi * r[:, None] ** 3)) * dD * (w_dot_gh - 3 * wp * gp / g)
    R = R_of_r(r)
    dR = np.gradient(R, r)
    rho = R[:, None] * rho0 - dR[:, None] * D * wp / (4 * np.pi * G)
    P2 = 0.5 * (3 * mu ** 2 - 1)
    Q2 = 3 * G * np.trapezoid(2 * np.pi * (rho * P2[None, :] * w[None, :]).sum(1) / r, r)
    iso = 2 * np.pi * (rho * w[None, :]).sum(1) * r ** 2
    Mph = np.concatenate([[0], np.cumsum(0.5 * (iso[1:] + iso[:-1]) * np.diff(r))])
    boost = {k: float(1 + np.interp(k * AU, r, Mph) / M) for k in seps_au}
    return float(Q2), boost


def blocker_profile(L_m, gd, ge, prep, K=None, n=200000, M=MSUN):
    """The release R_b(r) along the Sun's outward path (r in m): the load's excess over the local equilibrium, times the
    launch blocker's clearance. The field magnitude along the path is sqrt((GM/r^2)^2 + ge^2) (its average over
    directions, squared); x = g/g_d. The mean load obeys dm/dr = (x - m)/L; for K sites the occupied share p obeys
    dp/dr = ((x/K)(1 - p) - p)/L."""
    rr = np.geomspace(R_SUN, 1e9 * AU, n)
    x = np.sqrt((G * M / rr ** 2) ** 2 + ge ** 2) / gd
    s = np.exp(-(rr - R_SUN) / L_m)                                   # the launch blocker survives
    if prep == 'scalar':
        Rb = 1 - s
        return lambda r: np.interp(r, rr, Rb)
    # exact update over each step for piecewise-constant x: m -> x + (m - x) e^{-dr/L}
    dr = np.diff(rr); e = np.exp(-dr / L_m)
    if prep in ('equilibrium', 'path only'):
        m = np.empty_like(rr); m[0] = x[0] if prep == 'equilibrium' else 0.0
        xm = 0.5 * (x[1:] + x[:-1])
        for i in range(len(dr)):
            m[i + 1] = xm[i] + (m[i] - xm[i]) * e[i]
        excess = np.maximum(m - x, 0.0)
        Rb = (1 - s) * np.exp(-excess)
    elif prep == 'capacity':
        p = np.empty_like(rr); p[0] = 1.0                              # filled at launch
        # dp/dr = (x/K - p (1 + x/K))/L: relaxes to p_eq = x/(K + x) at rate (1 + x/K)/L
        for i in range(len(dr)):
            xi = 0.5 * (x[i] + x[i + 1]); rate = (1 + xi / K) / L_m; peq = xi / (K + xi)
            p[i + 1] = peq + (p[i] - peq) * np.exp(-rate * dr[i])
        peq_all = x / (K + x)
        Rb = (1 - s) * ((1 - p) / (1 - peq_all)) ** K                # open fraction over its equilibrium value
    else:
        raise ValueError(prep)
    return lambda r: np.interp(r, rr, Rb)


def switch_radius(Rfun):
    """Where the release reaches one half (AU)."""
    r = np.geomspace(10, 1e8, 4000) * AU
    R = Rfun(r)
    i = np.argmax(R >= 0.5)
    return float(r[i] / AU) if R[i] >= 0.5 else None


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    from law_config import load_law
    law = load_law('round12'); a, gd = law['a_SI'], law['g_d_SI']
    L0 = law['release_length_au'] * AU
    mw = json.loads((HERE.parent / 'run-mw-joint-v19/mw_joint_v19.json').read_text())['best']['ours']['sun']
    ge_suite, _ = galactic_pull(law)
    fields = {'round-19 Milky Way fit': (mw['g_N_SI'], mw['S_SI'], mw['g_hot_SI']), 'suite (230 km/s, no heat)': (ge_suite, 0.0, 0.0)}
    preps = [('scalar', None), ('equilibrium', None), ('path only', None), ('capacity', 10), ('capacity', 100)]
    out = dict(experiment='round 19: the blocker mechanism for screening and release, carried to the Solar System',
               law=dict(a_SI=a, g_d_SI=gd, L_pc=L0 / AU / PC_AU, u_kms=law['u_kms']), fields=fields, results={})
    # check: in a constant field the blocker release with the equilibrium preparation equals the scalar one
    Rc = blocker_profile(L0, gd, 3.0 * gd, 'equilibrium', M=0.0)         # a constant field: no Sun
    rr = np.geomspace(1e3, 1e8, 50) * AU
    out['check_constant_field_max_diff'] = float(np.max(np.abs(Rc(rr) - (-np.expm1(-(rr - R_SUN) / L0)))))
    print('check (constant field): blocker release minus 1 - exp(-r/L), max', out['check_constant_field_max_diff'], flush=True)
    for fname, (ge, S, H) in fields.items():
        rows = []
        for prep, K in preps:
            lab = prep if K is None else f'capacity K = {K}'
            Rf = blocker_profile(L0, gd, ge, prep, K)
            Q2, b = sun_in_galaxy_R(MSUN, ge, S, H, a, gd, Rf)
            row = dict(preparation=lab, lifetime_length_pc=L0 / AU / PC_AU, half_release_au=switch_radius(Rf), Q2=Q2,
                       z_2026=(Q2 - Q2_MEAS[0]) / Q2_MEAS[1], boost_3000=b[3000], boost_7000=b[7000], boost_20000=b[20000])
            # the lifetime that gives the adopted law's Q2 (the same Cassini standing), and its binaries
            if prep != 'scalar':
                Qtarget = rows[0]['Q2']
                f = lambda lnL: sun_in_galaxy_R(MSUN, ge, S, H, a, gd, blocker_profile(np.exp(lnL) * AU, gd, ge, prep, K), n_r=3000, n_mu=300)[0] / Qtarget - 1
                try:
                    lnL = brentq(f, np.log(10.0), np.log(L0 / AU), xtol=1e-3)
                    Lm = np.exp(lnL) * AU; Rf2 = blocker_profile(Lm, gd, ge, prep, K)
                    Q2b, b2 = sun_in_galaxy_R(MSUN, ge, S, H, a, gd, Rf2)
                    row['same_Q2'] = dict(lifetime_length_pc=float(Lm / AU / PC_AU), lifetime_years=float(Lm / (law['u_kms'] * 1e3) / 3.15576e7),
                                          half_release_au=switch_radius(Rf2), Q2=Q2b, boost_3000=b2[3000], boost_7000=b2[7000], boost_20000=b2[20000])
                except ValueError:
                    row['same_Q2'] = None
            rows.append(row)
            s2 = row.get('same_Q2')
            print(f"[{time.monotonic() - t0:4.0f} s] {fname:28s} {lab:18s}: half release at {row['half_release_au'] or float('nan'):10.0f} AU; Q2 {Q2:.2e} "
                  f"({row['z_2026']:+.2f} sigma); binaries +{100 * (b[7000] - 1):.2f}% / +{100 * (b[20000] - 1):.2f}% at 7,000 / 20,000 AU"
                  + (f" | same Q2 as scalar at lifetime length {s2['lifetime_length_pc']:.4f} pc ({s2['lifetime_years']:.0f} yr): half release "
                     f"{s2['half_release_au']:.0f} AU, binaries +{100 * (s2['boost_7000'] - 1):.2f}% / +{100 * (s2['boost_20000'] - 1):.2f}%" if s2 else ''), flush=True)
        out['results'][fname] = rows
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('wrote', args.output)


if __name__ == '__main__':
    main()

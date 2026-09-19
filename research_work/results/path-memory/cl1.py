"""CL-1: the written-track response at cluster scale and for light (protocol-cl1.md).

    python cl1.py [--output-dir DIR] [--canonical]

Runs the nine numerical gates with their negative controls, then the declared exploratory scan -- E1 the
X-COP clusters, E2 the six SLACS lenses under the light rule L1, E3 Coma under both rules, E4 the
compatibility map -- and writes cl1-results.json, comparing it with the archive (the first run creates it).
Exits non-zero if any gate fails, any control is not rejected, or any archived number moves. About four
minutes on one core; cl1_checks.py is the suite job.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
sys.path.insert(0, str(RESULTS/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(RESULTS/'radiation-polarized-gravity'))
import lensing as LZ  # noqa: E402
sys.path.insert(0, str(HERE))
import steady_field as SF  # noqa: E402
import rut1 as U  # noqa: E402

G = SF.G
INPUTS = HERE/'cl1-inputs-xcop.json'
KUBO = RESULTS/'cluster-observation-readiness/kubo-figure-data.json'
RPG_ARCHIVE = RESULTS/'radiation-polarized-gravity/rpg1-results.json'
H_COMA = .7                                   # adopted: the X-COP papers' H0 = 70, for the Kubo radii only
CLUSTER_WIDTHS = np.geomspace(10., 1e4, 41)   # kpc, E1
LENS_WIDTHS = np.geomspace(.1, 1e4, 31)       # kpc, E2 and the map
A_STAR_REPO, A_STAR_PM1 = 8.563335e-11, 6.54e-11      # m/s^2, the two reference scales of PM-1's law
FGAS_EARLIER = dict(A85=.150, A644=.132, A1644=.128, A1795=.139, A2142=.158, A2255=.153, A2319=.189, A3158=.145,
                    A3266=.132, RXC1825=.133, ZW1215=.106)   # companion_wave_test/cluster_test.py, fgas500
COMA_BRACKET = ((2.5e-3, .5e13), (4.5e-3, 2e13))            # CF-1's declared (n_e0 cm^-3, M_star Msun)
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:6.1f}s] {msg}', flush=True)


def d4(f, x, h):
    """Fourth-order central difference."""
    return (f(x - 2*h) - 8*f(x - h) + 8*f(x + h) - f(x + 2*h))/(12*h)


def scale_rel(a, b):
    """max |a - b| / max |b|: the relative difference on the scale of the quantity, defined where b crosses zero."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.max(np.abs(a - b))/np.max(np.abs(b)))


# ================================================================ the gates
def gate_g1():
    cases = ((1., 1., .1), (1.3, 1., .1), (.7, 1., .1), (2., 1., .3), (.5, 2., .25), (5., 1., .5), (1e-3, 1., .2))
    rows, ctrl = [], []
    for r, rp, w in cases:
        num = .5*quad(lambda mu: np.exp(-(r*r + rp*rp - 2*r*rp*mu)/(2*w*w)), -1, 1, epsabs=1e-30, epsrel=1e-13)[0]
        closed = SF.shell_kernel([r], [rp], w)[0, 0]
        planar = SF.ring_kernel([r], [rp], w)[0, 0]
        rows.append(dict(r=r, r_prime=rp, w=w, quadrature=num, closed_form=closed,
                         relative_difference=abs(closed/num - 1)))
        ctrl.append(abs(planar/num - 1))
    worst = max(x['relative_difference'] for x in rows)
    return dict(rows=rows, worst_relative_difference=worst, tolerance=1e-12, passed=bool(worst < 1e-12),
                control=dict(name='RUT-1 planar ring kernel in place of the shell kernel',
                             worst_relative_difference=max(ctrl), rejected=bool(max(ctrl) > 1e-3)))


def _plummer_source(n=4000, r_max=60.):
    r = np.geomspace(1e-3, r_max, n)
    return SF.SphericalSource(r, SF.plummer_mass(1e12, 5., r)), r


def gate_g2():
    src, r = _plummer_source()
    radii = np.array([.5, 2., 5., 10., 20.])
    rows, worst_g, worst_l, ctrl = [], 0., 0., 1.
    for w in (.5, 2., 8.):
        h = 1e-3*min(w, radii.min())
        grad_an = src.written_gradient(radii, w)
        grad_fd = np.array([d4(lambda x: src.written_potential([x], w)[0], x, h) for x in radii])
        lap_an = src.written_laplacian(radii, w)
        lap_fd = np.array([d4(lambda x: src.written_gradient([x], w)[0], x, h) for x in radii]) + 2*grad_an/radii
        eg, el = scale_rel(grad_an, grad_fd), scale_rel(lap_an, lap_fd)
        # the control: the local-limit gradient (2 pi)^{3/2} w^3 rho'(r)
        rho_p = d4(lambda x: SF.plummer_density(1e12, 5., x), radii, 1e-3*radii)
        local = (2*np.pi)**1.5*w**3*rho_p
        ec = scale_rel(local, grad_an)
        worst_g, worst_l, ctrl = max(worst_g, eg), max(worst_l, el), min(ctrl, ec)
        rows.append(dict(w=w, gradient_scale_relative=eg, laplacian_scale_relative=el,
                         local_limit_scale_relative=ec))
    return dict(rows=rows, worst_gradient=worst_g, worst_laplacian=worst_l, tolerance=1e-8,
                passed=bool(worst_g < 1e-8 and worst_l < 1e-8),
                control=dict(name='local-limit force -(2 pi)^{3/2} w^3 Lambda rho\' in place of the full gradient',
                             smallest_scale_relative_difference=ctrl, rejected=bool(ctrl > 1e-8)),
                note='scale-relative: max |a - b| / max |b| over the five radii, because the Laplacian crosses zero')


def gate_g3():
    cases = ((1., 1., .1), (1.3, 1., .1), (.7, 1., .1), (2., 1., .3), (.5, 2., .25))
    same = max(abs(SF.ring_kernel([r], [R], w)[0, 0]/(-U.phi_ring(r, R, w)) - 1) for r, R, w in cases)
    k, m = U.gate_kernel(), U.gate_mature_ring()
    A = .2/U.ring_depth(1., .1, 1.)                     # the control's writing amplitude q tau at G = M = 1
    w_over_ell = .1*A
    # (d) the ring kernel's second derivative
    Rp = np.array([1., 3., 7.])
    worst, ctrl = 0., 1.
    for w in (.5, 2.):
        for R in (.2, 1., 2.5, 6.):
            h = 1e-3*min(w, R)
            fd = d4(lambda x: SF.ring_kernel_dR([x], Rp, w)[0], R, h)
            an = SF.ring_kernel_d2R([R], Rp, w)[0]
            bad = SF.ring_kernel_d2R([R], Rp, w, drop_cross_term=True)[0]
            worst = max(worst, scale_rel(an, fd))
            ctrl = min(ctrl, scale_rel(bad, fd))
    return dict(ring_kernel_equals_rut1_phi_ring=same, rut1_gate_kernel_passed=bool(k['passed']),
                rut1_gate_mature_ring_passed=bool(m['passed']),
                mature_ring_control=dict(support_percent=m['percent_of_newtonian'], amplitude_q_tau=A,
                                         w_over_ell=w_over_ell, declared=.501),
                second_derivative_worst=worst, tolerance=1e-8,
                passed=bool(same < 1e-14 and k['passed'] and m['passed'] and abs(w_over_ell - .501) < 5e-4
                            and worst < 1e-8),
                control=dict(name='ring kernel second derivative with the cross term 2 E\' I0e\' dropped',
                             smallest_scale_relative_difference=ctrl, rejected=bool(ctrl > 1e-8)))


def gate_g4():
    src, r = _plummer_source()
    M_t = src.M[-1]
    rows, worst, ctrl = [], 0., 1.
    for w in (.5, 2., 8.):
        r_out = 60. + 8*w
        m_eff = src.equivalent_mass([r_out], w)[0]/(M_t/G)
        rr = np.geomspace(1e-3, r_out, 20000)
        rho_eff = src.equivalent_density(rr, w)
        integral = np.trapezoid(4*np.pi*rr*rr*rho_eff, rr)/(M_t/G)
        # Gauss consistency between the two kernels at interior radii
        interior = []
        for r_i in (5., 20.):                      # the integral ends exactly at r_i, not at the grid point below it
            ri = np.append(rr[rr < r_i], r_i)
            cum = np.trapezoid(4*np.pi*ri*ri*np.append(rho_eff[rr < r_i], src.equivalent_density([r_i], w)[0]), ri)
            interior.append(abs(cum/src.equivalent_mass([r_i], w)[0] - 1))
        newton = src.mass(r_out)/M_t
        rows.append(dict(w=w, r_out=r_out, M_eff_over_LambdaM_G=m_eff, integral_rho_eff_over_LambdaM_G=integral,
                         gauss_consistency_at_5_and_20=interior, newtonian_control=newton))
        worst = max(worst, abs(m_eff), abs(integral), max(interior))
        ctrl = min(ctrl, newton)
    return dict(rows=rows, worst=worst, tolerances=dict(M_eff=1e-8, integral=1e-6, gauss=1e-6),
                passed=bool(all(abs(x['M_eff_over_LambdaM_G']) < 1e-8 and abs(x['integral_rho_eff_over_LambdaM_G']) < 1e-6
                                and max(x['gauss_consistency_at_5_and_20']) < 1e-6 for x in rows)),
                control=dict(name='the Newtonian equivalent mass of the same source at r_out, which is M',
                             smallest=ctrl, rejected=bool(ctrl > 1e-8)),
                statement='the written field of a bounded source carries no net equivalent mass')


def _g5_routes(n3, n2, w, b):
    r3 = np.geomspace(1e-4, 300., n3)
    src3 = SF.SphericalSource(r3, SF.plummer_mass(1e12, 5., r3))
    R2 = np.geomspace(1e-3, 300., n2)
    proj = SF.ProjectedSource(R2, SF.plummer_projected_mass(1e12, 5., R2))
    a3 = SF.deflection_from_g(lambda x: src3.g_mem([x], w)[0], b)
    a2 = SF.deflection_from_projected_mass(proj.projected_equivalent_mass([b], w)[0], b)
    ac = SF.deflection_from_projected_mass(proj.projected_equivalent_mass([b], w, kernel_dR=SF.shell_kernel_dr)[0], b)
    return a3, a2, ac


def gate_g5():
    rows, worst, ctrl, fell = [], 0., 1., True
    for w in (1., 3., 10.):
        for b in (1., 4., 10., 25.):
            a3, a2, ac = _g5_routes(8000, 6000, w, b)
            a3d, a2d, _ = _g5_routes(16000, 12000, w, b)
            d, dd = abs(a3/a2 - 1), abs(a3d/a2d - 1)
            rows.append(dict(w=w, b=b, three_dimensional=a3, two_dimensional=a2, relative_difference=d,
                             relative_difference_doubled_grids=dd, control_relative_difference=abs(ac/a3 - 1)))
            worst, ctrl = max(worst, d), min(ctrl, abs(ac/a3 - 1))
            fell = fell and dd < d
    return dict(rows=rows, worst_relative_difference=worst, tolerance=1e-6, difference_falls_when_grids_doubled=fell,
                passed=bool(worst < 1e-6 and fell),
                control=dict(name='the two-dimensional route with the shell kernel in place of the ring kernel',
                             smallest_relative_difference=ctrl, rejected=bool(ctrl > 1e-6)))


class Lens:
    """One SLACS lens through RPG-1's inputs: deprojected light for the three-dimensional route, analytic
    projected light for the two-dimensional one."""

    def __init__(self, name, data, n_proj=8000):
        obs, geo, light, masses = data
        o, g = obs[name], geo[name]
        self.name, self.theta_obs = name, o['bSIE']
        self.Dl, self.ratio = g['conditional_Dl_Mpc']*1000, g['conditional_Dls_over_Ds']
        self.r, self.frac, self.parts = LZ.light_profile(light[name]['components'], self.Dl)
        self.lr = np.log(self.r)
        self.components = light[name]['components']
        self.Re = [c['R_arcsec']*self.Dl/SF.ARCSEC for c in self.components]
        self.Re_main = self.Re[int(np.argmax(self.parts))]
        self.masses = {imf: 10**masses[(name, imf)] for imf in ('Chabrier', 'Salpeter')}
        Rg = np.geomspace(1e-4, 1e5, n_proj)
        self.proj = SF.ProjectedSource(Rg, SF.sersic_projected_fraction(self.components, self.parts, self.Dl, Rg))

    def proj_frac(self, b):
        return SF.sersic_projected_fraction(self.components, self.parts, self.Dl, np.atleast_1d(b))[0]

    def theta_3d_newton(self, M, epsrel=1e-8, limit=200):
        """RPG-1's own quadrature settings, so the archive is reproduced exactly."""
        gN = lambda x: G*M*np.interp(np.log(x), self.lr, self.frac, left=0., right=self.frac[-1])/(x*x)
        bend = lambda b: 4/SF.C_KMS**2*quad(lambda t: gN(b/np.cos(t))*b/np.cos(t), 0, np.pi/2,
                                            limit=limit, epsabs=0, epsrel=epsrel)[0]
        return SF.einstein_radius(bend, self.Dl, self.ratio)

    def theta_2d(self, M, lam=0., w=1.):
        """theta_E under L1 with amplitude lam (0 is L0 = Newtonian) through the projected route."""
        def alpha(b):
            Mp = M*self.proj_frac(b)
            if lam:
                Mp = Mp + lam*M*self.proj.projected_equivalent_mass([b], w)[0]
            return SF.deflection_from_projected_mass(Mp, b)
        return SF.einstein_radius(alpha, self.Dl, self.ratio)

    def lambda_for_observed(self, M, w):
        """The Lambda giving the observed theta_E at width w: theta rises monotonically with Lambda."""
        lo, hi = 0., 1e-9
        while self.theta_2d(M, hi, w) < self.theta_obs:
            hi *= 4
            if hi > 1e6:
                return np.nan
        return brentq(lambda L: self.theta_2d(M, L, w) - self.theta_obs, lo, hi, xtol=1e-16, rtol=1e-10)


def load_lenses(n_proj=8000):
    data = LZ.load()
    return [Lens(name, data, n_proj) for name in LZ.LENSES]


def gate_g6(lenses):
    arch = {r['name']: r for r in json.loads(RPG_ARCHIVE.read_text(encoding='utf-8'))['T3']['lenses']}
    rows, worst3, worst2, ctrl = [], 0., 0., 1.
    for L in lenses:
        for imf in ('Chabrier', 'Salpeter'):
            M = L.masses[imf]
            a = arch[L.name][imf]['cases']['newtonian_baryons']['theta_arcsec']
            t3, t2 = L.theta_3d_newton(M), L.theta_2d(M)
            moved = L.theta_2d(M, 1e-6, 3.)
            rows.append(dict(lens=L.name, imf=imf, archived=a, three_dimensional=t3, two_dimensional=t2,
                             rel_3d=abs(t3/a - 1), rel_2d=abs(t2/a - 1), theta_with_lambda_1e_minus_6_w_3kpc=moved))
            worst3, worst2, ctrl = max(worst3, abs(t3/a - 1)), max(worst2, abs(t2/a - 1)), min(ctrl, abs(moved/a - 1))
    return dict(rows=rows, worst_3d=worst3, worst_2d=worst2, tolerances=dict(three_dimensional=1e-10, two_dimensional=1e-5),
                passed=bool(worst3 < 1e-10 and worst2 < 1e-5),
                control=dict(name='Lambda = 1e-6 (km/s)^2/Msun at w = 3 kpc must move every theta_E by more than 1e-3',
                             smallest_move=ctrl, rejected=bool(ctrl > 1e-3)))


def nfw_mass(M200, c, R200, r):
    m = lambda x: np.log1p(x) - x/(1 + x)
    return M200*m(np.asarray(r, float)/(R200/c))/m(c)


def gate_g7(inputs, swap=None):
    T1, cosmo = inputs['ettori_table1'], inputs['cosmology']
    H0 = cosmo['H0_km_s_Mpc']/1e3
    c200 = {n: v['c200'] for n, v in T1.items()}
    if swap:
        c200[swap[0]], c200[swap[1]] = c200[swap[1]], c200[swap[0]]
    rows, worst_m, worst_r = [], 0., 0.
    for n, v in T1.items():
        R200, R500 = v['R200_Mpc']*1e3, v['R500_Mpc']*1e3
        M = lambda r: nfw_mass(v['M200'], c200[n], R200, r)
        rho_c = 3*(H0**2*(cosmo['Omega_m']*(1 + v['z'])**3 + 1 - cosmo['Omega_m']))/(8*np.pi*G)
        R200_calc = (3*v['M200']*1e14/(4*np.pi*200*rho_c))**(1/3)
        R500_calc = (3*v['M500']*1e14/(4*np.pi*500*rho_c))**(1/3)
        dm = [M(500)/v['M_0p5Mpc'] - 1, M(1000)/v['M_1Mpc'] - 1, M(1500)/v['M_1p5Mpc'] - 1, M(R500)/v['M500'] - 1]
        dr = [R200_calc/R200 - 1, R500_calc/R500 - 1]
        rows.append(dict(cluster=n, mass_deviations=dm, radius_deviations=dr))
        worst_m, worst_r = max(worst_m, max(abs(x) for x in dm)), max(worst_r, max(abs(x) for x in dr))
    fgas_same = all(abs(inputs['eckert_table2'][n]['fgas_500'] - f) < 1e-12 for n, f in FGAS_EARLIER.items())
    out = dict(rows=rows, worst_mass_deviation=worst_m, worst_radius_deviation=worst_r,
               tolerances=dict(mass=5e-3, radius=1e-3), fgas500_agrees_with_earlier_transcription=fgas_same,
               passed=bool(worst_m < 5e-3 and worst_r < 1e-3 and fgas_same))
    if swap is None:
        ctrl = gate_g7(inputs, swap=('A85', 'A2255'))
        out['control'] = dict(name='A85 and A2255 concentrations swapped', worst_mass_deviation=ctrl['worst_mass_deviation'],
                              rejected=bool(ctrl['worst_mass_deviation'] > 5e-3))
    return out


def gate_g8(kubo):
    y = np.array([v['shear_t'] for v in kubo['rows']])
    e = np.array([v['plotted_sigma_t'] for v in kubo['rows']])
    chi = float(np.sum((y/e)**2))
    return dict(null_chi2=chi, paper=23.33, relative_difference=abs(chi/23.33 - 1), passed=bool(abs(chi/23.33 - 1) < 1e-3))


# ================================================================ E1: the clusters
def build_clusters(inputs, rc_frac=.15, fstar=.09, n_grid=3000):
    T1, T2 = inputs['ettori_table1'], inputs['eckert_table2']
    out = []
    for name, v in T1.items():
        if name not in T2:
            continue
        e = T2[name]
        R500, R200 = v['R500_Mpc']*1e3, v['R200_Mpc']*1e3
        Mg500, Mg200 = e['fgas_500']*e['M_HSE_500']*1e14, e['fgas_200']*e['M_HSE_200']*1e14
        r = np.geomspace(1., 2e4, n_grid)
        rc = rc_frac*R500
        ratio = lambda b: np.interp(R200, r, SF.beta_model_mass(r, 1., rc, b)[0])/np.interp(R500, r, SF.beta_model_mass(r, 1., rc, b)[0])
        beta = brentq(lambda b: ratio(b) - Mg200/Mg500, .3, 1.5, xtol=1e-12)
        M1 = SF.beta_model_mass(r, 1., rc, beta)[0]
        Mb = (1 + fstar)*Mg500/np.interp(R500, r, M1)*M1
        radii = np.array([500., 1000., 1500., R500, R200])
        Mreq = nfw_mass(v['M200']*1e14, v['c200'], R200, radii)
        err = np.array([v['M_0p5Mpc_err'], v['M_1Mpc_err'], v['M_1p5Mpc_err'], v['M500_err'], v['M200_err']])*1e14
        src = SF.SphericalSource(r, Mb)
        out.append(dict(name=name, source=src, radii=radii, R500=R500, R200=R200, beta=beta,
                        g_req=G*Mreq/radii**2, sigma=G*err/radii**2, g_N=src.g_newton(radii),
                        M_req=Mreq, M_b=src.mass(radii), fb_500=float(src.mass(R500)/(v['M500']*1e14))))
    return out


def cluster_shapes(cls, w):
    return [c['source'].g_mem(c['radii'], w) for c in cls]


def cluster_fit(cls, w, shapes=None):
    """Weighted linear least squares for Lambda at fixed w; returns Lambda, chi2, log-rms (dex)."""
    shapes = cluster_shapes(cls, w) if shapes is None else shapes
    num = sum(np.sum((c['g_req'] - c['g_N'])*s/c['sigma']**2) for c, s in zip(cls, shapes))
    den = sum(np.sum(s*s/c['sigma']**2) for c, s in zip(cls, shapes))
    lam = num/den
    chi2 = sum(np.sum(((c['g_N'] + lam*s - c['g_req'])/c['sigma'])**2) for c, s in zip(cls, shapes))
    model = np.concatenate([(c['g_N'] + lam*s)/c['g_req'] for c, s in zip(cls, shapes)])
    with np.errstate(invalid='ignore', divide='ignore'):
        logs = np.log10(model[model > 0])
    return float(lam), float(chi2), float(np.sqrt(np.mean(logs**2))) if len(logs) == len(model) else float('nan')


def cluster_chi2_fixed(cls, w, lam):
    shapes = cluster_shapes(cls, w)
    return float(sum(np.sum(((c['g_N'] + lam*s - c['g_req'])/c['sigma'])**2) for c, s in zip(cls, shapes)))


def cluster_scan(cls, widths):
    rows = [dict(w=float(w), Lambda=lam, ell=G/lam if lam else np.inf, w_over_ell=float(w*lam/G), chi2=chi2, log_rms_dex=lr)
            for w in widths for lam, chi2, lr in [cluster_fit(cls, w)]]
    best = min(rows, key=lambda x: x['chi2'])
    return rows, best


def cluster_references(cls):
    def chi2(law):
        return float(sum(np.sum(((law(c['g_N']) - c['g_req'])/c['sigma'])**2) for c in cls))
    def logrms(law):
        return float(np.sqrt(np.mean(np.concatenate([np.log10(law(c['g_N'])/c['g_req'])**2 for c in cls]))))
    laws = dict(newtonian_baryons=lambda g: g,
                pm1_law_repo_a_star=lambda g: g + np.sqrt(A_STAR_REPO*SF.SI_TO_CODE*g),
                pm1_law_fitted_a_star=lambda g: g + np.sqrt(A_STAR_PM1*SF.SI_TO_CODE*g))
    return {k: dict(chi2=chi2(f), log_rms_dex=logrms(f)) for k, f in laws.items()}


def e1(inputs):
    base = build_clusters(inputs)
    rows, best = cluster_scan(base, CLUSTER_WIDTHS)
    shapes = cluster_shapes(base, best['w'])
    per = [dict(cluster=c['name'], beta=c['beta'], baryon_fraction_R500=c['fb_500'],
                radii_kpc=c['radii'].tolist(), required_over_baryonic_mass=(c['M_req']/c['M_b']).tolist(),
                model_over_required=((c['g_N'] + best['Lambda']*s)/c['g_req']).tolist(),
                baryons_over_required=(c['g_N']/c['g_req']).tolist())
           for c, s in zip(base, shapes)]
    fac = np.concatenate([c['M_req']/c['M_b'] for c in base])
    sens = {}
    for label, kw in (('rc_0.10', dict(rc_frac=.10)), ('rc_0.25', dict(rc_frac=.25)),
                      ('fstar_0.07', dict(fstar=.07)), ('fstar_0.12', dict(fstar=.12))):
        _, b = cluster_scan(build_clusters(inputs, **kw), CLUSTER_WIDTHS)
        sens[label] = b
    # G9's cluster half: the base scan on a doubled radial grid
    rows_d, _ = cluster_scan(build_clusters(inputs, n_grid=6000), CLUSTER_WIDTHS)
    conv = max(max(abs(a['chi2']/b['chi2'] - 1), abs(a['Lambda']/b['Lambda'] - 1)) for a, b in zip(rows, rows_d))
    return dict(clusters=len(base), points=5*len(base), widths_kpc=CLUSTER_WIDTHS.tolist(), scan=rows, best=best,
                sensitivities=sens, per_cluster_at_best=per,
                required_over_baryonic_mass=dict(median=float(np.median(fac)), min=float(fac.min()), max=float(fac.max())),
                references=cluster_references(base), grid_doubling_max_relative_change=conv,
                reading='exploratory; exposed data; the five values of one cluster come from one NFW fit'), base, conv


# ================================================================ E2: the lenses under L1
def e2(lenses, cluster_best):
    arch = {r['name']: r for r in json.loads(RPG_ARCHIVE.read_text(encoding='utf-8'))['T3']['lenses']}
    out, conv = [], 0.
    doubled = {L.name: Lens(L.name, LZ.load(), n_proj=16000) for L in lenses}
    for L in lenses:
        row = dict(lens=L.name, theta_observed=L.theta_obs, Dl_kpc=L.Dl, Re_kpc=L.Re, Re_main_kpc=L.Re_main)
        for imf in ('Chabrier', 'Salpeter'):
            M = L.masses[imf]
            lam = np.array([L.lambda_for_observed(M, w) for w in LENS_WIDTHS])
            lam_d = np.array([doubled[L.name].lambda_for_observed(M, w) for w in LENS_WIDTHS])
            conv = max(conv, float(np.nanmax(np.abs(lam/lam_d - 1))))
            mono = all(L.theta_2d(M, .5*l, w) < L.theta_2d(M, l, w) < L.theta_2d(M, 2*l, w)
                       for l, w in zip(lam[::5], LENS_WIDTHS[::5]))
            w_over_ell = LENS_WIDTHS*lam/G
            i = int(np.nanargmin(w_over_ell))
            row[imf] = dict(
                log10_stellar_mass=float(np.log10(M)),
                theta_newtonian_L0=arch[L.name][imf]['cases']['newtonian_baryons']['theta_arcsec'],
                L0_ratio_to_observed=arch[L.name][imf]['cases']['newtonian_baryons']['ratio_to_observed'],
                Lambda=lam.tolist(), w_over_ell=w_over_ell.tolist(), w3_over_ell_kpc2=(LENS_WIDTHS**3*lam/G).tolist(),
                monotone_in_Lambda=bool(mono),
                min_w_over_ell=dict(value=float(w_over_ell[i]), w_kpc=float(LENS_WIDTHS[i]),
                                    w_over_Re_main=float(LENS_WIDTHS[i]/L.Re_main), Lambda=float(lam[i])),
                theta_at_cluster_best=L.theta_2d(M, cluster_best['Lambda'], cluster_best['w']),
                ratio_at_cluster_best=L.theta_2d(M, cluster_best['Lambda'], cluster_best['w'])/L.theta_obs)
        out.append(row)
        log(f'E2 {L.name} done')
    spread = {imf: [float(np.nanmax([r[imf]['Lambda'][k] for r in out])/np.nanmin([r[imf]['Lambda'][k] for r in out]))
                    for k in range(len(LENS_WIDTHS))] for imf in ('Chabrier', 'Salpeter')}
    return dict(widths_kpc=LENS_WIDTHS.tolist(), lenses=out, spread_max_over_min=spread,
                grid_doubling_max_relative_change=conv, geometry='PF-1 static Euclidean, as RPG-1',
                reading='exploratory; exposed data; L0 is the archived Newtonian bracket restated'), conv


# ================================================================ E3: Coma under both rules
def coma_model(ne0, mstar):
    r = np.geomspace(.1, 3000., 4000)
    rho0 = 1.17*SF.MP_G*ne0*(SF.KPC_M*100)**3/(SF.MSUN_KG*1000)
    Mg, rho_g = SF.beta_model_mass(r, rho0, 296., .75, r_trunc=3000.)
    rho = rho_g*(1 + mstar/Mg[-1])
    Rg = np.geomspace(.1, 3200., 1500)
    Sigma = SF.project_density(r, rho, Rg)
    Mp = cumulative_trapezoid(2*np.pi*Rg*Sigma, Rg, initial=0.) + np.pi*Rg[0]**2*Sigma[0]
    return dict(r=r, rho=rho, M_gas=float(Mg[-1]), M_b=float(Mg[-1]*(1 + mstar/Mg[-1])), Rg=Rg, Sigma=Sigma, Mp=Mp,
                proj=SF.ProjectedSource(Rg, Mp), projected_fraction=float(Mp[-1]/(Mg[-1]*(1 + mstar/Mg[-1]))))


def e3(kubo, cluster_best):
    Rk = np.array([v['published_radius_h_inverse_Mpc'] for v in kubo['rows']])*1e3/H_COMA
    y = np.array([v['shear_t'] for v in kubo['rows']])
    e = np.array([v['plotted_sigma_t'] for v in kubo['rows']])
    out = []
    for ne0, mstar in COMA_BRACKET:
        m = coma_model(ne0, mstar)
        Sb = np.interp(Rk, m['Rg'], m['Sigma'], right=0.)
        Mpb = m['proj'].projected_mass(Rk)
        dS0 = SF.delta_sigma(Mpb, Sb, Rk)
        A0 = float(np.sum(dS0*y/e**2)/np.sum((dS0/e)**2))
        chi0 = float(np.sum(((A0*dS0 - y)/e)**2))
        kap0 = A0*Sb
        red0 = float(np.max(np.abs(A0*dS0/(1 - kap0) - A0*dS0)/e))
        lam, w = cluster_best['Lambda'], cluster_best['w']
        Mpe = lam*m['proj'].projected_equivalent_mass(Rk, w)
        Se = lam*m['proj'].equivalent_surface_density(Rk, w)
        dS1 = dS0 + SF.delta_sigma(Mpe, Se, Rk)
        A1 = float(np.sum(dS1*y/e**2)/np.sum((dS1/e)**2))
        chi1 = float(np.sum(((A1*dS1 - y)/e)**2))
        kap1 = A1*(Sb + Se)
        red1 = float(np.max(np.abs(A1*dS1/(1 - kap1) - A1*dS1)/e))
        # the free two-nuisance fit: amplitude and Lambda, at the cluster-best width
        unit = SF.delta_sigma(m['proj'].projected_equivalent_mass(Rk, w), m['proj'].equivalent_surface_density(Rk, w), Rk)
        X = np.vstack([dS0, unit]).T/e[:, None]
        coef, *_ = np.linalg.lstsq(X, y/e, rcond=None)
        chif = float(np.sum((X@coef - y/e)**2))
        out.append(dict(n_e0_cm3=ne0, M_star=mstar, M_gas=m['M_gas'], M_b=m['M_b'], projected_fraction=m['projected_fraction'],
                        bins_kpc=Rk.tolist(),
                        L0=dict(delta_sigma_Msun_kpc2=dS0.tolist(), amplitude_inverse_sigma_crit=A0, implied_sigma_crit=1/A0,
                                shape_chi2=chi0, reduced_minus_weak_max_sigma=red0),
                        L1_at_cluster_best=dict(Lambda=lam, w=w, delta_sigma_Msun_kpc2=dS1.tolist(),
                                                ratio_L1_over_L0=(dS1/dS0).tolist(), amplitude_inverse_sigma_crit=A1,
                                                shape_chi2=chi1, reduced_minus_weak_max_sigma=red1),
                        free_fit=dict(w=w, Lambda=float(coef[1]/coef[0]), shape_chi2=chif, sign_of_Lambda=int(np.sign(coef[1]/coef[0])))))
    return dict(h_adopted=H_COMA, brackets=out, references=dict(nfw=3.855, plummer=3.727,
                note='the inverse fits have more freedom and are not a benchmark to beat'),
                reading='exploratory shape test; Sigma_crit is a nuisance; exposed figure-reconstructed bins')


# ================================================================ E4: the compatibility map
def e4(inputs, base_clusters, lens_rows, lenses, cluster_best):
    lam_c = [cluster_fit(base_clusters, w)[0] for w in LENS_WIDTHS]
    table = []
    for k, w in enumerate(LENS_WIDTHS):
        row = dict(w=float(w), Lambda_cluster=lam_c[k])
        for imf in ('Chabrier', 'Salpeter'):
            lams = [r[imf]['Lambda'][k] for r in lens_rows]
            allv = lams + [lam_c[k]]
            row[imf] = dict(Lambda_lenses=lams, lens_over_cluster=[l/lam_c[k] for l in lams],
                            max_over_min_including_cluster=float(np.nanmax(allv)/np.nanmin(allv)))
        table.append(row)
    compatible = {imf: [r['w'] for r in table if r[imf]['max_over_min_including_cluster'] < 2.] for imf in ('Chabrier', 'Salpeter')}
    lens_pairs = []
    for L, r in zip(lenses, lens_rows):
        mn = r['Chabrier']['min_w_over_ell']
        lens_pairs.append(dict(lens=L.name, w=mn['w_kpc'], Lambda=mn['Lambda'],
                               cluster_chi2_with_this_pair=cluster_chi2_fixed(base_clusters, mn['w_kpc'], mn['Lambda'])))
    return dict(widths_kpc=LENS_WIDTHS.tolist(), table=table, criterion='max/min of Lambda over the six lenses and the '
                'clusters below 2 at one w', widths_where_compatible=compatible,
                cluster_best_pair=dict(w=cluster_best['w'], Lambda=cluster_best['Lambda'], w_over_ell=cluster_best['w_over_ell']),
                lens_preferred_pairs_applied_to_clusters=lens_pairs,
                rut1_reference=dict(w_over_ell=.501, w_over_R=.1, support_percent=10.03,
                                    note='the mature-ring control of RUT-1 stage 1, in the same units'),
                reading='exploratory; a statement of the region, not a fit')


# ================================================================ driver
def gates(inputs, kubo, lenses):
    g = dict(G1_shell_kernel=gate_g1(), G2_derivatives=gate_g2(), G3_link_to_rut1=gate_g3(), G4_no_net_equivalent_mass=gate_g4(),
             G5_two_lensing_routes=gate_g5(), G6_lens_anchor=gate_g6(lenses), G7_xcop_ingestion=gate_g7(inputs),
             G8_coma_ingestion=gate_g8(kubo))
    for k, v in g.items():
        log(f"{k}: {'pass' if v['passed'] else 'FAIL'}" + (f", control {'rejected' if v['control']['rejected'] else 'NOT REJECTED'}" if 'control' in v else ''))
    return g


def main():
    args = evidence_io.parse(__doc__)
    inputs = json.loads(INPUTS.read_text(encoding='utf-8'))
    kubo = json.loads(KUBO.read_text(encoding='utf-8'))
    lenses = load_lenses()
    g = gates(inputs, kubo, lenses)
    E1, base, conv_c = e1(inputs)
    log(f"E1 best w {E1['best']['w']:.1f} kpc, Lambda {E1['best']['Lambda']:.4e}, chi2 {E1['best']['chi2']:.1f}")
    E2, conv_l = e2(lenses, E1['best'])
    E3 = e3(kubo, E1['best'])
    E4 = e4(inputs, base, E2['lenses'], lenses, E1['best'])
    g['G9_convergence'] = dict(cluster_scan_grid_doubling=conv_c, lens_lambda_grid_doubling=conv_l,
                               tolerances=dict(clusters=1e-5, lenses=1e-4), passed=bool(conv_c < 1e-5 and conv_l < 1e-4))
    log(f"G9: {'pass' if g['G9_convergence']['passed'] else 'FAIL'} ({conv_c:.2e}, {conv_l:.2e})")
    verified = all(v['passed'] for v in g.values()) and all(v['control']['rejected'] for v in g.values() if 'control' in v)
    j37 = E2['lenses'][0]
    result = dict(
        experiment='CL-1: the written-track response at cluster scale and for light',
        protocol='protocol-cl1.md',
        statuses=dict(numerical_verification='passed' if verified else 'FAILED',
                      scientific_outcome='exploratory readings on exposed data; see report-cl1.md'),
        gates=g, numerical_verification_passed=bool(verified),
        E1_clusters=E1, E2_lenses=E2, E3_coma=E3, E4_map=E4,
        input_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (HERE/'protocol-cl1.md', INPUTS, KUBO, HERE/'steady_field.py')},
        checks_short_run=dict(
            G1_worst=g['G1_shell_kernel']['worst_relative_difference'], G2_worst=max(g['G2_derivatives']['worst_gradient'], g['G2_derivatives']['worst_laplacian']),
            G3_w_over_ell=g['G3_link_to_rut1']['mature_ring_control']['w_over_ell'], G4_worst=g['G4_no_net_equivalent_mass']['worst'],
            G5_worst=g['G5_two_lensing_routes']['worst_relative_difference'], G6_worst_2d=g['G6_lens_anchor']['worst_2d'],
            G7_worst_mass=g['G7_xcop_ingestion']['worst_mass_deviation'],
            E1_best_w=E1['best']['w'], E1_best_Lambda=E1['best']['Lambda'], E1_best_chi2=E1['best']['chi2'],
            E1_newtonian_chi2=E1['references']['newtonian_baryons']['chi2'],
            E2_J0037_Chabrier_Lambda_at_10kpc=j37['Chabrier']['Lambda'][12], E2_width_index_10kpc=12,
            E3_low_bracket_L0_chi2=E3['brackets'][0]['L0']['shape_chi2']),
        runtime_seconds=round(time.time() - T0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    canonical = HERE/'cl1-results.json'
    if not canonical.exists():
        canonical.write_text(text, encoding='utf-8', newline='\n')
        log('first run: archive created')
    status = evidence_io.finish(args, 'path-memory-cl1', text, canonical, ignore={'/runtime_seconds'})
    print(json.dumps(dict(numerical_verification=result['statuses']['numerical_verification'], E1_best=E1['best'],
                          E1_references=E1['references'], E4_compatible_widths=E4['widths_where_compatible']), indent=1, default=float))
    return status if verified else 1


if __name__ == '__main__':
    raise SystemExit(main())

"""TF-1 stage 3: the six lenses under the well (W) and the steering law (S), two force families, through CL-2's
LensSystem (CR-2's measurement interface), static geometry PF1 primary."""
import sys
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.optimize import nnls, minimize_scalar
from common import PM, save, read
from support import support_factor

sys.path.insert(0, str(PM))
import cl2_sources as CS      # noqa: E402
import cl2_response as CR     # noqa: E402
import steady_field as SF     # noqa: E402
from scipy.integrate import quad   # noqa: E402
sys.path.insert(0, str(PM.parent/'radiation-polarized-gravity'))
import lensing as LZ          # noqa: E402  (RPG-1's deprojection, the route CL-1 reproduced to 1e-10)

G, C = SF.G, SF.C_KMS
WIDTHS = CR.WIDTHS[CR.WIDTHS >= CR.GALAXY_MIN_WIDTH]      # 21 widths, the same for every lens
SHELLS = np.geomspace(.2, 50., 16)
PRIMARY, SENSITIVITY = 'PF1_static_euclidean', 'G1_flat_FLRW'
EINSTEIN_WEIGHT = 1e3
CHI2_LIMIT, EINSTEIN_LIMIT = 128., .03


class Lens:
    def __init__(self, name, scenario):
        self.name, self.scenario = name, scenario
        self.L = L = CS.LensSystem(name, scenario, widths=WIDTHS)
        r = L.model.r
        b = L.bE
        self.shell_force = np.array([np.where(r >= rk, G*1e11/r**2, 0.) for rk in SHELLS])
        frac = np.array([1 - np.sqrt(max(0., 1 - (b/rk)**2)) if b < rk else 1. for rk in SHELLS])
        self.shell_bend = 4*G*1e11*frac/(C**2*b)
        self.families = dict(F1=(L.memforce, L.membend), F2=(self.shell_force, self.shell_bend))
        self.cache = {}
        D = np.diag(2*L.y)
        self.Lc = cholesky(D@L.cov@D, lower=True)
        self.n_bins = len(L.y)

    def coefficients(self, family, beta):
        key = (family, round(float(beta), 12))
        if key not in self.cache:
            self.L.model.forces = np.vstack([self.L.starforce, self.families[family][0]])
            self.cache[key] = self.L.model.coefficients(beta)
        return self.cache[key]

    def rows(self, family, coupling, imf, beta):
        m = self.L.pop[imf]/1e11
        co = self.coefficients(family, beta)
        star = 1. if coupling == 'W' else support_factor(beta)
        light = 1. if coupling == 'W' else .5
        A = (m*star*co[1:]).T
        y = self.L.y**2 - m*co[0]
        return dict(Ak=solve_triangular(self.Lc, A, lower=True), yk=solve_triangular(self.Lc, y, lower=True),
                    ae=m*light*self.families[family][1], ye=self.L.need - m*self.L.starbend, se=EINSTEIN_LIMIT*self.L.need,
                    m=m, co=co, star=star)

    def evaluate(self, family, coupling, imf, beta, amps):
        R = self.rows(family, coupling, imf, beta)
        v2 = R['m']*(R['co'][0] + (R['star']*np.asarray(amps))@R['co'][1:])
        return dict(chi2=self.L.chi2_v(v2), einstein_residual=float((R['ae']@amps - R['ye'])/self.L.need), beta=float(beta))

    def fit_beta(self, family, coupling, imf, amps):
        ev = lambda b: self.evaluate(family, coupling, imf, b, amps)['chi2']
        lo, hi = self.L.beta_bounds
        grid = np.linspace(lo, hi, 11)
        vals = [ev(b) for b in grid]
        i = int(np.argmin(vals))
        res = minimize_scalar(ev, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 10)]), method='bounded', options={'xatol': 1e-6})
        return (float(res.fun), float(res.x)) if res.fun < vals[i] else (float(vals[i]), float(grid[i]))

    def solve_at(self, family, coupling, imf, beta):
        R = self.rows(family, coupling, imf, beta)
        A = np.vstack([R['Ak'], (R['ae']/R['se']*EINSTEIN_WEIGHT)[None, :]])
        y = np.concatenate([R['yk'], [R['ye']/R['se']*EINSTEIN_WEIGHT]])
        amps, _ = nnls(A, y, maxiter=100000)
        grad = A.T@(A@amps - y)
        scale = float(np.max(np.abs(A.T@y))) or 1.
        act = amps > 0
        kkt = max(float(np.max(np.abs(grad[act]))/scale) if act.any() else 0., float(max(0., -np.min(grad[~act]))/scale) if (~act).any() else 0.)
        ev = self.evaluate(family, coupling, imf, beta, amps)
        return dict(**ev, amplitudes=amps.tolist(), active=int(act.sum()), kkt=kkt)

    def per_lens(self, family, coupling, imf):
        lo, hi = self.L.beta_bounds
        grid = np.linspace(lo, hi, 11)
        sols = [self.solve_at(family, coupling, imf, b) for b in grid]
        vals = [s['chi2'] for s in sols]
        i = int(np.argmin(vals))
        res = minimize_scalar(lambda b: self.solve_at(family, coupling, imf, b)['chi2'],
                              bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 10)]), method='bounded', options={'xatol': 1e-6})
        best = self.solve_at(family, coupling, imf, float(res.x)) if res.fun < vals[i] else sols[i]
        best['beta_on_bound'] = bool(abs(best['beta'] - lo) < 1e-6 or abs(best['beta'] - hi) < 1e-6)
        return best


def joint(lenses, family, coupling, imf, betas0, iters=6):
    betas = list(betas0)
    history, best = [], None
    for it in range(iters):
        R = [L.rows(family, coupling, imf, b) for L, b in zip(lenses, betas)]
        E = CR.Block('lens_einstein', np.array([r['ae'] for r in R]), np.array([r['ye'] for r in R]), np.array([r['se'] for r in R]))
        K = CR.Block('lens_kinematics', np.vstack([r['Ak'] for r in R]), np.concatenate([r['yk'] for r in R]), np.ones(sum(len(r['yk']) for r in R)))
        sol = CR.solve([E, K])
        amps = sol['amplitudes']
        ev = [L.evaluate(family, coupling, imf, b, amps) for L, b in zip(lenses, betas)]
        objective = float(sum(c/b.N for c, b in zip(sol['chi2'], (E, K))))
        history.append(dict(iteration=it, objective=objective, total_kinematics_chi2=float(sum(e['chi2'] for e in ev)),
                            max_abs_einstein_residual=float(max(abs(e['einstein_residual']) for e in ev)), betas=list(betas),
                            kkt=sol['kkt']))
        if best is not None and abs(objective - best['objective']) < 1e-6*max(1., abs(best['objective'])):
            best = dict(objective=objective, amplitudes=amps.tolist(), per_lens=ev, betas=list(betas), kkt=sol['kkt'], widths_active=sol['widths_active'])
            break
        best = dict(objective=objective, amplitudes=amps.tolist(), per_lens=ev, betas=list(betas), kkt=sol['kkt'], widths_active=sol['widths_active'])
        betas = [L.fit_beta(family, coupling, imf, amps)[1] for L in lenses]
    objs = [h['objective'] for h in history]
    total = float(sum(e['chi2'] for e in best['per_lens']))
    worst = float(max(abs(e['einstein_residual']) for e in best['per_lens']))
    return dict(**best, history=history, iterations=len(history), monotone=bool(all(np.diff(objs) <= 1e-9)),
                total_kinematics_chi2=total, max_abs_einstein_residual=worst,
                described=bool(total <= CHI2_LIMIT and worst <= EINSTEIN_LIMIT))


def reproduction(lenses_pf1, lenses_flrw):
    cl1 = {L['lens']: L for L in read(PM/'cl1-results.json')['E2_lenses']['lenses']}
    cl2 = read(PM/'cl2-results.json')['E2_lenses']['cr2_benchmarks']
    rows, worst = {}, 0.
    for Lp, Lf in zip(lenses_pf1, lenses_flrw):
        L = Lp.L
        mass = L.pop['Chabrier']
        comps, parts, dl = L.components, L.parts, L.geo['Dl_kpc']
        # CL-1's three-dimensional route: RPG-1's deprojected light profile and quadrature settings
        r3, frac3, _ = LZ.light_profile(comps, dl)
        lr = np.log(r3)
        gN = lambda x: G*mass*np.interp(np.log(x), lr, frac3, left=0., right=frac3[-1])/(x*x)
        bend = lambda b: 4/C**2*quad(lambda t: gN(b/np.cos(t))*b/np.cos(t), 0, np.pi/2, limit=200, epsabs=0, epsrel=1e-8)[0]
        theta = SF.einstein_radius(bend, dl, L.ratio)
        # the projected route (CL-1's two-dimensional check, 1e-5 there) and the component-model route, recorded only
        alpha = lambda b: SF.deflection_from_projected_mass(mass*SF.sersic_projected_fraction(comps, parts, dl, np.array([b]))[0], b)
        theta_projected = SF.einstein_radius(alpha, dl, L.ratio)
        theta_component = L.model.angle(mass, dl, L.ratio, False)
        ref = cl1[L.name]['Chabrier']
        so = Lf.L.stars_only()
        bench = cl2[L.name]['stars_only']
        r = dict(log10_mass=float(np.log10(mass)), cl1_log10_mass=ref['log10_stellar_mass'],
                 theta_newtonian=float(theta), cl1_theta_newtonian_L0=ref['theta_newtonian_L0'], theta_projected_route=float(theta_projected), theta_component_route=float(theta_component),
                 stars_only_flrw=so, cl2_stars_only=bench)
        r['errors'] = dict(mass=abs(r['log10_mass'] - r['cl1_log10_mass']), theta=abs(theta/ref['theta_newtonian_L0'] - 1),
                           chi2=abs(so['chi2']/bench['chi2'] - 1), beta=abs(so['beta'] - bench['beta']),
                           stellar_mass=abs(so['stellar_mass_Msun']/bench['stellar_mass_Msun'] - 1))
        worst = max(worst, max(r['errors'].values()))
        rows[L.name] = r
    return dict(rows=rows, worst=worst, routes_max_difference=float(max(abs(r['theta_component_route']/r['theta_newtonian'] - 1) for r in rows.values())))


def run(out):
    names = list(CS.LENSES)
    lenses = {s: [Lens(n, s) for n in names] for s in (PRIMARY, SENSITIVITY)}
    print('lenses built', flush=True)
    repro = reproduction(lenses[PRIMARY], lenses[SENSITIVITY])
    print('reproduction worst', repro['worst'], flush=True)
    betas0 = [L.L.stars_only()['beta'] for L in lenses[PRIMARY]]
    results = dict(per_lens={}, joint={})
    for imf in ('Chabrier', 'Salpeter'):
        for family in ('F1', 'F2'):
            for coupling in ('W', 'S'):
                key = f'{imf}-{family}-{coupling}'
                per = {L.name: L.per_lens(family, coupling, imf) for L in lenses[PRIMARY]}
                results['per_lens'][key] = dict(lenses=per, total_kinematics_chi2=float(sum(p['chi2'] for p in per.values())),
                                                max_abs_einstein_residual=float(max(abs(p['einstein_residual']) for p in per.values())),
                                                betas_on_bound=int(sum(p['beta_on_bound'] for p in per.values())),
                                                worst_kkt=float(max(p['kkt'] for p in per.values())))
                print('per-lens', key, results['per_lens'][key]['total_kinematics_chi2'], flush=True)
        for coupling in ('W', 'S'):
            key = f'{imf}-F1-{coupling}'
            results['joint'][key] = joint(lenses[PRIMARY], 'F1', coupling, imf, betas0)
            print('joint', key, results['joint'][key]['total_kinematics_chi2'], results['joint'][key]['max_abs_einstein_residual'], flush=True)
    betas_f = [L.L.stars_only()['beta'] for L in lenses[SENSITIVITY]]
    results['sensitivity_flrw'] = {c: joint(lenses[SENSITIVITY], 'F1', c, 'Chabrier', betas_f) for c in ('W', 'S')}
    results['support_factors'] = {str(b): support_factor(b) for b in (-2., -1., 0., .45)}
    results['vector_channel'] = dict(fitted=False, reason='G7: zero mean support for a non-rotating population; G8: no focusing',
                                     light_to_matter_ratio={f'{v:.0f}_km_s': C/v for v in (200., 300.)})
    save(out/'lens-details.json', results)
    jw, js = results['joint']['Chabrier-F1-W'], results['joint']['Chabrier-F1-S']
    cw, cs = results['per_lens']['Chabrier-F2-W']['total_kinematics_chi2'], results['per_lens']['Chabrier-F2-S']['total_kinematics_chi2']
    if cw > 0 and cs > 0:
        ratio = cs/cw
        ceiling = 'disfavoured' if ratio > 2 else ('favoured' if ratio < .5 else 'indistinguishable')
    else:
        ratio, ceiling = None, 'indistinguishable (both ceilings at zero)'
    if js['described'] and not jw['described']:
        outcome = 'a: the steering coupling describes the lenses under the universal family where the well does not'
    elif js['described'] == jw['described']:
        outcome = 'b: the lens verdict is unchanged by the transverse coupling'
    else:
        outcome = 'b-: the well describes the lenses and the steering coupling does not'
    gates = dict(reproduction=repro['worst'] < 1e-6, joint_monotone=all(j['monotone'] for j in results['joint'].values()),
                 kkt=all(v['worst_kkt'] < 1e-8 for v in results['per_lens'].values()) and all(max(j['kkt']['max_abs_gradient_active'], max(0., -j['kkt']['min_gradient_inactive'])) < 1e-8 for j in results['joint'].values()))
    summary = dict(stage='lenses', gates=gates, numerical_pass=all(gates.values()), reproduction=repro, outcome=outcome,
                   joint={k: {kk: v[kk] for kk in ('total_kinematics_chi2', 'max_abs_einstein_residual', 'described', 'iterations', 'monotone', 'betas', 'widths_active', 'objective')} for k, v in results['joint'].items()},
                   per_lens={k: {kk: v[kk] for kk in ('total_kinematics_chi2', 'max_abs_einstein_residual', 'betas_on_bound', 'worst_kkt')} for k, v in results['per_lens'].items()},
                   per_lens_chi2={k: {n: dict(chi2=p['chi2'], beta=p['beta'], einstein_residual=p['einstein_residual'], active=p['active']) for n, p in v['lenses'].items()} for k, v in results['per_lens'].items()},
                   ceiling_ratio_S_over_W=ratio, ceiling_reading=ceiling,
                   sensitivity_flrw={k: dict(total_kinematics_chi2=v['total_kinematics_chi2'], max_abs_einstein_residual=v['max_abs_einstein_residual'], described=v['described']) for k, v in results['sensitivity_flrw'].items()},
                   support_factors=results['support_factors'], vector_channel=results['vector_channel'],
                   acceptance=dict(chi2_limit=CHI2_LIMIT, einstein_limit=EINSTEIN_LIMIT, bins=int(sum(L.n_bins for L in lenses[PRIMARY]))),
                   geometry=PRIMARY, widths_kpc=WIDTHS.tolist(), shells_kpc=SHELLS.tolist())
    return summary

"""CR-2: one repulsive support law against joint lensing and stellar motions (see protocol.md).

    python cr2.py [--canonical] [--output-dir DIR]        (CR2_SMOKE=1 runs one lens on a short grid)

A Thomas–Fermi condensate (P = K rho^2/2 with one shared K, quoted as R_TF) in equilibrium with each lens's stars
is fitted jointly to the KCWI stellar motions and the exact lens constraint. It is fitted in one geometry at a
time, with the stellar mass either lens-determined (M1) or fixed to population values (M2). It uses 6 workers,
fails fast and logs its progress.
"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import brentq, minimize, minimize_scalar

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
for _p in (RESULTS/'capture-to-orbit', RESULTS/'companion-extensions', HERE, RESULTS/'slacs-component-refit'):
    sys.path.insert(0, str(_p))
import evidence_io  # noqa: E402
import tf  # noqa: E402
from model import ComponentModel, G, C, ARCSEC  # noqa: E402  (slacs-component-refit)

ALPHA = 0.0002488993286382367                          # per Mpc, the repository's redshift rate
MPC_M, KPC_M, C_MS = 3.0856775814913673e22, 3.0856775814913673e19, 299792458.
LSUN_W, MSUN_KG, GYR_S, EV_J, HBAR = 3.828e26, 1.98847e30, 3.15576e16, 1.602176634e-19, 1.054571817e-34
CONSTITUENTS_EV = (1.34e-24, 1e-22)                     # CR-1's two masses, for the healing-length check
GEOMETRIES = ('G1_flat_FLRW', 'G2_coscaling_coasting')
SMOKE = os.environ.get('CR2_SMOKE') == '1'
R_GRID = np.geomspace(.3, 300, 5 if SMOKE else 13)      # kpc
BUDGET_S = 90*60
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:7.1f}s] {msg}', flush=True)


def read(folder, file='results.json'):
    return json.loads((RESULTS/folder/file).read_text(encoding='utf-8'))


DATA = read('slacs-resolved-input-audit')['systems']
PROFILES = {r['Name']: r for r in read('slacs-light-profile-audit')['rows']}
PILOT = {r['Name']: r for r in read('slacs-motion-lensing-pilot')['rows'] if r['model'] == 'empirical_extra' and r['cutoff_in_a'] == 20}
GEO = {r['Name']: r for r in read('lensing-data-readiness', 'conditional-geometry.json')}
OBS = {r['Name']: r for r in read('lensing-data-readiness', 'lens-observations-and-image-models.json')}
CFG = read('slacs-outer-bin-check', 'protocol.json')
_allowed = {r['Name'] for r in read('slacs-outer-bin-check')['rows'] if r['model'] == 'empirical_extra'}
LENSES = [d['Name'] for d in DATA if d['Name'] in _allowed][:1 if SMOKE else None]
PUB = {}
for _r in read('lens-photometric-audit', 'normalization-sensitivity.json'):
    PUB.setdefault((_r['Name'], _r['imf']), _r['published_log10_stellar_mass'])
ARCHIVE = [r for r in json.loads((RESULTS/'isotropic-galaxy-transfer/nfw-geometry-results.json').read_text())['rows']
           if r['geometry'] == 'standard_flat_FLRW']


def chi_flrw(z):
    """Comoving distance in Mpc, flat FLRW H0 = 70, Omega_m = 0.3 (as nfw-geometry.py)."""
    return C/70*quad(lambda t: 1/np.sqrt(.3*(1 + t)**3 + .7), 0, z, epsabs=1e-9)[0]


def geometry(name, which):
    zl, zs = OBS[name]['zFG'], OBS[name]['zBG']
    DM = chi_flrw if which == 'G1_flat_FLRW' else (lambda z: math.log1p(z)/ALPHA)
    return dict(zl=zl, zs=zs, Dl_kpc=DM(zl)/(1 + zl)*1000, ratio=1 - DM(zl)/DM(zs), DL_Mpc=DM(zl)*(1 + zl),
                DL_FLRW_Mpc=chi_flrw(zl)*(1 + zl))


def pmass(x):
    """Projected NFW mass inside x = R/rs, in units of the halo normalization (as nfw-geometry.py)."""
    if x < 1:
        return np.log(x/2) + np.arccosh(1/x)/np.sqrt(1 - x*x)
    if x > 1:
        return np.log(x/2) + np.arccos(1/x)/np.sqrt(x*x - 1)
    return 1 - np.log(2)


class Lens:
    def __init__(self, name, which):
        self.name, self.which = name, which
        g = self.g = geometry(name, which)
        dl = g['Dl_kpc']
        item = next(d for d in DATA if d['Name'] == name)
        a = PILOT[name]['scale_a_kpc']*dl/(GEO[name]['conditional_Dl_Mpc']*1000)
        edges = np.r_[item['inner_arcsec'], item['outer_arcsec'][-1]]*dl/ARCSEC
        psf = item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
        comps = [dict(R=q['R_arcsec']*dl/ARCSEC, n=q['n'], amp=q['amp_at_R'], bn=q['bn']) for q in PROFILES[name]['components']]
        self.model = model = ComponentModel(a, edges, psf, 0, .5, 1, 20, comps)
        r = model.r
        self.starforce = model.forces[0].copy()                      # G 1e11 frac/r^2
        self.frac = self.starforce*r*r/(G*1e11)
        self.F1 = self.frac/r**2
        self.Re = PROFILES[name]['computed_equal_area_half_light_arcsec']*dl/ARCSEC
        self.theta = PILOT[name]['catalog_SIE_arcsec']/ARCSEC
        self.bE, self.ratio, self.need = self.theta*dl, g['ratio'], self.theta/g['ratio']
        self.S = quad(lambda t: model.mass_fraction(self.bE/np.cos(t))/(self.bE/np.cos(t)), 0, np.pi/2, epsabs=1e-9, epsrel=1e-8)[0]
        self.starbend = 4*G*1e11/C**2*self.S                            # radians per 1e11 Msun, as nfw-geometry.py
        self.y = np.array(item['vrms_kms'])
        self.fac = cho_factor(np.array(item['covariance_kms_squared']))
        self.mb = np.array(CFG['mass_Msun_bounds'])
        self.pop = {imf: 10**PUB[name, imf]*(g['DL_Mpc']/g['DL_FLRW_Mpc'])**2 for imf in ('Chabrier', 'Salpeter')}

    def chi2(self, s2):
        if np.any(s2 <= 0):
            return 1e30
        e = np.sqrt(s2) - self.y
        return float(e @ cho_solve(self.fac, e))

    def fit_beta(self, forces, weights):
        """Constant anisotropy: a coarse grid, then bounded Brent around the best node."""
        self.model.forces = np.asarray(forces)
        w = np.asarray(weights, float)
        ev = lambda beta: self.chi2(w @ self.model.coefficients(beta))
        lo, hi = CFG['constant_beta_bounds']
        grid = np.linspace(lo, hi, 11)
        vals = [ev(b) for b in grid]
        i = int(np.argmin(vals))
        res = minimize_scalar(ev, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 10)]), method='bounded', options={'xatol': 1e-7})
        return (float(res.fun), float(res.x)) if res.fun < vals[i] else (float(vals[i]), float(grid[i]))

    # Benchmarks, replicating nfw-geometry.py exactly (same starts, bounds and optimizer settings).
    def stars_only(self):
        mstar = self.need/self.starbend
        if not self.mb[0]/1e11 <= mstar <= self.mb[1]/1e11:
            return None
        self.model.forces = np.array([self.starforce, np.zeros_like(self.starforce)])
        ev = lambda beta: self.chi2(self.model.coefficients(beta)[0]*mstar)
        opts = [minimize(lambda p: ev(p[0]), [b], bounds=[CFG['constant_beta_bounds']], method='L-BFGS-B') for b in CFG['starts_beta']]
        best = min((o for o in opts if o.success), key=lambda o: o.fun)
        return dict(chi2=float(best.fun), beta=float(best.x[0]), stellar_mass_Msun=float(mstar*1e11))

    def nfw(self):
        r = self.model.r
        flo, fhi = max(0., 1 - self.mb[1]/1e11*self.starbend/self.need), min(1., 1 - self.mb[0]/1e11*self.starbend/self.need)

        def ev(p):
            f, logscale, beta = p
            rs = self.Re*10**logscale
            x = r/rs
            halo = 1e11*(np.log1p(x) - x/(1 + x))
            bend = 4*G*1e11/C**2/self.bE*pmass(self.bE/rs)
            self.model.forces = np.array([self.starforce, G*halo/r**2])
            cb, ch = self.model.coefficients(beta)
            s2 = cb*(1 - f)*self.need/self.starbend + ch*f*self.need/bend
            e = np.sqrt(np.maximum(s2, 1e-100)) - self.y
            return float(e @ cho_solve(self.fac, e))
        starts = [[flo + (fhi - flo)*f, s, b] for f in (.1, .7) for s in (-1., 0., 1.) for b in (-.5, .3)]
        opts = [minimize(ev, p, bounds=[(flo, fhi), (-2, 2), CFG['constant_beta_bounds']], method='L-BFGS-B',
                         options={'ftol': 1e-11, 'maxiter': 400}) for p in starts]
        best = min((o for o in opts if o.success and np.isfinite(o.fun)), key=lambda o: o.fun)
        return dict(chi2=float(best.fun), share=float(best.x[0]), rs_over_Re=float(10**best.x[1]), beta=float(best.x[2]))

    # The declared support law.
    def terms(self, fam, lam):
        """Condensate solution, its deflection integral A = int m(k b sec t)/(b sec t) dt and its unit force."""
        s = fam.solve(lam)
        k, b, m = fam.k, self.bE, s['m']
        re = s['x_edge']/k
        nodes, weights = np.polynomial.legendre.leggauss(64)

        def panel(t0, t1):
            t = .5*(t1 - t0)*nodes + .5*(t1 + t0)
            rr = b/np.cos(t)
            return .5*(t1 - t0)*float(np.sum(weights*m(k*rr)/rr))
        A = panel(0., math.acos(b/re)) + panel(math.acos(b/re), math.pi/2) if re > b else panel(0., math.pi/2)
        return s, A, m(k*self.model.r)/self.model.r**2

    def m1_at(self, fam, lam, beta_fit=True):
        s, A, F2 = self.terms(fam, lam)
        k = fam.k
        psi = self.need*C**2*k/(4*(lam*self.S + A))                     # exact lens: stars plus condensate
        row = dict(lam=lam, share=A/(lam*self.S + A), stellar_mass_Msun=lam*psi/(G*k), condensate_mass_Msun=psi*s['m_edge']/(G*k),
                   psi_c=psi, r_edge_kpc=s['x_edge']/k)
        row['admissible'] = bool(self.mb[0] <= row['stellar_mass_Msun'] <= self.mb[1])
        if beta_fit and row['admissible']:
            row['chi2'], row['beta'] = self.fit_beta([self.F1, F2], [psi*lam/k, psi/k])
        return row

    def m1_fit(self, R_tf):
        fam = tf.Family(self.model.r, self.frac, tf.k_of(R_tf))
        lams = np.geomspace(1e-5, 1e7, 25)
        shares = np.array([self.m1_at(fam, l, False)['share'] for l in lams])
        order = np.argsort(shares)
        targets = np.r_[.002, .01, .03, np.linspace(.06, .94, 12), .97, .99, .998]
        trial = [math.exp(np.interp(t, shares[order], np.log(lams[order]))) for t in targets if shares.min() < t < shares.max()]
        rows = [self.m1_at(fam, l) for l in trial]
        rows = [r for r in rows if r['admissible']]
        if not rows:
            return dict(R_tf_kpc=R_tf, admissible=False)
        i = min(range(len(rows)), key=lambda j: rows[j]['chi2'])
        lo = math.log(rows[i + 1]['lam']) if i + 1 < len(rows) else math.log(rows[i]['lam']) - .5
        hi = math.log(rows[i - 1]['lam']) if i > 0 else math.log(rows[i]['lam']) + .5
        lo, hi = min(lo, hi), max(lo, hi)

        def ev(ll):
            q = self.m1_at(fam, math.exp(ll))
            return q['chi2'] if q['admissible'] else 1e30
        res = minimize_scalar(ev, bounds=(lo, hi), method='bounded', options={'xatol': 1e-4})
        best = self.m1_at(fam, math.exp(res.x))
        if not best.get('admissible') or best['chi2'] > rows[i]['chi2']:
            best = rows[i]
        best['R_tf_kpc'], best['samples'] = R_tf, len(rows)
        return best

    def m2_fit(self, R_tf, imf):
        mpop = self.pop[imf]
        rem = self.need - 4*G*mpop*self.S/C**2
        out = dict(R_tf_kpc=R_tf, imf=imf, population_mass_Msun=mpop, star_share_of_lens=1 - rem/self.need)
        if rem <= 0:
            out['admissible'] = False
            return out
        fam = tf.Family(self.model.r, self.frac, tf.k_of(R_tf))
        target = rem*C**2/(4*G*mpop)                                     # A(lam)/lam required
        g = lambda ll: math.log(self.terms(fam, math.exp(ll))[1]/math.exp(ll)/target)
        ll = brentq(g, math.log(1e-8), math.log(1e8), xtol=1e-10)
        lam = math.exp(ll)
        s, A, F2 = self.terms(fam, lam)
        psi = G*fam.k*mpop/lam
        out['chi2'], out['beta'] = self.fit_beta([self.F1, F2], [psi*lam/fam.k, psi/fam.k])
        out.update(admissible=True, lam=lam, psi_c=psi, condensate_mass_Msun=psi*s['m_edge']/(G*fam.k), r_edge_kpc=s['x_edge']/fam.k)
        return out

    def closure(self, R_tf, lam):
        """V3: independent quadrature of the total deflection at the catalogue radius, split at the edge."""
        fam = tf.Family(self.model.r, self.frac, tf.k_of(R_tf))
        q = self.m1_at(fam, lam, False)
        s = fam.solve(lam)
        k, b, re = fam.k, self.bE, s['x_edge']/fam.k

        def integrand(t):
            rr = b/math.cos(t)
            return (q['stellar_mass_Msun']*self.model.mass_fraction(rr) + q['psi_c']*float(s['m'](np.array([k*rr]))[0])/(G*k))/rr
        alpha = 4*G/C**2*quad(integrand, 0, math.pi/2, points=[math.acos(b/re)] if re > b else None, epsabs=0, epsrel=1e-10, limit=400)[0]
        return float(alpha*self.ratio/self.theta - 1)


_CACHE = {}


def _lens(name, which):
    if (name, which) not in _CACHE:
        _CACHE[name, which] = Lens(name, which)
    return _CACHE[name, which]


def _init():
    np.seterr(over='raise', invalid='raise', divide='raise')


def task(job):
    kind, name, which = job[:3]
    L = _lens(name, which)
    if kind == 'bench':
        return dict(stars_only=L.stars_only(), nfw=L.nfw())
    if kind == 'm1':
        return L.m1_fit(job[3])
    if kind == 'm2':
        return L.m2_fit(job[3], job[4])
    if kind == 'final':
        R_tf, lam = job[3], job[4]
        return dict(closure=L.closure(R_tf, lam), pop=L.pop, Re_kpc=L.Re, bE_kpc=L.bE, Dl_kpc=L.g['Dl_kpc'], ratio=L.ratio)
    raise ValueError(kind)


def milky_way(R_tf):
    """The same R_TF in the Milky Way: spherical baryon model I for the equilibrium, the archive's baseline-I speeds."""
    import inputs as I
    base = next(r for r in I.milky_way_runs() if r['baryons'] == 'I' and abs(r['rd'] - 2.6) < 1e-9 and abs(r['lf'] - 1) < 1e-9)
    M = I.milky_way_receivers('I')['mass']
    fam = tf.Family(I.GRID, M/M[-1], tf.k_of(R_tf))
    R, y, vb = base['R'], base['y'], base['vb']
    rmse = lambda v, n=None: float(np.sqrt(np.mean((v[:n] - y[:n])**2)))

    def speeds(ll):
        s = fam.solve(math.exp(ll))
        return np.sqrt(vb**2 + G*M[-1]*s['m'](fam.k*R)/math.exp(ll)/R), s
    grid = np.linspace(math.log(1e-4), math.log(1e4), 33)
    vals = [rmse(speeds(v)[0]) for v in grid]
    i = int(np.argmin(vals))
    res = minimize_scalar(lambda v: rmse(speeds(v)[0]), bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 32)]), method='bounded', options={'xatol': 1e-5})
    v, s = speeds(res.x)
    lam = math.exp(res.x)
    return dict(R_tf_kpc=R_tf, lam=lam, condensate_mass_Msun=M[-1]*s['m_edge']/lam, edge_kpc=s['x_edge']/fam.k,
                rmse_38=rmse(v), rmse_inner20=rmse(v, 20), baryons_rmse_38=rmse(vb), baryons_rmse_inner20=rmse(vb, 20),
                predicted_kms=v.tolist(), observed_kms=y.tolist(), R_kpc=R.tolist())


def main():
    args = evidence_io.parse(__doc__)
    _init()
    fam = tf.Family(np.geomspace(1e-5, 1e5, 101), np.geomspace(1e-5, 1e5, 101)**3/(1 + np.geomspace(1e-5, 1e5, 101)**3), 1.)
    s0 = fam.solve(0.)
    v1 = dict(edge_minus_pi=s0['x_edge'] - math.pi, mass_minus_pi=s0['m_edge'] - math.pi)
    v1['passed'] = bool(abs(v1['edge_minus_pi']) < 1e-8 and abs(v1['mass_minus_pi']) < 1e-8)
    log(f'V1 {v1}')
    out = dict(validation=dict(V1=v1), geometries={})
    with Pool(min(6, len(LENSES)), initializer=_init) as pool:
        for which in GEOMETRIES:
            log(f'{which}: benchmarks')
            bench = dict(zip(LENSES, pool.map(task, [('bench', n, which) for n in LENSES])))
            g = out['geometries'][which] = dict(benchmarks=bench,
                                                 benchmark_totals=dict(stars_only=sum(b['stars_only']['chi2'] for b in bench.values() if b['stars_only']),
                                                                       nfw=sum(b['nfw']['chi2'] for b in bench.values())))
            if which == 'G1_flat_FLRW':
                diffs = []
                for n, b in bench.items():
                    for kind, key in (('stars_only_exact_lens', 'stars_only'), ('NFW_exact_lens', 'nfw')):
                        arc = next(r for r in ARCHIVE if r['Name'] == n and r['model'] == kind)['stellar_chi2']
                        diffs.append(abs(b[key]['chi2']/arc - 1))
                out['validation']['V2'] = dict(max_relative_difference=float(max(diffs)), passed=bool(max(diffs) < 1e-4))
                log(f"V2 {out['validation']['V2']}")
            log(f"{which}: benchmark totals {g['benchmark_totals']}")
            profile = []
            for R in R_GRID:
                rows = pool.map(task, [('m1', n, which, float(R)) for n in LENSES])
                total = sum(r['chi2'] for r in rows) if all(r.get('admissible') for r in rows) else float('inf')
                profile.append(dict(R_tf_kpc=float(R), total_chi2=total, rows=rows))
                log(f'{which} M1 R_TF={R:.3g} kpc total chi2={total:.4g}')
                if time.time() - T0 > BUDGET_S:
                    raise RuntimeError('wall-clock budget exhausted')
            i = int(np.argmin([p['total_chi2'] for p in profile]))
            lo, hi = math.log(R_GRID[max(i - 1, 0)]), math.log(R_GRID[min(i + 1, len(R_GRID) - 1)])
            memo = {}

            def total(lr):
                rows = pool.map(task, [('m1', n, which, math.exp(lr)) for n in LENSES])
                memo[lr] = rows
                t = sum(r['chi2'] for r in rows) if all(r.get('admissible') for r in rows) else 1e30
                log(f'{which} M1 refine R_TF={math.exp(lr):.4g} kpc total chi2={t:.6g}')
                return t
            res = minimize_scalar(total, bounds=(lo, hi), method='bounded', options={'xatol': .01})
            best_rows = memo.get(res.x) or pool.map(task, [('m1', n, which, math.exp(res.x)) for n in LENSES])
            best_total = sum(r['chi2'] for r in best_rows)
            if best_total > profile[i]['total_chi2']:
                best_rows, best_R, best_total = profile[i]['rows'], float(R_GRID[i]), profile[i]['total_chi2']
            else:
                best_R = math.exp(res.x)
            finals = pool.map(task, [('final', n, which, best_R, r['lam']) for n, r in zip(LENSES, best_rows)])
            for n, r, f in zip(LENSES, best_rows, finals):
                r.update(name=n, lens_closure=f['closure'], Re_kpc=f['Re_kpc'], bE_kpc=f['bE_kpc'],
                         stellar_over_population=dict((imf, r['stellar_mass_Msun']/f['pop'][imf]) for imf in f['pop']),
                         supply_upper_bound_Msun=ALPHA/MPC_M*f['pop']['Chabrier']*LSUN_W*r['r_edge_kpc']*KPC_M*10*GYR_S/C_MS**2/MSUN_KG,
                         healing_length_over_edge={f'{m:g}': HBAR*C_MS**2/(m*EV_J)/(KPC_M*1e3)/math.sqrt(2*r['psi_c'])/r['r_edge_kpc']
                                                   for m in CONSTITUENTS_EV})
                r['needed_over_supply'] = r['condensate_mass_Msun']/r['supply_upper_bound_Msun']
            g['M1'] = dict(profile=[dict(R_tf_kpc=p['R_tf_kpc'], total_chi2=p['total_chi2']) for p in profile],
                           best_R_tf_kpc=best_R, best_total_chi2=best_total, rows=best_rows)
            log(f'{which} M1 best R_TF={best_R:.4g} kpc total chi2={best_total:.6g}; NFW {g["benchmark_totals"]["nfw"]:.6g}')
            g['M2'] = {}
            for imf in ('Chabrier', 'Salpeter'):
                prof = []
                for R in R_GRID:
                    rows = pool.map(task, [('m2', n, which, float(R), imf) for n in LENSES])
                    ok = [r for r in rows if r['admissible']]
                    prof.append(dict(R_tf_kpc=float(R), admissible=len(ok), total_chi2_admissible=sum(r['chi2'] for r in ok), rows=rows))
                j = int(np.argmin([p['total_chi2_admissible'] if p['admissible'] == len(LENSES) else 1e30 + p['total_chi2_admissible'] for p in prof]))
                g['M2'][imf] = dict(profile=[{k: v for k, v in p.items() if k != 'rows'} for p in prof], best=prof[j])
                log(f"{which} M2 {imf}: best R_TF={prof[j]['R_tf_kpc']:.3g} admissible {prof[j]['admissible']}/{len(LENSES)} chi2={prof[j]['total_chi2_admissible']:.5g}")
            g['milky_way'] = milky_way(best_R)
            mw = g['milky_way']
            log(f"{which} Milky Way at R_TF={best_R:.4g}: RMSE {mw['rmse_38']:.3f} (baryons {mw['baryons_rmse_38']:.3f}), inner {mw['rmse_inner20']:.3f} (baryons {mw['baryons_rmse_inner20']:.3f})")
            admissible = all(r.get('admissible') and 0 <= r['share'] < 1 for r in best_rows)
            g['verdict'] = dict(admissible=bool(admissible),
                                chi2_within_nfw_plus_10=bool(best_total <= g['benchmark_totals']['nfw'] + 10),
                                milky_way=bool(mw['rmse_38'] <= 20 and mw['rmse_inner20'] <= mw['baryons_rmse_inner20']))
            g['verdict']['passed'] = all(g['verdict'].values())
    out['validation']['V3'] = dict(max_abs_closure=float(max(abs(r['lens_closure']) for g in out['geometries'].values() for r in g['M1']['rows'])))
    out['validation']['V3']['passed'] = out['validation']['V3']['max_abs_closure'] < 1e-8
    out['all_validation_passed'] = all(v['passed'] for v in out['validation'].values())
    out['runtime_seconds'] = time.time() - T0
    text = json.dumps(out, indent=1, default=float)
    log(f"validation {'passed' if out['all_validation_passed'] else 'FAILED'}; verdicts " +
        str({w: g['verdict'] for w, g in out['geometries'].items()}))
    if SMOKE:
        print(text[:3000])
        return 0
    status = evidence_io.finish(args, 'supported-reservoir-cr2', text, HERE/'cr2-results.json', ignore={'/runtime_seconds'})
    return status if out['all_validation_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

"""CL-2's observation blocks with their anchors (protocol-cl2.md, section 1): the X-COP forward-pressure block,
the SPARC and Milky Way galaxy blocks, and the SLACS lens block through CR-2's measurement interface.
Everything is per unit Lambda per width; nothing here chooses a law."""
import json
import math
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
for _p in ('capture-to-orbit', 'companion-extensions', 'radiation-polarized-gravity', 'slacs-component-refit', 'supported-reservoir'):
    sys.path.insert(0, str(RESULTS/_p))
import steady_field as SF  # noqa: E402
import cl1 as C1  # noqa: E402
import cl2_response as CR  # noqa: E402
import inputs as I  # noqa: E402
import baryons as BAR  # noqa: E402
from model import ComponentModel, ARCSEC  # noqa: E402  (slacs-component-refit)

G = SF.G
C_KMS = SF.C_KMS
MU, MU_E = .6, 1.14                                   # the release's molecular weights
MP_KG, KEV_CM3_PA, MSUN_KG = 1.67262192e-27, 1.602176634e-10, 1.98847e30
A_STAR_REPO, A_STAR_PM1 = 8.563335e-11*SF.SI_TO_CODE, 6.54e-11*SF.SI_TO_CODE
ALPHA0 = 0.0002488993286382367                        # per Mpc, the repository's redshift rate
XCOP = HERE/'cl2-inputs-xcop-profiles.json'
# Eckert et al. 2019 Table 2: non-thermal pressure fractions at R500 and R200 in per cent; upper limits taken as 0
ALPHA_NT = dict(A1644=(0., 0.), A1795=(2.2, 6.7), A2029=(6.0, 10.4), A2142=(15.8, 18.6), A2255=(5.6, 6.1), A2319=(43.6, 52.3),
                A3158=(8.5, 12.5), A3266=(0., 0.), A644=(3.2, 5.6), A85=(10.2, 11.5), RXC1825=(5.1, 15.2), ZW1215=(0., 0.))
LENSES = ('J0037-0942', 'J1112+0826', 'J1204+0358', 'J1402+6321', 'J1621+3931', 'J1630+4520')
HELD_OUT_LENS = 'J1538+5817'


def _loglog_interp(x, xp, fp, left_power=None, right='hold'):
    """Interpolate log f against log x; inside the first node scale as x^left_power (or hold), beyond hold."""
    x, xp, fp = np.asarray(x, float), np.asarray(xp, float), np.asarray(fp, float)
    out = np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))
    if left_power is not None:
        inner = x < xp[0]
        out[inner] = fp[0]*(x[inner]/xp[0])**left_power
    return out


# ================================================================ clusters: the X-COP forward pressure
def load_xcop():
    return json.loads(XCOP.read_text(encoding='utf-8'))


def stellar_fraction_profile(extract):
    """Median over the clusters with measured stars of M_star(<r)/M_gas(<r) against r/R500, for the others."""
    xs = np.geomspace(.01, 2., 60)
    ratios = []
    for name, c in extract['clusters'].items():
        if 'stellar_mass' not in c:
            continue
        R500 = c['header']['R500_kpc']
        sm = c['stellar_mass']
        gm = c['gas_mass']
        ms = _loglog_interp(xs*R500, sm['radius_kpc'], sm['Mstar'], left_power=3.)
        mg = _loglog_interp(xs, gm['RADIUS'], gm['MGAS'], left_power=3.)
        ratios.append(ms/mg)
    return xs, np.median(np.array(ratios), axis=0), len(ratios)


def build_cluster(name, c, widths, n_grid=1500, mu=MU, mu_e=MU_E, alpha_nt=False, stars='release', frac_profile=None,
                  truncate_gas_at=None):
    """One cluster's pressure block: grid, baryons, pressure points and the per-width pressure operators."""
    R500 = c['header']['R500_kpc']
    d = c['density']
    rm = np.sqrt(np.array(d['r_in_kpc'])*np.array(d['r_out_kpc']))
    ne = np.array(d['ne_cm3'])
    r_out = min(d['r_out_kpc'][-1], c['pressure']['sz']['r_over_R500'][-1]*R500)
    r = np.geomspace(1., r_out, n_grid)
    lne = np.interp(np.log(r), np.log(rm), np.log(ne), left=np.log(ne[0]))
    slope = (np.log(ne[-1]) - np.log(ne[-2]))/(np.log(rm[-1]) - np.log(rm[-2]))
    beyond = r > rm[-1]
    lne[beyond] = np.log(ne[-1]) + slope*(np.log(r[beyond]) - np.log(rm[-1]))
    n_e = np.exp(lne)
    rho_gas = mu_e*MP_KG*n_e*1e6*SF.KPC_M**3/MSUN_KG                     # Msun/kpc^3
    f = 4*np.pi*r*r*rho_gas
    M_int = np.concatenate([[0.], np.cumsum(.5*(f[1:] + f[:-1])*np.diff(r))]) + 4/3*np.pi*r[0]**3*rho_gas[0]
    gm = c['gas_mass']
    r_rel, M_rel = np.array(gm['RADIUS'])*R500, np.array(gm['MGAS'])
    M_gas = _loglog_interp(r, r_rel, M_rel)
    inner = r < r_rel[0]
    M_gas[inner] = M_int[inner]*(M_rel[0]/np.interp(r_rel[0], r, M_int))     # the release's value carried inward with the measured shape
    if truncate_gas_at is not None:
        M_gas = np.minimum(M_gas, np.interp(truncate_gas_at*R500, r, M_gas))
    if stars == 'release' and 'stellar_mass' in c:
        sm = c['stellar_mass']
        M_star = _loglog_interp(r, sm['radius_kpc'], sm['Mstar'], left_power=3.)
        star_source = 'release'
    elif stars == 'none':
        M_star = np.zeros_like(r)
        star_source = 'none'
    else:
        xs, frac, n = frac_profile
        M_star = np.interp(np.log(r/R500), np.log(xs), frac, left=frac[0], right=frac[-1])*M_gas
        star_source = f'median fraction profile of {n} clusters'
    M_b = M_gas + M_star
    P500 = c['pressure']['P500_keV_cm3']
    pts = []
    for key in ('xray', 'sz'):
        t = c['pressure'][key]
        for x, p, e in zip(t['r_over_R500'], t['P_over_P500'], t['err']):
            if r[0] <= x*R500 <= r_out:
                pts.append((x*R500, p*P500, e*P500, key))
    pts.sort()
    rp = np.array([p[0] for p in pts])
    Pobs = np.array([p[1] for p in pts])
    eP = np.array([p[2] for p in pts])
    kinds = [p[3] for p in pts]
    conv = mu*MP_KG*1e6*(1e6/SF.KPC_M)*SF.KPC_M/KEV_CM3_PA           # n_e cm^-3, g code units, ds kpc -> keV cm^-3
    thermal = np.ones_like(rp)
    if alpha_nt:
        a500, a200 = ALPHA_NT[name]
        x = rp/R500
        thermal = 1 - np.minimum(a500*x, a200)/100.        # linear from zero at the centre, capped at the R200 value

    def pressure_of(g):
        integrand = n_e*g*conv
        cum = np.concatenate([[0.], np.cumsum(.5*(integrand[1:] + integrand[:-1])*np.diff(r))])
        return np.interp(rp, r, cum[-1] - cum)*thermal

    src = SF.SphericalSource(r, M_b)
    ops = np.zeros((len(rp), len(widths)))
    for j, w in enumerate(widths):
        if w >= CR.CLUSTER_MIN_WIDTH:
            ops[:, j] = pressure_of(src.g_mem(r, w))
    gN = G*M_b/r**2
    h = c['hydro_mass']
    M_nfw = _loglog_interp(r, h['RADIUS'], h['M_NFW'], left_power=2.)
    return dict(name=name, r=r, n_e=n_e, M_gas=M_gas, M_star=M_star, M_b=M_b, star_source=star_source, R500=R500, r_out=float(r_out),
                rp=rp, Pobs=Pobs, eP=eP, kinds=kinds, thermal=thermal, ops=ops, source=src,
                P_N=pressure_of(gN), P_nfw=pressure_of(G*M_nfw/r**2),
                P_pm1_repo=pressure_of(gN + np.sqrt(A_STAR_REPO*gN)), P_pm1_fit=pressure_of(gN + np.sqrt(A_STAR_PM1*gN)),
                gas_mass_check=float(np.interp(R500, r, M_int)/np.interp(R500, r, M_gas)), pressure_of=pressure_of)


def cluster_block(extract, widths, names=None, **kw):
    frac = stellar_fraction_profile(extract)
    cls = [build_cluster(n, c, widths, frac_profile=frac, **kw) for n, c in extract['clusters'].items() if names is None or n in names]
    A = np.vstack([c['ops'] for c in cls])
    y = np.concatenate([c['Pobs'] - c['P_N'] for c in cls])
    s = np.concatenate([c['eP'] for c in cls])
    nu = np.zeros((len(y), len(cls)))
    off = 0
    for k, c in enumerate(cls):
        nu[off:off + len(c['rp']), k] = c['thermal']
        off += len(c['rp'])
    return CR.Block('clusters', A, y, s, nuisance=nu, meta=dict(clusters=cls, frac_profile=frac))


def fit_boundary(c, P_model):
    """chi2 of a fixed model with the boundary pressure fitted (nonnegative weighted mean)."""
    w = 1/c['eP']**2
    pout = max(float(np.sum(w*(c['Pobs'] - P_model*c['thermal']))/np.sum(w*c['thermal']**2)), 0.)
    return float(np.sum(((P_model*c['thermal'] + pout*c['thermal'] - c['Pobs'])/c['eP'])**2)), pout


def cluster_references(block):
    out = {}
    for label, key in (('release_nfw_floor', 'P_nfw'), ('newtonian_baryons', 'P_N'), ('pm1_law_repo_a_star', 'P_pm1_repo'), ('pm1_law_fitted_a_star', 'P_pm1_fit')):
        per = {c['name']: fit_boundary(c, c[key])[0] for c in block.meta['clusters']}
        tot = float(sum(per.values()))
        out[label] = dict(chi2=tot, chi2_per_point=tot/block.N, per_cluster_chi2_per_point={k: v/len(c['rp']) for (k, v), c in zip(per.items(), block.meta['clusters'])})
    return out


def shell_decomposition(c, widths, amplitudes, edges=(.5, 1., 2.)):
    """Fraction of the model's memory force at each pressure radius coming from source shells at r' < r/2,
    r/2..r, r..2r and > 2r."""
    src = c['source']
    nodes = src.nodes
    force = np.zeros((len(c['rp']), len(nodes)))
    for w, a in zip(widths, amplitudes):
        if a > 0 and w >= CR.CLUSTER_MIN_WIDTH:
            force += -a*SF.shell_kernel_dr(c['rp'], nodes, w)*src.dM[None, :]
    total = force.sum(axis=1)
    out = []
    for i, ri in enumerate(c['rp']):
        bins = [nodes < edges[0]*ri, (nodes >= edges[0]*ri) & (nodes < edges[1]*ri), (nodes >= edges[1]*ri) & (nodes < edges[2]*ri), nodes >= edges[2]*ri]
        out.append([float(force[i, b].sum()/total[i]) if total[i] != 0 else 0. for b in bins])
    return dict(radii_kpc=c['rp'].tolist(), fractions_inner_half_to_twice_beyond=out, total_memory_force=total.tolist())


# ================================================================ galaxies: SPARC and the Milky Way
def _galaxy_arrays(gal, upsilon_disk=I.UPSILON_DISK):
    rm = gal['rotmod']
    keep = rm[:, 0] > 0
    R, vobs, ev = rm[keep, 0], rm[keep, 1], rm[keep, 2]
    vg, vd, vb = rm[keep, 3], rm[keep, 4], rm[keep, 5]
    gN = np.maximum(np.sign(vg)*vg**2 + upsilon_disk*vd**2 + I.UPSILON_BULGE*vb**2, 0.)/R
    return R, vobs, ev, gN


def galaxy_operator(gal, widths, R, n_grid=3000, r_max=300., exclude_narrow=True):
    comp = BAR.sparc_components(gal)
    Rg = np.geomspace(1e-3, r_max, n_grid)
    sig = comp['sigma_star'](Rg) + (comp['sigma_gas'](Rg) if comp['sigma_gas'] is not None else 0.)
    f = 2*np.pi*Rg*sig
    Mp = np.concatenate([[0.], np.cumsum(.5*(f[1:] + f[:-1])*np.diff(Rg))]) + np.pi*Rg[0]**2*sig[0]
    disk = SF.ProjectedSource(Rg, Mp)
    bul = SF.SphericalSource(Rg, comp['m_bulge'](Rg)) if comp['m_bulge'] is not None else None
    out = np.array([-disk.convolved_dR(R, w) + (bul.g_mem(R, w) if bul is not None else 0.) for w in widths]).T
    if exclude_narrow:
        out[:, np.asarray(widths, float) < CR.GALAXY_MIN_WIDTH] = 0.       # amendment 2
    return out


def sparc_block(widths, split='train', n_grid=3000, upsilon_disk=I.UPSILON_DISK):
    gals = [g for g in I.sparc_galaxies() if g['split'] == split]
    rows, ys, ss, per = [], [], [], []
    for gal in gals:
        R, vobs, ev, gN = _galaxy_arrays(gal, upsilon_disk)
        rows.append(galaxy_operator(gal, widths, R, n_grid))
        ys.append(vobs**2/R - gN)
        ss.append(np.maximum(2*vobs*ev/R, 1e-3*vobs**2/R))
        per.append(dict(name=gal['name'], R=R, vobs=vobs, ev=ev, gN=gN, n=len(R)))
    return CR.Block(f'sparc_{split}', np.vstack(rows), np.concatenate(ys), np.concatenate(ss), meta=dict(galaxies=per, split=split))


def galaxy_scores(block, g_extra):
    """Equal-galaxy RMSE (km/s) and the outer mass-speed slope for a model g = gN + g_extra on the block's points."""
    per, off, rmse, M, vo, vm = block.meta['galaxies'], 0, [], [], [], []
    for p in per:
        sl = slice(off, off + p['n'])
        off += p['n']
        v = np.sqrt(np.maximum(p['R']*(p['gN'] + g_extra[sl]), 0.))
        rmse.append(np.sqrt(np.mean((v - p['vobs'])**2)))
        if p['n'] >= 5 and p['gN'][-1] > 0:
            M.append(p['R'][-1]**2*p['gN'][-1]/G)
            vo.append(p['vobs'][-5:].mean())
            vm.append(v[-5:].mean())
    M, vo, vm = (np.array(x) for x in (M, vo, vm))
    return dict(equal_galaxy_rmse_km_s=float(np.sqrt(np.mean(np.square(rmse)))), galaxies=len(per),
                mass_speed_slope_model=float(np.polyfit(np.log10(M), np.log10(vm), 1)[0]),
                mass_speed_slope_observed=float(np.polyfit(np.log10(M), np.log10(vo), 1)[0]),
                chi2=float(np.sum(((g_extra - block.y)/block.s)**2)))


def galaxy_references(block):
    gN = np.concatenate([p['gN'] for p in block.meta['galaxies']])
    laws = dict(baryons=np.zeros_like(gN), pm1_law_repo_a_star=np.sqrt(A_STAR_REPO*gN), pm1_law_fitted_a_star=np.sqrt(A_STAR_PM1*gN),
                simple_mond_repo_a_star=gN*(.5 + np.sqrt(.25 + A_STAR_REPO/np.maximum(gN, 1e-300))) - gN)
    return {k: galaxy_scores(block, v) for k, v in laws.items()}


def milky_way_source(variant):
    """Projected surface density of the archived Milky Way disks (Miyamoto-Nagai stellar disks integrated over z,
    holed exponential gas) and the Plummer bulge's enclosed mass; a razor-thin treatment of thick disks, stated."""
    components, gas, _ = I._milky_way_definitions()
    disks = [(M, a, b) for M, a, b in components[variant] if a > 0]
    bulges = [(M, b) for M, a, b in components[variant] if a == 0]
    z = np.concatenate([[0.], np.geomspace(1e-3, 100., 600)])

    def sigma(R):
        R = np.atleast_1d(np.asarray(R, float))
        out = np.zeros_like(R)
        for M, a, b in disks:
            zb = np.sqrt(z[None, :]**2 + b*b)
            rho = b*b*M/(4*np.pi)*(a*R[:, None]**2 + (a + 3*zb)*(a + zb)**2)/((R[:, None]**2 + (a + zb)**2)**2.5*zb**3)
            out += 2*np.trapezoid(rho, z, axis=1)
        x = np.maximum(R, 1e-9)
        out += sum(S*1e6*np.exp(-hole/x - x/h) for S, h, hole in gas)
        return out

    m_bulge = (lambda r: sum(M*np.asarray(r, float)**3/(np.asarray(r, float)**2 + b*b)**1.5 for M, b in bulges)) if bulges else None
    return sigma, m_bulge


def milky_way_blocks(widths, n_grid=3000):
    """The 38 Eilers bins for baselines I and II; the disks are treated as razor-thin with their projected surface
    density (a stated approximation for widths below the scale heights)."""
    out = []
    for variant in ('I', 'II'):
        run = next(x for x in I.milky_way_runs() if x['baryons'] == variant and abs(x['rd'] - 2.6) < 1e-9 and abs(x['lf'] - 1) < 1e-9)
        R, y, vb = np.asarray(run['R'], float), np.asarray(run['y'], float), np.asarray(run['vb'], float)
        Rg = np.geomspace(1e-3, 300., n_grid)
        sigma, m_bulge = milky_way_source(variant)
        f = 2*np.pi*Rg*sigma(Rg)
        Mp = np.concatenate([[0.], np.cumsum(.5*(f[1:] + f[:-1])*np.diff(Rg))]) + np.pi*Rg[0]**2*sigma(Rg[:1])[0]
        disk = SF.ProjectedSource(Rg, Mp)
        bul = SF.SphericalSource(Rg, m_bulge(Rg)) if m_bulge is not None else None
        A = np.array([-disk.convolved_dR(R, w) + (bul.g_mem(R, w) if bul is not None else 0.) for w in widths]).T
        gN = vb**2/R
        out.append(CR.Block(f'milky_way_{variant}', A, y**2/R - gN, np.maximum(2*y*5./R, 1e-3*y**2/R), meta=dict(R=R, y=y, gN=gN, variant=variant)))
    return out


# ================================================================ lenses: CR-2's measurement interface
def read(folder, file='results.json'):
    return json.loads((RESULTS/folder/file).read_text(encoding='utf-8'))


def chi_flrw(z):
    return C_KMS/70*quad(lambda t: 1/np.sqrt(.3*(1 + t)**3 + .7), 0, z, epsabs=1e-9)[0]


def geometry(zl, zs, scenario):
    """One record per scenario: angular-diameter distance to the lens (kpc), D_ls/D_s, and the luminosity-distance
    factor (D_L/D_L,FLRW)^2 that rescales a published stellar mass."""
    if scenario == 'G1_flat_FLRW':
        DM = chi_flrw
        Dl, ratio, DL = DM(zl)/(1 + zl)*1000, 1 - DM(zl)/DM(zs), DM(zl)*(1 + zl)
    elif scenario == 'G2_coscaling_coasting':
        DM = lambda z: math.log1p(z)/ALPHA0
        Dl, ratio, DL = DM(zl)/(1 + zl)*1000, 1 - DM(zl)/DM(zs), DM(zl)*(1 + zl)
    elif scenario == 'PF1_static_euclidean':
        D = lambda z: math.log1p(z)/ALPHA0
        Dl, ratio, DL = D(zl)*1000, 1 - D(zl)/D(zs), D(zl)*(1 + zl)
    else:
        raise ValueError(scenario)
    return dict(scenario=scenario, Dl_kpc=Dl, Dls_over_Ds=ratio, DL_Mpc=DL, mass_factor=(DL/(chi_flrw(zl)*(1 + zl)))**2)


class LensSystem:
    """A lens through CR-2's data and model construction, with the geometry a parameter and the memory forces
    added as rows per width (per unit Lambda per 1e11 Msun of stars)."""
    DATA = PILOT = PROFILES = OBS = CFG = PUB = COND = None

    @classmethod
    def load(cls):
        if cls.DATA is None:
            cls.DATA = {d['Name']: d for d in read('slacs-resolved-input-audit')['systems']}
            cls.PILOT = {r['Name']: r for r in read('slacs-motion-lensing-pilot')['rows'] if r['model'] == 'empirical_extra' and r['cutoff_in_a'] == 20}
            cls.PROFILES = {r['Name']: r for r in read('slacs-light-profile-audit')['rows']}
            cls.OBS = {r['Name']: r for r in read('lensing-data-readiness', 'lens-observations-and-image-models.json')}
            cls.COND = {r['Name']: r for r in read('lensing-data-readiness', 'conditional-geometry.json')}
            cls.CFG = read('slacs-outer-bin-check', 'protocol.json')
            cls.PUB = {}
            for r in read('lens-photometric-audit', 'normalization-sensitivity.json'):
                cls.PUB.setdefault((r['Name'], r['imf']), r['published_log10_stellar_mass'])

    def __init__(self, name, scenario='G1_flat_FLRW', widths=CR.WIDTHS, psf_scale=1.):
        self.load()
        self.name, self.scenario, self.widths = name, scenario, np.asarray(widths, float)
        item, obs = self.DATA[name], self.OBS[name]
        g = self.geo = geometry(obs['zFG'], obs['zBG'], scenario)
        dl = g['Dl_kpc']
        a = self.PILOT[name]['scale_a_kpc']*dl/(self.COND[name]['conditional_Dl_Mpc']*1000)
        edges = np.r_[item['inner_arcsec'], item['outer_arcsec'][-1]]*dl/ARCSEC
        psf = psf_scale*item['psf_fwhm_arcsec']*dl/ARCSEC/np.sqrt(8*np.log(2))
        comps = [dict(R=q['R_arcsec']*dl/ARCSEC, n=q['n'], amp=q['amp_at_R'], bn=q['bn']) for q in self.PROFILES[name]['components']]
        self.components = self.PROFILES[name]['components']
        self.parts = [q['total_light_fraction'] for q in self.components]
        self.model = model = ComponentModel(a, edges, psf, 0, .5, 1, 20, comps)
        r = model.r
        self.starforce = model.forces[0].copy()                       # G 1e11 frac / r^2
        self.frac = self.starforce*r*r/(G*1e11)
        self.theta = self.PILOT[name]['catalog_SIE_arcsec']/ARCSEC
        self.bE, self.ratio = self.theta*dl, g['Dls_over_Ds']
        self.need = self.theta/g['Dls_over_Ds']
        self.S = quad(lambda t: model.mass_fraction(self.bE/np.cos(t))/(self.bE/np.cos(t)), 0, np.pi/2, epsabs=1e-9, epsrel=1e-8)[0]
        self.starbend = 4*G*1e11/C_KMS**2*self.S                       # radians per 1e11 Msun
        self.y = np.array(item['vrms_kms'])
        self.cov = np.array(item['covariance_kms_squared'])
        self.fac = cho_factor(self.cov)
        self.pop = {imf: 10**self.PUB[name, imf]*g['mass_factor'] for imf in ('Chabrier', 'Salpeter')}
        self.beta_bounds = tuple(self.CFG['constant_beta_bounds'])
        # the memory rows: the written field of the deprojected stars, per unit Lambda per 1e11 Msun
        src = SF.SphericalSource(r, 1e11*self.frac)
        self.memforce = np.array([src.g_mem(r, w) for w in self.widths])
        Rg = np.geomspace(1e-4, 1e5, 8000)
        proj = SF.ProjectedSource(Rg, 1e11*SF.sersic_projected_fraction(self.components, self.parts, dl, Rg))
        self.membend = np.array([4*G*proj.projected_equivalent_mass([self.bE], w)[0]/(C_KMS**2*self.bE) for w in self.widths])
        # linearised covariance of V^2: D C D with D = diag(2 V)
        D = np.diag(2*self.y)
        self.fac2 = cho_factor(D@self.cov@D)
        self.rows_cache = {}

    # -- the measured second moments and the Einstein bend, both linear in the amplitudes at fixed beta and mass
    def coefficients(self, beta):
        if beta not in self.rows_cache:
            self.model.forces = np.vstack([self.starforce, self.memforce])
            self.rows_cache[beta] = self.model.coefficients(beta)      # rows: stars, then one per width; per bin
        return self.rows_cache[beta]

    def chi2_v(self, v2_model):
        """chi2 of predicted V_rms against the measured, exactly as CR-2 (in V, with the release covariance)."""
        if np.any(v2_model <= 0):
            return 1e30
        e = np.sqrt(v2_model) - self.y
        return float(e@cho_solve(self.fac, e))

    def stars_only(self):
        """CR-2's benchmark: the stellar mass from the Einstein radius, beta fitted, no extra field."""
        mstar = self.need/self.starbend
        self.model.forces = np.array([self.starforce, np.zeros_like(self.starforce)])
        ev = lambda beta: self.chi2_v(self.model.coefficients(beta)[0]*mstar)
        best = None
        for b0 in self.CFG['starts_beta']:
            from scipy.optimize import minimize
            o = minimize(lambda p: ev(p[0]), [b0], bounds=[self.beta_bounds], method='L-BFGS-B')
            if o.success and (best is None or o.fun < best.fun):
                best = o
        self.rows_cache = {}
        return dict(chi2=float(best.fun), beta=float(best.x[0]), stellar_mass_Msun=float(mstar*1e11))

    def kinematic_rows(self, m11, beta):
        """Whitened linear rows: V^2 = m (c_star + sum_j Lambda_j c_j). Returns (A_rows over widths, target, whitening)."""
        co = self.coefficients(beta)
        A = (m11*co[1:]).T                                             # bins x widths
        y = self.y**2 - m11*co[0]
        return A, y

    def einstein_row(self, m11):
        """need = m starbend + sum_j Lambda_j m membend_j: one linear row with a 3% error on need."""
        return m11*self.membend, self.need - m11*self.starbend, .03*self.need

    def fit_beta(self, m11, amplitudes):
        ev = lambda beta: self.chi2_v(m11*(self.coefficients(beta)[0] + amplitudes@self.coefficients(beta)[1:]))
        lo, hi = self.beta_bounds
        grid = np.linspace(lo, hi, 11)
        vals = [ev(b) for b in grid]
        i = int(np.argmin(vals))
        res = minimize_scalar(ev, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 10)]), method='bounded', options={'xatol': 1e-6})
        return (float(res.fun), float(res.x)) if res.fun < vals[i] else (float(vals[i]), float(grid[i]))

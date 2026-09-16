"""PM-2A: which shared field reproduces PM-1's acceleration law, from an independently specified source.

    python pm2a.py [--output-dir DIR] [--canonical]

Stage A separates the four masses and gates the source construction before any fit:
    M_force(R) = R^2 g_N/G   reconstructed from a force (what PM-1 used)
    M_cyl(<R)                the source density integrated over a cylinder
    M_sph(<r)                the same density integrated over a sphere
    M_total                  including the source model's stated outer continuation

The source is RPG-1's archived construction, kept intact: ln Sigma interpolated linearly in radius over
the retained SBdisk samples, continued beyond the final sample as ln Sigma(R_end) - (x - R_end)/rd, a
separately normalized gas disk and the existing bulge. Because ln Sigma is piecewise linear in radius the
stellar cylindrical mass is analytic segment by segment, and the continuation closes to
2 pi Sigma(R_end) rd (R_end + rd), so Stage A uses exact integrals, checks them against an independently
converged Gauss-Legendre integration, and measures the archived quadrature against both rather than
assuming it.

The primary source keeps the continuation; the declared sensitivity removes only the extrapolated stellar
tail, renormalizing nothing. Neither is chosen by which fits better.

Regenerates pm2a-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
sys.path.insert(0, str(HERE.parent/'capture-to-orbit'))
import inputs as I  # noqa: E402
sys.path.insert(0, str(HERE.parent/'radiation-polarized-gravity'))
import baryons as BAR  # noqa: E402
sys.path.insert(0, str(HERE))
import pm1 as P  # noqa: E402
import fields as F  # noqa: E402

G = P.G
TOL_MASS = 1e-10          # annular refinement conserves the source mass
TOL_ROWS = 1e-10          # end-to-end invariance of the prediction under velocity-row removal
TOL_INTEGRAL = 1e-8       # closed form against an independently converged quadrature
EXPECTED_SOURCES = 149    # a run that built fewer has not passed, whatever its gates report


# ---------------------------------------------------------------- analytic sources with known answers
def exponential_disk(sigma0=1e8, h=3.):
    """Sigma = sigma0 exp(-R/h); M_cyl(<R) = 2 pi sigma0 h^2 [1 - (1 + R/h) exp(-R/h)]."""
    return dict(name='analytic exponential disk', h=h, sigma0=sigma0,
                sigma=lambda R: sigma0*np.exp(-np.asarray(R, float)/h),
                m_cyl=lambda R: 2*np.pi*sigma0*h**2*(1 - (1 + np.asarray(R, float)/h)
                                                     * np.exp(-np.asarray(R, float)/h)),
                m_total=2*np.pi*sigma0*h**2)


# ---------------------------------------------------------------- exact integrals of the archived profile
def _segment_mass(a, b, s_a, k):
    """int_a^b 2 pi x Sigma dx for Sigma(x) = s_a exp(k (x - a)), in closed form."""
    a, b, s_a, k = (np.asarray(v, float) for v in (a, b, s_a, k))
    flat = np.abs(k) < 1e-13
    kk = np.where(flat, 1., k)
    e = np.exp(np.clip(kk*(b - a), -700, 700))
    exact = 2*np.pi*s_a*(e*(b/kk - 1/kk**2) - (a/kk - 1/kk**2))
    return np.where(flat, np.pi*s_a*(b**2 - a**2), exact)


def stellar_profile(d):
    """The archived construction's own samples: radius, ln Sigma at the retained SBdisk rows, and rd."""
    rows, cat = d['rotmod'], d['catalog']
    R = rows[:, 0]
    use = (R > 0) & np.isfinite(rows[:, 6]) & (rows[:, 6] > 0)
    Rs, lnS = R[use], np.log(I.UPSILON_DISK*rows[use, 6]*1e6)
    assert np.all(np.diff(Rs) > 0), 'the SBdisk samples must be strictly increasing'
    return Rs, lnS, float(cat['rd'])


def _cumulative(Rs, lnS):
    """Mass at each sample radius, accumulated exactly segment by segment."""
    S = np.exp(lnS)
    k = np.diff(lnS)/np.diff(Rs)
    inner = np.pi*S[0]*Rs[0]**2
    return np.concatenate(([inner], inner + np.cumsum(_segment_mass(Rs[:-1], Rs[1:], S[:-1], k)))), S, k


def stellar_cyl_mass(Rs, lnS, rd, R):
    """M_cyl(<R) of the archived stellar profile: flat inside Rs[0] because np.interp clamps there,
    piecewise exponential between the samples, endpoint-anchored continuation beyond Rs[-1]."""
    R = np.atleast_1d(np.asarray(R, float))
    C, S, k = _cumulative(Rs, lnS)
    out = np.where(R <= Rs[0], np.pi*S[0]*np.maximum(R, 0.)**2, 0.)
    mid = (R > Rs[0]) & (R <= Rs[-1])
    if mid.any():
        j = np.clip(np.searchsorted(Rs, R[mid], side='right') - 1, 0, len(Rs) - 2)
        out[mid] = C[j] + _segment_mass(Rs[j], R[mid], S[j], k[j])
    far = R > Rs[-1]
    if far.any():
        out[far] = C[-1] + _segment_mass(Rs[-1], R[far], S[-1], -1./rd)
    return out


def stellar_tail_mass(Rs, lnS, rd):
    """The whole extrapolated tail beyond R_star_source_end: 2 pi Sigma(R_end) rd (R_end + rd)."""
    return float(2*np.pi*np.exp(lnS[-1])*rd*(Rs[-1] + rd))


def converged_interior(Rs, lnS, nodes=(8, 16)):
    """An independent Gauss-Legendre integration of the interpolated stellar interior, at two node counts
    so its own convergence is shown rather than asserted."""
    S = np.exp(lnS)
    k = np.diff(lnS)/np.diff(Rs)
    totals = []
    for n in nodes:
        x, w = np.polynomial.legendre.leggauss(n)
        a, b = Rs[:-1, None], Rs[1:, None]
        xx = .5*(b - a)*(x[None, :] + 1) + a
        f = 2*np.pi*xx*S[:-1, None]*np.exp(k[:, None]*(xx - a))
        totals.append(float(np.pi*S[0]*Rs[0]**2 + np.sum(.5*(b - a)*w[None, :]*f)))
    return totals


# ---------------------------------------------------------------- spherical mass inside a cylinder
def cylinder_from_sphere(m_sph, dm_dr, knots, R, nodes=32):
    """M_cyl(R) = M_sph(R) + int_R^r_max M'(r)[1 - sqrt(1 - R^2/r^2)] dr for a spherical component.

    A shell at r > R puts the fraction 1 - sqrt(1 - R^2/r^2) of itself inside the cylinder. The
    substitution u = sqrt(r^2 - R^2) turns the integrand into u(r - u)/r^2 times M'(r), removing the
    square-root singularity of its derivative at r = R (correction 8).
    """
    R = np.atleast_1d(np.asarray(R, float))
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    a, b = np.asarray(knots[:-1], float)[None, :], np.asarray(knots[1:], float)[None, :]
    lo = np.maximum(a, R[:, None])
    live = b > lo
    ua = np.sqrt(np.maximum(lo*lo - R[:, None]**2, 0.))
    ub = np.sqrt(np.maximum(b*b - R[:, None]**2, 0.))
    half = np.where(live, .5*(ub - ua), 0.)
    u = half[..., None]*(xg + 1) + ua[..., None]
    r = np.sqrt(R[:, None, None]**2 + u*u)
    w = dm_dr(r)*u*(r - u)/(r*r)
    return np.asarray(m_sph(R), float) + np.sum(half*np.sum(wg*w, axis=-1), axis=1)


def plummer(M=1., b=1.):
    """M_sph = M r^3/(r^2+b^2)^{3/2}; its cylindrical projection is M R^2/(R^2+b^2) in closed form."""
    return dict(name='Plummer sphere', M=M, b=b,
                m_sph=lambda r: M*np.asarray(r, float)**3/(np.asarray(r, float)**2 + b*b)**1.5,
                dm_dr=lambda r: 3*M*b*b*np.asarray(r, float)**2/(np.asarray(r, float)**2 + b*b)**2.5,
                m_cyl=lambda R: M*np.asarray(R, float)**2/(np.asarray(R, float)**2 + b*b),
                g_N=lambda r: P.G*M*np.asarray(r, float)/(np.asarray(r, float)**2 + b*b)**1.5)


def bulge_profile(d):
    """The archived bulge's own samples, exactly as baryons.py builds them."""
    rows = d['rotmod']
    R = rows[:, 0]
    bul = (R > 0) & (rows[:, 5] > 0)
    if not bul.any():
        return None, None
    return R[bul], np.maximum.accumulate(I.UPSILON_BULGE*R[bul]*rows[bul, 5]**2/I.G)


def bulge_cyl_mass(src, R, nodes=32):
    """The bulge's mass inside a cylinder of radius R. Its enclosed mass is interpolated linearly in
    ln r between the samples, as r^3 inside the first and constant beyond the last, so M'(r) is c/r on
    each segment, 3 m0 r^2/Rb0^3 in the core, and zero beyond the support."""
    Rb, mb = src['Rb'], src['mb']
    if Rb is None:
        return np.zeros_like(np.atleast_1d(np.asarray(R, float)))
    c = np.diff(mb)/np.diff(np.log(Rb))

    def dm_dr(r):
        r = np.asarray(r, float)
        out = np.where(r < Rb[0], 3*mb[0]*r**2/Rb[0]**3, 0.)
        mid = (r >= Rb[0]) & (r <= Rb[-1])
        if mid.any():
            j = np.clip(np.searchsorted(Rb, r[mid], side='right') - 1, 0, len(Rb) - 2)
            out[mid] = c[j]/r[mid]
        return out

    return cylinder_from_sphere(lambda r: src['comp']['m_bulge'](np.maximum(r, 1e-12)), dm_dr,
                                np.concatenate(([0.], Rb)), R, nodes)


def thick_disk_inner_mass(sigma, r_in, h, nodes=64):
    """Mass inside the sphere of radius r_in for rho = Sigma(R) exp(-|z|/h)/(2h), exact for a varying
    Sigma: 2 pi int_0^{r_in} R Sigma(R) [1 - exp(-sqrt(r_in^2 - R^2)/h)] dR (correction 9). The inherited
    builder returns Sigma(0) pi r_in^2 instead, which is a cylinder's mass at constant Sigma."""
    xg, wg = np.polynomial.legendre.leggauss(nodes)
    Rq = .5*r_in*(xg + 1)
    z = np.sqrt(np.maximum(r_in*r_in - Rq*Rq, 0.))
    return float(2*np.pi*np.sum(.5*r_in*wg*Rq*sigma(Rq)*(1 - np.exp(-z/h))))


# ---------------------------------------------------------------- one frozen source per galaxy
def source_fingerprint(src):
    """A hash of the live source specification, recomputed rather than read back from a cached field, and
    covering the bulge's whole radial profile rather than only its total (plus the vertical prescription
    once one exists). Comparing a stored string with itself could not detect a later mutation."""
    spec = dict(name=src['name'], Rs=np.asarray(src['Rs']).tolist(), lnS=np.asarray(src['lnS']).tolist(),
                rd=src['rd'], m_gas=src['m_gas'], h_gas=src['h_gas'],
                Rb=None if src['Rb'] is None else np.asarray(src['Rb']).tolist(),
                mb=None if src['mb'] is None else np.asarray(src['mb']).tolist())
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


def build_source(d):
    """Build the source once, from the full source inputs, and hash the numbers that define it."""
    comp = BAR.sparc_components(d)
    Rs, lnS, rd = stellar_profile(d)
    r_end = float(Rs[-1])
    m_star_interior = float(stellar_cyl_mass(Rs, lnS, rd, np.array([r_end]))[0])
    m_tail = stellar_tail_mass(Rs, lnS, rd)
    m_gas = float(I.HELIUM*d['catalog']['MHI9']*1e9) if comp['sigma_gas'] is not None else 0.
    m_bulge = float(comp['m_bulge'](np.array([1e6]))[0]) if comp['m_bulge'] is not None else 0.
    Rb, mb = bulge_profile(d)
    src = dict(name=d['name'], comp=comp, Rs=Rs, lnS=lnS, rd=rd, R_star_source_end=r_end,
               m_star_interior=m_star_interior, m_star_tail=m_tail, m_gas=m_gas,
               h_gas=float(comp['h_gas']), m_bulge=m_bulge, Rb=Rb, mb=mb)
    src['hash'] = source_fingerprint(src)
    src['m_total_primary'] = m_star_interior + m_tail + m_gas + m_bulge
    src['m_total_tail_removed'] = m_star_interior + m_gas + m_bulge
    src['extrapolated_stellar_fraction'] = m_tail/(m_star_interior + m_tail)
    src['extrapolated_total_fraction'] = m_tail/src['m_total_primary']
    return src


def m_cyl(src, R, model='primary'):
    """The source's cylindrical mass inside R: stars (continued or truncated), gas and bulge."""
    R = np.atleast_1d(np.asarray(R, float))
    x = np.minimum(R, src['R_star_source_end']) if model == 'tail_removed' else R
    out = stellar_cyl_mass(src['Rs'], src['lnS'], src['rd'], x)
    if src['h_gas'] > 0:
        out = out + src['m_gas']*(1 - (1 + R/src['h_gas'])*np.exp(-R/src['h_gas']))
    if src['comp']['m_bulge'] is not None:
        out = out + bulge_cyl_mass(src, R)      # projected onto the cylinder, not M_sph (correction 8)
    return out


def sigma_of(src, model='primary'):
    """The source's surface density, for quadrature diagnostics; the bulge is not a surface density and is
    added as its own enclosed mass wherever a cumulative mass is wanted."""
    def sigma(x):
        x = np.atleast_1d(np.asarray(x, float))
        out = np.array(src['comp']['sigma_star'](x), float)
        if model == 'tail_removed':
            out[x > src['R_star_source_end']] = 0.
        if src['h_gas'] > 0:
            out = out + src['comp']['sigma_gas'](x)
        return out
    return sigma


def m_total(src, model='primary'):
    return src['m_total_primary'] if model == 'primary' else src['m_total_tail_removed']


def predict_collective(src, v_bar, radii, beta, model='primary'):
    """PM-1's collective rule, now from the source's own mass: v^2 = v_bar^2 + beta M_cyl/sqrt(M_total)."""
    return np.sqrt(np.maximum(np.asarray(v_bar, float)**2
                              + beta*m_cyl(src, radii, model)/np.sqrt(m_total(src, model)), 0.))


# ---------------------------------------------------------------- gate 1: mass-preserving refinement
def _refinement_rows(name, cum, sigma, r_out, splits):
    exact = float(cum(np.array([r_out]))[0])
    rows, quad = [], []
    for n in splits:
        edges = np.linspace(0, r_out, 40*n + 1)
        m = np.diff(cum(edges))
        rows.append(dict(sub_rings=n, annuli=int(len(m)), relative_error=float(abs(m.sum()/exact - 1)),
                         any_negative=bool(np.any(m < -1e-9*exact))))
        mid = .5*(edges[1:] + edges[:-1])
        quad.append(dict(sub_rings=n, midpoint_relative_error=float(
            abs(np.sum(2*np.pi*mid*sigma(mid)*np.diff(edges))/exact - 1))))
    ratios = [quad[i]['midpoint_relative_error']/max(quad[i + 1]['midpoint_relative_error'], 1e-300)
              for i in range(len(quad) - 1)]
    return dict(source=name, r_out=float(r_out), rows=rows,
                worst_relative_error=max(r['relative_error'] for r in rows),
                clipped=any(r['any_negative'] for r in rows),
                midpoint_quadrature=dict(rows=quad, refinement_ratios=ratios))


def gate_mass_preserving(sources, splits=(1, 2, 4, 8, 16)):
    """Refining the annuli must conserve the source mass and must not clip it. The midpoint rule's own
    convergence is reported alongside as a diagnostic: telescoping conserves by construction, which is not
    evidence that the underlying cumulative mass is accurate."""
    a = exponential_disk()
    analytic = _refinement_rows(a['name'], a['m_cyl'], a['sigma'], 10*a['h'], splits)
    gal = [_refinement_rows(name, lambda R, s=src: m_cyl(s, R), sigma_of(src),
                            2*src['R_star_source_end'], splits)
           for name, src in sources.items()]
    worst = max([analytic['worst_relative_error']] + [g['worst_relative_error'] for g in gal])
    clipped = analytic['clipped'] or any(g['clipped'] for g in gal)
    return dict(analytic=analytic, galaxies=gal[:5], galaxy_count=len(gal), worst_relative_error=worst,
                clipped=clipped, tolerance=TOL_MASS, passed=bool(worst <= TOL_MASS and not clipped),
                note='the galaxy rows carry the bulge too, whose enclosed mass is interpolated, so the gate '
                     'covers the whole source and not only the analytic disks')


# ---------------------------------------------------------------- gate 2: end-to-end row removal
def gate_row_removal_end_to_end(sources, gals, data, keep_fraction=.8, control_galaxies=8):
    """The predicted force, not a stored mass function, must be invariant when velocity rows are masked.

    The source is built once from the full frozen source inputs and passed separately from the mask; only
    the evaluation radii change, and the source's hash is checked across the call. A positive control
    rebuilds the source from the shortened rotmod array -- the mistake this gate exists to catch -- and
    must fail the same comparison, so that passing is evidence about the pipeline and not about the test.
    """
    rows, control = [], []
    for g in gals:
        src = sources.get(g['name'])
        if src is None or len(g['R']) < 8:
            continue
        before = source_fingerprint(src)            # recomputed from the live spec, not read back
        keep = g['R'] <= keep_fraction*g['R'][-1]
        if keep.sum() < 3:
            continue
        full = predict_collective(src, g['v_bar'], g['R'], .1)
        masked = predict_collective(src, g['v_bar'][keep], g['R'][keep], .1)
        rows.append(dict(galaxy=g['name'], common_radii=int(keep.sum()),
                         max_relative_change=float(np.max(np.abs(masked/full[keep] - 1))),
                         hash_unchanged=bool(source_fingerprint(src) == before == src['hash'])))
        if len(control) < control_galaxies:
            d = data[g['name']]
            cut = keep_fraction*d['rotmod'][:, 0].max()
            short = dict(d, rotmod=d['rotmod'][d['rotmod'][:, 0] <= cut])
            try:
                rebuilt = build_source(short)
                bad = predict_collective(rebuilt, g['v_bar'][keep], g['R'][keep], .1)
                control.append(dict(galaxy=g['name'], hash_changed=bool(rebuilt['hash'] != src['hash']),
                                    max_relative_change=float(np.max(np.abs(bad/full[keep] - 1))),
                                    total_mass_ratio=float(rebuilt['m_total_primary']
                                                           / src['m_total_primary'])))
            except Exception as exc:
                control.append(dict(galaxy=g['name'], rebuild_failed=f'{type(exc).__name__}: {exc}',
                                    hash_changed=True, max_relative_change=float('inf')))
    worst = max(r['max_relative_change'] for r in rows)
    ctrl_min = min(c['max_relative_change'] for c in control)
    detected = bool(ctrl_min > TOL_ROWS and all(c['hash_changed'] for c in control))
    return dict(galaxies=len(rows), keep_fraction=keep_fraction, worst_relative_change=worst,
                all_hashes_unchanged=bool(all(r['hash_unchanged'] for r in rows)),
                positive_control=dict(
                    galaxies=len(control), rows=control, smallest_relative_change=ctrl_min,
                    all_hashes_changed=bool(all(c['hash_changed'] for c in control)), detected=detected,
                    note='rebuilding sparc_components from a shortened rotmod array moves the stellar '
                         'profile, the gas scale and the bulge at once; that is what must not happen, and '
                         'the control shows this comparison would see it'),
                tolerance=TOL_ROWS,
                passed=bool(worst <= TOL_ROWS and all(r['hash_unchanged'] for r in rows) and detected),
                note="the force-proxy contrast reported in PM-1's terms stays a separate diagnostic: a "
                     'model can show that contrast while its prediction path still depends on where the '
                     'velocity curve stops')


# ---------------------------------------------------------------- gate 3: integral agreement
def gate_integral_agreement(sources, cell_galaxies=4):
    """The closed-form masses against an independently converged Gauss-Legendre integration, with the
    grid's own cell-mass construction measured against the same totals and reported, not gated."""
    a = exponential_disk()
    x, w = np.polynomial.legendre.leggauss(400)
    r_out = 60*a['h']
    xx = .5*r_out*(x + 1)
    analytic = float(abs(np.sum(.5*r_out*w*2*np.pi*xx*a['sigma'](xx))/a['m_total'] - 1))
    rows = []
    for name, src in sources.items():
        lo, hi = converged_interior(src['Rs'], src['lnS'])
        rows.append(dict(galaxy=name, closed_form=src['m_star_interior'], quadrature=hi,
                         relative_difference=float(abs(hi/src['m_star_interior'] - 1)),
                         quadrature_self_convergence=float(abs(hi/lo - 1))))
    worst = max(r['relative_difference'] for r in rows)
    worst_self = max(r['quadrature_self_convergence'] for r in rows)
    cells = []
    for name, src in list(sources.items())[:cell_galaxies]:
        h = .1*src['comp']['rd']
        grid, _ = BAR.sparc_grid(dict(r=np.array([src['R_star_source_end']])), src['comp'],
                                 src['m_total_primary'])
        m, m_in = BAR.sparc_masses(grid, src['comp'], h)
        total = 2*float(m.sum()) + m_in                       # the grid is a hemisphere; m_in is the full core
        r_in = float(grid.rf[0])
        corrected = thick_disk_inner_mass(sigma_of(src), r_in, h)
        fixed = total - m_in + corrected
        cells.append(dict(galaxy=name, cell_construction=total, exact=src['m_total_primary'],
                          relative_difference=float(abs(total/src['m_total_primary'] - 1)),
                          cells=int(m.size), r_in=r_in, h=h,
                          inner_mass_inherited=m_in, inner_mass_spherical=corrected,
                          inner_ratio=float(m_in/corrected), leading_ratio=float(3*h/(2*r_in)),
                          inner_share_of_total=float(m_in/total),
                          cell_total_with_spherical_inner=fixed,
                          relative_difference_after_correction_9=float(abs(fixed/src['m_total_primary'] - 1))))
    pl = plummer()
    Rt = np.array([.5, 1., 2., 5.])*pl['b']
    proj = cylinder_from_sphere(pl['m_sph'], pl['dm_dr'], np.concatenate(([0.], np.geomspace(1e-3, 1e4, 60))), Rt)
    coarse = cylinder_from_sphere(pl['m_sph'], pl['dm_dr'],
                                  np.concatenate(([0.], np.geomspace(1e-3, 1e4, 60))), Rt, nodes=16)
    bulge = dict(source=pl['name'], rows=[
        dict(R_over_b=float(r/pl['b']), m_sph=float(pl['m_sph'](np.array([r]))[0]),
             m_cyl_closed_form=float(pl['m_cyl'](np.array([r]))[0]), m_cyl_projected=float(p),
             relative_difference=float(abs(p/pl['m_cyl'](np.array([r]))[0] - 1)),
             spherical_underestimates_by=float(1 - pl['m_sph'](np.array([r]))[0]/pl['m_cyl'](np.array([r]))[0]))
        for r, p in zip(Rt, proj)],
        node_refinement=float(np.max(np.abs(proj/coarse - 1))),
        note='correction 8: substituting M_sph for M_cyl understates the cylinder by 29.2893% at R = b; '
             'in spherical symmetry M_force equals M_sph, and M_cyl >= M_sph whenever the source extends '
             'past the evaluation radius')
    bulge['worst_relative_difference'] = max(r['relative_difference'] for r in bulge['rows'])
    return dict(analytic_disk_relative_error=analytic, bulge_projection=bulge,
                galaxies=len(rows), rows=rows[:8],
                worst_relative_difference=worst, worst_quadrature_self_convergence=worst_self,
                cell_construction=dict(
                    rows=cells, worst=float(max(c['relative_difference'] for c in cells)),
                    inherited_inner_over_spherical=float(np.median([c['inner_ratio'] for c in cells])),
                    worst_after_correction_9=float(max(c['relative_difference_after_correction_9']
                                                       for c in cells)),
                    note='the AQUAL cells spread a thickened disk over spherical shells, so only the total '
                         "is comparable with a cylindrical integral; the difference measures that "
                         "construction's quadrature and its outer truncation and is reported, not gated",
                    correction_9='the inherited builder hands the solver Sigma(0) pi r_in^2 for its inner '
                                 'SPHERE, which is a cylinder mass at constant Sigma; the exact spherical '
                                 'integral is smaller by about 3h/(2 r_in). The two are reported side by '
                                 'side here so this source-volume mismatch is not folded into the '
                                 'quadrature difference above, and PM-2A substitutes the spherical value '
                                 'before any disk comparison'),
                tolerance=TOL_INTEGRAL,
                passed=bool(worst <= TOL_INTEGRAL and analytic <= TOL_INTEGRAL
                            and bulge['worst_relative_difference'] <= TOL_INTEGRAL))


# ---------------------------------------------------------------- diagnostics
def quadrature_diagnostic(sources):
    """The archived _cylinder_mass measured against the exact integrals it approximates."""
    rows = []
    for name, src in sources.items():
        archived = float(BAR.sparc_total_mass(src['comp']))
        rows.append(dict(galaxy=name, archived=archived, exact=src['m_total_primary'],
                         relative_difference=float(abs(archived/src['m_total_primary'] - 1))))
    d = np.array([r['relative_difference'] for r in rows])
    h, unit = 3., []
    s0 = 1/(2*np.pi*h**2)
    for r_out in (h, 3*h, 10*h, 30*h):
        got = BAR._cylinder_mass(lambda z: s0*np.exp(-np.asarray(z, float)/h), r_out)
        ex = 1 - (1 + r_out/h)*np.exp(-r_out/h)
        unit.append(dict(r_out=r_out, relative_error=float(abs(got/ex - 1))))
    return dict(galaxies=len(rows), median_relative_difference=float(np.median(d)),
                max_relative_difference=float(d.max()),
                unit_exponential=dict(rows=unit, worst=float(max(u['relative_error'] for u in unit))),
                tolerance_it_cannot_meet=TOL_INTEGRAL,
                note='a regression diagnostic, not an exact reference: 20,001 geometrically spaced '
                     'trapezoid points carry 1e-8 to 1e-7 relative error on a unit exponential, so two '
                     'quantities produced by that route agreeing cannot establish 1e-8 accuracy')


def collective_tail_algebra(sources):
    """Inside the endpoint, removing the tail leaves M_cyl(<r) untouched and raises the collective term by
    1/sqrt(1 - f), with f the extrapolated fraction of the total: algebra of the prescription."""
    rows = []
    for name, src in sources.items():
        r = np.array([.5*src['R_star_source_end']])
        same = float(abs(m_cyl(src, r, 'tail_removed')[0]/m_cyl(src, r, 'primary')[0] - 1))
        f = src['extrapolated_total_fraction']
        predicted, actual = 1/np.sqrt(1 - f), np.sqrt(m_total(src, 'primary')/m_total(src, 'tail_removed'))
        rows.append(dict(galaxy=name, inner_mass_change=same, extrapolated_total_fraction=f,
                         predicted_force_ratio=float(predicted), actual_force_ratio=float(actual),
                         difference=float(abs(predicted/actual - 1))))
    worst_same = max(r['inner_mass_change'] for r in rows)
    worst_diff = max(r['difference'] for r in rows)
    return dict(galaxies=len(rows), worst_inner_mass_change=worst_same, worst_ratio_difference=worst_diff,
                ten_percent_tail_raises_the_term_by=float(1/np.sqrt(.9) - 1),
                passed=bool(worst_same < 1e-12 and worst_diff < 1e-12),
                note='algebra of the collective prescription, not a substitute for the field calculation')


def four_masses(sources, gals):
    """M_force at the last sampled radius against the source's own cylindrical and total masses."""
    rows = []
    for g in gals:
        src = sources.get(g['name'])
        if src is None or g['M_force'][-1] <= 0:
            continue
        R = float(g['R'][-1])
        rows.append(dict(galaxy=g['name'], R_last=R, M_force=float(g['M_force'][-1]),
                         M_cyl=float(m_cyl(src, np.array([R]))[0]), M_total=src['m_total_primary']))
    cyl = np.array([r['M_force']/r['M_cyl'] for r in rows])
    tot = np.array([r['M_force']/r['M_total'] for r in rows])

    def spread(v):
        return dict(median=float(np.median(v)), tenth=float(np.percentile(v, 10)),
                    ninetieth=float(np.percentile(v, 90)), minimum=float(v.min()), maximum=float(v.max()))
    bulged = [s for s in sources.values() if s['m_bulge'] > 0]
    return dict(galaxies=len(rows), force_over_cylindrical=spread(cyl), force_over_total=spread(tot),
                spherical_mass=dict(
                    exactly_spherical_component='the bulge, whose enclosed mass is spherical by '
                                                'construction', galaxies_with_a_bulge=len(bulged),
                    note='M_sph(<r) of the disks is not defined by the surface density alone: it depends on '
                         'the declared vertical profile, which the field stages fix. It is therefore left '
                         'to stage B, where the spherical law is verified on a genuinely spherical source, '
                         'rather than quoted here from a cylindrical integral'),
                note='M_force is a force-equivalent mass and M_cyl a source integral. In spherical symmetry '
                     'M_force equals the SPHERICAL enclosed mass, by Newton\'s theorems, not the '
                     'cylindrical one: a cylinder also contains material outside the sphere, so M_cyl >= '
                     'M_sph whenever the source extends past the evaluation radius. A ratio away from one '
                     'is therefore a statement about geometry and about whether the force law holds, not a '
                     'defect of the source model (correction 8 replaced an earlier note that said the two '
                     'coincide in spherical symmetry)',
                caveat='v_bar is the archived tabulated ordinary-matter contribution while M_cyl is the '
                       'mass of the source reconstructed here, and RPG-1 documents that the two differ '
                       'materially at some radii, so this ratio mixes disk flattening with that '
                       'reconstruction difference in unknown proportions; separating them needs the '
                       "source's own Newtonian field on the same grid, which is stage C's recomputed "
                       'baseline')


# ---------------------------------------------------------------- stage B: the field equations
def scalar_identities():
    """Each equation's mu, spherical inverse and functional, checked as algebra before any solve."""
    y = np.geomspace(1e-10, 1e10, 20001)
    rows = []
    for eq in (F.NEWTONIAN, F.COMPLETION_I, F.COMPLETION_II_AUX, F.SIMPLE):
        x = eq.x_of_y(y)
        rows.append(dict(equation=eq.name,
                         mu_x_equals_y=float(np.max(np.abs(eq.mu(x)*x/y - 1))),
                         nu_consistent=float(np.max(np.abs(eq.nu(y)*y/x - 1))),
                         weight_matches_numeric=_weight_check(eq)))
    xs = np.geomspace(1e-9, 1e6, 4001)
    h = xs*1e-6
    deriv = [dict(equation=eq.name,
                  max_relative_error=float(np.max(np.abs(((eq.F((xs+h)**2) - eq.F((xs-h)**2))/(2*h))
                                                         / (2*xs*eq.mu(xs)) - 1))))
             for eq in (F.NEWTONIAN, F.COMPLETION_I, F.COMPLETION_II_AUX, F.SIMPLE)]
    yy = np.geomspace(1e-10, 1e10, 20001)
    matched = float(np.max(np.abs(F.COMPLETION_I.F((yy + np.sqrt(yy))**2)
                                  / (F.NEWTONIAN.F(yy**2) + F.COMPLETION_II_AUX.F(yy)) - 1)))
    tiny = float(F.COMPLETION_I.F(1e-24))
    written = 1e-16 - 2e-8 + 2*np.log1p(1e-8)
    return dict(
        constitutive=rows, functional_derivative=deriv,
        matched_spherical_functional=dict(
            max_relative_difference=matched, tolerance=1e-12,
            statement='F_I(y+sqrt y) == y^2 + (2/3) y^(3/2) == F_N + F_aux on the matched solution',
            passed=bool(matched < 1e-12)),
        cancellation=dict(
            completion_I_at_x_1e_12=tiny, expected_two_thirds_x_cubed=2/3*1e-36,
            simple_mu_as_written_at_x_1e_8=float(written),
            simple_mu_series_at_x_1e_8=float(F.SIMPLE.F(1e-16)),
            note='the inherited form returns exactly 0.0 at x = 1e-8 against a true 6.667e-25; '
                 "Completion I's t-form has no subtraction of nearly equal numbers at all"),
        passed=bool(all(r['mu_x_equals_y'] < 1e-12 and r['nu_consistent'] < 1e-12
                        and r['weight_matches_numeric'] < 1e-6 for r in rows)
                    and all(d['max_relative_error'] < 1e-6 for d in deriv) and matched < 1e-12))


def _weight_check(eq, xs=(1e-4, 1e-2, 1., 1e2, 1e4), h=1e-6):
    """w = 1/(1 + dln mu/dln x), against a numerical log-derivative of that equation's own mu."""
    x = np.asarray(xs, float)
    d = (np.log(eq.mu(x*(1+h))) - np.log(eq.mu(x*(1-h))))/(np.log(x*(1+h)) - np.log(x*(1-h)))
    return float(np.max(np.abs(eq.weight(x)/(1/(1 + d)) - 1)))


def spherical_solves(compactness=(1e-2, 1., 1e2), b=1., a_star=1e-2, nr=(300, 600)):
    """Each equation solved on a Plummer sphere, then checked against its own flux identity."""
    out = []
    for c in compactness:
        M = c*a_star*b*b/P.G
        pl = plummer(M=M, b=b)
        r_in, r_out = 1e-3*b, 1e3*b
        runs, grids = {}, {}
        for eq in (F.NEWTONIAN, F.COMPLETION_I, F.COMPLETION_II_AUX, F.SIMPLE):
            runs[eq.name] = F.spherical_solve(eq, pl['m_sph'], a_star, r_in, r_out, nr=nr[0])
            grids[eq.name] = runs[eq.name].grid
        row = dict(compactness=c, a_star=a_star, M=M, equations={})
        for name, s in runs.items():
            eq = F.EQUATIONS[name]
            g = s.grid
            gr = s.gradients(s.phi)[0][1:-1]                 # interior faces only; 0 and -1 are imposed
            rf = g.rf[1:-1][:, None]
            lhs = rf**2*eq.mu(np.maximum(gr, 1e-300)/s.a if s.a > 0 else np.ones_like(gr))*gr
            rhs = P.G*pl['m_sph'](rf)
            flux = float(np.max(np.abs(lhs/rhs - 1)))
            ang = float(np.max(np.ptp(gr, axis=1)/np.maximum(np.abs(np.mean(gr, axis=1)), 1e-300)))
            row['equations'][name] = dict(
                interior_flux_identity=flux, angular_spread=ang, residual_off_pin=s.residual_off_pin,
                residual_including_pin=s.residual, iteration_change=s.iteration_change,
                iterations=len(s.history), converged=bool(s.converged))
        # the two completions must AGREE on the spherical benchmark; g_II adds vectors, not magnitudes
        sI, sN, sA = runs['completion_I'], runs['newtonian'], runs['completion_II_aux']
        rf = sI.grid.rf[1:-1]
        gI = sI.gradients(sI.phi)[0][1:-1][:, -1]
        gII = (sN.gradients(sN.phi)[0][1:-1][:, -1] + sA.gradients(sA.phi)[0][1:-1][:, -1])
        row['completions_agree'] = float(np.max(np.abs(gII/gI - 1)))
        analytic = pl['g_N'](rf) + np.sqrt(a_star*pl['g_N'](rf))
        row['against_analytic_law'] = dict(
            completion_I=float(np.max(np.abs(gI/analytic - 1))),
            completion_II=float(np.max(np.abs(gII/analytic - 1))))
        # off-grid readout is a different numerical path from the face flux, and must be measured as one
        r_off = np.sqrt(sI.grid.rf[1:-1][:-1]*sI.grid.rf[1:-1][1:])
        law = pl['g_N'](r_off) + np.sqrt(a_star*pl['g_N'](r_off))
        err = lambda p: float(np.median(np.abs(np.asarray(p)/law - 1)))
        row['off_grid_readout'] = dict(
            inherited_gradient_at=err(sI.gradient_at(r_off, np.zeros_like(r_off))[0]),
            inherited_midplane_speed=err(sI.midplane_speed(r_off)**2/r_off),
            face_spline=err(sI.radial_gradient_at(r_off)),
            note='the two inherited paths interpolate cell-centred averages, or face values linearly, and '
                 'converge at O(h^2); the face spline keeps the solution\'s own accuracy. Only the face '
                 'spline is gated, and stage C must not read galaxies through the inherited paths')
        # resolution refinement, reported rather than assumed
        fine = F.spherical_solve(F.COMPLETION_I, pl['m_sph'], a_star, r_in, r_out, nr=nr[1])
        rf_f = fine.grid.rf[1:-1]
        gf = fine.gradients(fine.phi)[0][1:-1][:, -1]
        row['refinement'] = dict(
            nr=list(nr), max_relative_change=float(np.max(np.abs(
                np.interp(np.log(rf), np.log(rf_f), gf)/gI - 1))))
        out.append(row)
    worst_flux = max(e['interior_flux_identity'] for r in out for e in r['equations'].values())
    worst_law = max(max(r['against_analytic_law'].values()) for r in out)
    return dict(
        rows=out, worst_interior_flux=worst_flux, worst_against_analytic_law=worst_law,
        worst_completion_disagreement=max(r['completions_agree'] for r in out),
        off_grid_readout=dict(
            face_spline=max(r['off_grid_readout']['face_spline'] for r in out),
            inherited_gradient_at=max(r['off_grid_readout']['inherited_gradient_at'] for r in out),
            inherited_midplane_speed=max(r['off_grid_readout']['inherited_midplane_speed'] for r in out),
            finding='both inherited readouts are three to four orders worse than the field they read, '
                    'and converge only at second order; this is a bottleneck for stage C, where a '
                    'difference between completions smaller than the readout error is unresolved'),
        all_converged=bool(all(e['converged'] for r in out for e in r['equations'].values())),
        tolerance=1e-6,
        passed=bool(worst_flux < 1e-6 and worst_law < 1e-6
                    and max(r['off_grid_readout']['face_spline'] for r in out) < 1e-6
                    and all(e['converged'] for r in out for e in r['equations'].values())))


def equation_distinction():
    """At g_N = a*, the completions give 2.000000 and the simple-mu law 1.618034. A test that treated the
    three as one equation must fail here; requiring the two completions to differ would be a wrong test."""
    y = np.array([1.])
    got = {e.name: float(e.x_of_y(y)[0]) for e in (F.COMPLETION_I, F.COMPLETION_II_AUX, F.SIMPLE)}
    got['completion_II_total'] = 1. + got['completion_II_aux']
    golden = .5*(1 + np.sqrt(5))
    return dict(values=got, expected=dict(completions=2., simple_mu=float(golden)),
                completions_agree=abs(got['completion_I'] - got['completion_II_total']) < 1e-12,
                simple_is_distinct=bool(abs(got['simple_mu']/golden - 1) < 1e-12
                                        and abs(got['simple_mu'] - 2.) > .3),
                passed=bool(abs(got['completion_I'] - 2.) < 1e-12
                            and abs(got['completion_II_total'] - 2.) < 1e-12
                            and abs(got['simple_mu'] - golden) < 1e-12))


def stage_b():
    ident = scalar_identities()
    solves = spherical_solves()
    dist = equation_distinction()
    return dict(
        scalar_identities=ident, spherical_solves=solves, equation_distinction=dist,
        gauge='the outer equatorial cell is pinned at a finite radius; these isolated models have a '
              'logarithmic far field, so Phi(infinity) = 0 is not a valid normalization and every '
              'functional comparison uses the same convention',
        passed=bool(ident['passed'] and solves['passed'] and dist['passed']),
        what_this_is_not='an equation implementation verified on a spherical source. It is not path '
                         'memory, and the solver\'s relaxation time is a numerical quantity that is never '
                         'a physical memory time')


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    gals = P.sparc()
    data = {d['name']: d for d in I.sparc_galaxies()}
    sources, failures = {}, []
    for g in gals:
        try:
            sources[g['name']] = build_source(data[g['name']])
        except Exception as exc:
            failures.append(dict(galaxy=g['name'], error=f'{type(exc).__name__}: {exc}'))
    frac_star = np.array([s['extrapolated_stellar_fraction'] for s in sources.values()])
    frac_tot = np.array([s['extrapolated_total_fraction'] for s in sources.values()])
    mass_gate = gate_mass_preserving(sources)
    row_gate = gate_row_removal_end_to_end(sources, gals, data)
    int_gate = gate_integral_agreement(sources)
    tail = collective_tail_algebra(sources)
    stage = stage_b()
    result = dict(
        experiment='PM-2A: shared-field completions from an independently specified source',
        protocol='protocol-pm2a.md',
        stage='A, the source audit, and B, the field equations; C and D follow',
        source_model=dict(
            primary='the archived stellar interpolation with its endpoint-anchored exponential '
                    'continuation, a separately normalized gas disk, and the existing bulge',
            sensitivity='the extrapolated stellar tail removed, with nothing renormalized',
            endpoint='R_star_source_end, the final retained positive SBdisk sample in the archived source '
                     'construction; not assumed to be the final original photometric measurement',
            independence='independent of the observed rotation speeds being fitted: the gas scale is fitted '
                         'to the tabulated gas contribution and the bulge is reconstructed from the '
                         'tabulated bulge contribution, so this is not a source built from no force inputs',
            selection='neither prescription is selected by which fits better',
            galaxies=len(sources), failures=failures),
        extrapolated=dict(stellar_fraction_median=float(np.median(frac_star)),
                          stellar_fraction_90th=float(np.percentile(frac_star, 90)),
                          stellar_fraction_max=float(frac_star.max()),
                          total_fraction_median=float(np.median(frac_tot)),
                          total_fraction_90th=float(np.percentile(frac_tot, 90)),
                          total_fraction_max=float(frac_tot.max()),
                          largest=sorted(({'galaxy': s['name'],
                                           'total_fraction': s['extrapolated_total_fraction'],
                                           'stellar_fraction': s['extrapolated_stellar_fraction'],
                                           'R_star_source_end': s['R_star_source_end'], 'rd': s['rd']}
                                          for s in sources.values()),
                                         key=lambda r: -r['total_fraction'])[:5]),
        four_masses=four_masses(sources, gals),
        gates=dict(mass_preserving=mass_gate, row_removal_end_to_end=row_gate,
                   integral_agreement=int_gate, collective_tail_algebra=tail),
        stage_B=stage,
        quadrature_diagnostic=quadrature_diagnostic(sources),
        gates_passed=bool(mass_gate['passed'] and row_gate['passed'] and int_gate['passed']
                          and tail['passed'] and not failures and len(sources) == EXPECTED_SOURCES),
        stage_B_passed=stage['passed'],
        input_sha256={'protocol-pm2a.md': hashlib.sha256((HERE/'protocol-pm2a.md').read_bytes()).hexdigest()},
        checks_short_run=dict(mass_worst=mass_gate['worst_relative_error'],
                              rows_worst=row_gate['worst_relative_change'],
                              rows_control_smallest=row_gate['positive_control']['smallest_relative_change'],
                              integral_worst=int_gate['worst_relative_difference'],
                              tail_worst=tail['worst_ratio_difference'],
                              extrapolated_total_median=float(np.median(frac_tot)),
                              stage_B_flux=stage['spherical_solves']['worst_interior_flux'],
                              stage_B_law=stage['spherical_solves']['worst_against_analytic_law'],
                              stage_B_readout=stage['spherical_solves']['off_grid_readout']['face_spline']),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:2200])
    status = evidence_io.finish(args, 'path-memory-pm2a', text, HERE/'pm2a-results.json',
                                ignore={'/runtime_seconds'})
    return status if (result['gates_passed'] and result['stage_B_passed']) else 1


if __name__ == '__main__':
    raise SystemExit(main())

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

G = P.G
TOL_MASS = 1e-10          # annular refinement conserves the source mass
TOL_ROWS = 1e-10          # end-to-end invariance of the prediction under velocity-row removal
TOL_INTEGRAL = 1e-8       # closed form against an independently converged quadrature


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


# ---------------------------------------------------------------- one frozen source per galaxy
def build_source(d):
    """Build the source once, from the full source inputs, and hash the numbers that define it."""
    comp = BAR.sparc_components(d)
    Rs, lnS, rd = stellar_profile(d)
    r_end = float(Rs[-1])
    m_star_interior = float(stellar_cyl_mass(Rs, lnS, rd, np.array([r_end]))[0])
    m_tail = stellar_tail_mass(Rs, lnS, rd)
    m_gas = float(I.HELIUM*d['catalog']['MHI9']*1e9) if comp['sigma_gas'] is not None else 0.
    m_bulge = float(comp['m_bulge'](np.array([1e6]))[0]) if comp['m_bulge'] is not None else 0.
    payload = json.dumps(dict(name=d['name'], Rs=Rs.tolist(), lnS=lnS.tolist(), rd=rd, m_gas=m_gas,
                              h_gas=float(comp['h_gas']), m_bulge=m_bulge), sort_keys=True).encode()
    src = dict(name=d['name'], comp=comp, Rs=Rs, lnS=lnS, rd=rd, R_star_source_end=r_end,
               m_star_interior=m_star_interior, m_star_tail=m_tail, m_gas=m_gas,
               h_gas=float(comp['h_gas']), m_bulge=m_bulge, hash=hashlib.sha256(payload).hexdigest())
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
        out = out + src['comp']['m_bulge'](np.maximum(R, 1e-12))   # its interpolation is in log r
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
        before = src['hash']
        keep = g['R'] <= keep_fraction*g['R'][-1]
        if keep.sum() < 3:
            continue
        full = predict_collective(src, g['v_bar'], g['R'], .1)
        masked = predict_collective(src, g['v_bar'][keep], g['R'][keep], .1)
        rows.append(dict(galaxy=g['name'], common_radii=int(keep.sum()),
                         max_relative_change=float(np.max(np.abs(masked/full[keep] - 1))),
                         hash_unchanged=bool(src['hash'] == before)))
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
        grid, _ = BAR.sparc_grid(dict(r=np.array([src['R_star_source_end']])), src['comp'],
                                 src['m_total_primary'])
        m, m_in = BAR.sparc_masses(grid, src['comp'], .1*src['comp']['rd'])
        total = 2*float(m.sum()) + m_in                       # the grid is a hemisphere; m_in is the full core
        cells.append(dict(galaxy=name, cell_construction=total, exact=src['m_total_primary'],
                          relative_difference=float(abs(total/src['m_total_primary'] - 1)),
                          cells=int(m.size)))
    return dict(analytic_disk_relative_error=analytic, galaxies=len(rows), rows=rows[:8],
                worst_relative_difference=worst, worst_quadrature_self_convergence=worst_self,
                cell_construction=dict(
                    rows=cells, worst=float(max(c['relative_difference'] for c in cells)),
                    note='the AQUAL cells spread a thickened disk over spherical shells, so only the total '
                         "is comparable with a cylindrical integral; the difference measures that "
                         "construction's quadrature and its outer truncation and is reported, not gated"),
                tolerance=TOL_INTEGRAL, passed=bool(worst <= TOL_INTEGRAL and analytic <= TOL_INTEGRAL))


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
                note='M_force is a force-equivalent mass and M_cyl a source integral; they coincide only in '
                     'spherical symmetry, so a ratio away from one is a statement about disk geometry and '
                     'about whether the force law holds, not a defect of the source model',
                caveat='v_bar is the archived tabulated ordinary-matter contribution while M_cyl is the '
                       'mass of the source reconstructed here, and RPG-1 documents that the two differ '
                       'materially at some radii, so this ratio mixes disk flattening with that '
                       'reconstruction difference in unknown proportions; separating them needs the '
                       "source's own Newtonian field on the same grid, which is stage C's recomputed "
                       'baseline')


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
    result = dict(
        experiment='PM-2A: shared-field completions from an independently specified source',
        protocol='protocol-pm2a.md',
        stage='A, the source audit; stages B-D follow once these gates hold',
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
        quadrature_diagnostic=quadrature_diagnostic(sources),
        gates_passed=bool(mass_gate['passed'] and row_gate['passed'] and int_gate['passed']
                          and tail['passed']),
        input_sha256={'protocol-pm2a.md': hashlib.sha256((HERE/'protocol-pm2a.md').read_bytes()).hexdigest()},
        checks_short_run=dict(mass_worst=mass_gate['worst_relative_error'],
                              rows_worst=row_gate['worst_relative_change'],
                              rows_control_smallest=row_gate['positive_control']['smallest_relative_change'],
                              integral_worst=int_gate['worst_relative_difference'],
                              tail_worst=tail['worst_ratio_difference'],
                              extrapolated_total_median=float(np.median(frac_tot))),
        runtime_seconds=round(time.time() - t0, 1))
    text = json.dumps(result, indent=1, default=float) + '\n'
    print(text[:2200])
    status = evidence_io.finish(args, 'path-memory-pm2a', text, HERE/'pm2a-results.json',
                                ignore={'/runtime_seconds'})
    return status if result['gates_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

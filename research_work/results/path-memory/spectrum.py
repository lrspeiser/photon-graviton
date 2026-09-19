"""RUT-1 stage 7, part E: the cold-ring spectrum certified down to the axis (protocol-rut7.md).

Stage 6 located roots with Beyn's method on a rectangle whose left edge sat at Re s = 0.001, and called the
result complete. It was complete only above that edge: the two-stage ring has unstable roots at every m
from 1 to 8 below it, the moments had located them, and the rectangle filter threw them away. Its
robustness gate compared Newton's fixed points rather than what the contour located, its rank rule dropped
a genuine eigenvalue in the solver's own test, and it filtered and deduplicated before polishing.

What changes here:

* **The strip's left edge lies to the left of the imaginary axis**, at Re s = -0.008, halfway to the
  nearest poles of B_m(s) at -1/tau_keep. Marginal and slowly growing roots are inside the contour, so they
  are located and classified rather than excluded, and no growth rate is too slow to be seen.
* **Two independent answers from one set of contour evaluations.** Beyn's moments locate the eigenvalues;
  the argument principle counts the zeros of det T. The winding number counts every zero whatever its
  residue, which the moment method cannot promise, so a disagreement between them is a finding.
* Moments are shifted and scaled to the contour, the rank is chosen at the largest gap of the singular
  spectrum, and located roots are polished FIRST and tested afterwards.
"""
import numpy as np

import ring_modes as RM


# ---------------------------------------------------------------- T(s) on many points at once
def M_many(ring, m, S):
    """The 2x2 mode matrix at an array of complex s, vectorised over the contour."""
    S = np.asarray(S, complex)
    I2 = np.eye(2)
    base = (S*S)[:, None, None]*I2 - 2*ring.omega*S[:, None, None]*RM.JMAT - (ring.free + ring.A)[None]
    if ring.strength == 0.:
        return base
    Z = ring.Z[m % ring.n]                                        # (p, 2, 2)
    if not ring.lagged:
        return base + (ring.gain*ring.pw[0])*Z[0][None]
    T = 2*np.pi/ring.omega
    W = 0.
    for c, lam in ring.terms:
        z = lam + S[:, None]
        W = W + c/(1. - np.exp(-z*T))*ring.pw[None, :]*np.exp(-z*ring.p[None, :]/ring.omega)/ring.omega
    return base + np.einsum('sp,pab->sab', W, Z)


def det_many(M):
    return M[:, 0, 0]*M[:, 1, 1] - M[:, 0, 1]*M[:, 1, 0]


# ---------------------------------------------------------------- the contour
def rectangle(rect, panel=(.02, .05, .25, .05), order=16):
    """Counter-clockwise boundary of (re_lo, re_hi, im_lo, im_hi): ordered nodes and complex weights.

    `panel` is the panel length on the left, bottom, right and top edges. The left edge runs 0.008 from
    the axis on one side and from the poles on the other, so its panels are short; Gauss-Legendre on a
    panel of half-length h converges like rho^(-2n) with rho = d/h + sqrt((d/h)^2 + 1) for a singularity at
    distance d, which for d = 0.008, h = 0.01 and n = 16 is about 6e-11.
    """
    re_lo, re_hi, im_lo, im_hi = rect
    corners = [complex(re_lo, im_lo), complex(re_hi, im_lo), complex(re_hi, im_hi), complex(re_lo, im_hi)]
    lengths = (panel[1], panel[2], panel[3], panel[0])             # bottom, right, top, left
    xg, wg = np.polynomial.legendre.leggauss(order)
    z, dz = [], []
    for (a, b), ell in zip(zip(corners, corners[1:] + corners[:1]), lengths):
        n = max(int(np.ceil(abs(b - a)/ell)), 1)
        edges = np.linspace(0., 1., n + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            t = .5*(hi - lo)*xg + .5*(lo + hi)
            z.append(a + (b - a)*t)
            dz.append(.5*(hi - lo)*wg*(b - a))
    return np.concatenate(z), np.concatenate(dz)


def winding(det_fn, z, values, max_phase, depth=40):
    """Zeros of an analytic determinant inside the ordered closed contour, by the argument principle.

    The contour is refined between consecutive nodes until no phase increment exceeds `max_phase`, so the
    count cannot be corrupted by an under-resolved sweep -- which is what made stage 5's counter report a
    root at every mode of the neutral control.
    """
    zs = list(z) + [z[0]]
    vs = list(values) + [values[0]]
    total, worst, inserted = 0., 0., 0

    def sweep(za, zb, va, vb, level):
        nonlocal total, worst, inserted
        d = np.angle(vb/va)
        if abs(d) > max_phase and level < depth:
            zm = .5*(za + zb)
            vm = det_fn(np.array([zm]))[0]
            inserted += 1
            sweep(za, zm, va, vm, level + 1)
            sweep(zm, zb, vm, vb, level + 1)
        else:
            total += d
            worst = max(worst, abs(d))

    for k in range(len(zs) - 1):
        sweep(zs[k], zs[k + 1], vs[k], vs[k + 1], 0)
    n = total/(2*np.pi)
    return dict(count=int(round(n)), raw=float(n), remainder=float(abs(n - round(n))),
                worst_phase_step=float(worst), points_inserted=int(inserted))


# ---------------------------------------------------------------- Beyn, scaled, gap-ranked
def located(M_fn, rect, n=2, moments=12, panel=(.02, .05, .25, .05), order=16, seed=0,
            gap_required=1e3, noise=1e-11):
    """Raw eigenvalues inside the rectangle from scaled moments, with the singular spectrum that decided
    the rank, and the winding number from the same evaluations."""
    z, dz = rectangle(rect, panel, order)
    M = M_fn(z)
    centre = complex(.5*(rect[0] + rect[1]), .5*(rect[2] + rect[3]))
    rho = .5*abs(complex(rect[1] - rect[0], rect[3] - rect[2]))
    rng = np.random.default_rng(seed)
    V = rng.standard_normal((n, n)) + 1j*rng.standard_normal((n, n))
    X = np.linalg.solve(M, np.broadcast_to(V, (len(z), n, n)))     # T(z)^-1 V at every node
    zeta = (z - centre)/rho
    K = moments//2
    A = np.empty((moments, n, n), complex)
    pw = np.ones(len(z), complex)
    for p in range(moments):
        A[p] = np.einsum('s,sab->ab', dz*pw, X)/(2j*np.pi*rho)
        pw = pw*zeta
    B0 = np.block([[A[i + j] for j in range(K)] for i in range(K)])
    B1 = np.block([[A[i + j + 1] for j in range(K)] for i in range(K)])
    U, sv, Wh = np.linalg.svd(B0)
    scale = float(np.sum(np.abs(dz)*np.linalg.norm(X, axis=(1, 2)))/(2*np.pi*rho))
    floor = noise*max(scale, 1e-300)
    if sv[0] <= floor:
        k, gap = 0, float('inf')
    else:
        live = int(np.sum(sv > floor))
        ratios = sv[:live]/np.maximum(np.append(sv[1:live], floor), 1e-300)
        k = int(np.argmax(ratios)) + 1
        gap = float(ratios[k - 1])
    if k:
        red = U[:, :k].conj().T @ B1 @ Wh[:k].conj().T @ np.diag(1./sv[:k])
        raw = centre + rho*np.linalg.eigvals(red)
    else:
        raw = np.array([], complex)
    dets = det_many(M)
    return dict(raw=raw, rank=k, gap=gap, saturated=bool(k >= K*n), singular_values=sv, floor=floor,
                nodes=len(z), contour=z, det_on_contour=dets,
                min_det_on_contour=float(np.min(np.abs(dets))), median_det_on_contour=float(np.median(np.abs(dets))))


def polish(det_scalar, s, limit, iters=80, tol=1e-15):
    """Newton on the determinant, refused if it would carry the root further than `limit`."""
    s0 = complex(s)
    s = s0
    for _ in range(iters):
        f = det_scalar(s)
        h = 1e-7*max(1., abs(s))
        df = (det_scalar(s + h) - det_scalar(s - h))/(2*h)
        if df == 0:
            break
        ds = f/df
        if abs(s - ds - s0) > limit:
            break
        s = s - ds
        if abs(ds) <= tol*max(1., abs(s)):
            break
    return s


def solve_generic(M_fn, rect, thresholds, **kw):
    """Locate, polish, THEN test: region, deduplication and residual are all evaluated on polished values.
    `M_fn` maps an array of complex s to the stack of matrices T(s)."""
    n = M_fn(np.array([complex(.5*(rect[0] + rect[1]), .5*(rect[2] + rect[3]))])).shape[-1]
    det_fn = (lambda zz: det_many(M_fn(zz))) if n == 2 else (lambda zz: np.linalg.det(M_fn(zz)))
    det1 = lambda s: det_fn(np.array([s]))[0]
    loc = located(M_fn, rect, n=n, **kw)
    wind = winding(det_fn, loc['contour'], loc['det_on_contour'], thresholds['winding_phase_step'])
    raw = list(loc['raw'])
    roots = []
    for i, s in enumerate(raw):
        others = [abs(s - q) for j, q in enumerate(raw) if j != i]
        to_edge = min(abs(s.real - rect[0]), abs(s.real - rect[1]), abs(s.imag - rect[2]), abs(s.imag - rect[3]))
        limit = max(min([10*thresholds['polish_movement'], max(to_edge, 1e-9)]
                        + [.5*d for d in others if d > 1e-5]), 1e-9)
        sharp = polish(det1, s, limit)
        roots.append(dict(raw=complex(s), s=complex(sharp), moved=float(abs(sharp - s)),
                          residual=float(abs(det1(sharp)))))
    is_in = lambda z: rect[0] < z.real < rect[1] and rect[2] < z.imag < rect[3]
    inside = [r for r in roots if is_in(r['s'])]
    outside = [r for r in roots if not is_in(r['s'])]
    merged = []
    for r in sorted(inside, key=lambda q: (-q['s'].real, q['s'].imag)):
        for q in merged:
            if abs(q['s'] - r['s']) <= 1e-7*max(1., abs(r['s'])):
                q['multiplicity'] += 1
                q['raws'].append(r['raw'])
                q['moved'] = max(q['moved'], r['moved'])
                q['residual'] = max(q['residual'], r['residual'])
                break
        else:
            merged.append(dict(r, multiplicity=1, raws=[r['raw']]))
    for r in merged:
        r['distance_to_edge'] = float(min(abs(r['s'].real - rect[0]), abs(r['s'].real - rect[1]),
                                          abs(r['s'].imag - rect[2]), abs(r['s'].imag - rect[3])))
    count = int(sum(r['multiplicity'] for r in merged))
    return dict(rect=[float(v) for v in rect], roots=merged, count_with_multiplicity=count,
                winding=wind, counts_agree=bool(count == wind['count']),
                rank=loc['rank'], gap=loc['gap'], saturated=loc['saturated'],
                singular_values=[float(v) for v in loc['singular_values']], floor=loc['floor'],
                nodes=loc['nodes'], min_det_on_contour=loc['min_det_on_contour'],
                median_det_on_contour=loc['median_det_on_contour'],
                discarded_outside=[[float(r['s'].real), float(r['s'].imag)] for r in outside])


def solve_rectangle(ring, m, rect, thresholds, **kw):
    return solve_generic(lambda zz: M_many(ring, m, zz), rect, thresholds, **kw)


def partition(omega, im_max, shifted=False):
    """Cuts in Im s. Roots gather near the memory-pole branches at k*Omega, near the epicyclic frequency
    and near zero, so the cuts sit between those: at (k + 1/2) Omega for the base partition and at
    (k + 0.3) Omega for the shifted one, which shares no cut with it."""
    off = .3 if shifted else .5
    k_max = int(np.ceil(im_max/omega)) + 1
    pos = [(k + off)*omega for k in range(k_max) if (k + off)*omega < im_max - .05*omega]
    cuts = sorted(set([-c for c in pos] + pos))
    return [-im_max] + cuts + [im_max]


NEAR_EDGE, EDGE_MOVE = 1e-3, 2e-3            # the protocol's declared edge rule, not a tolerance of the code


def _near_edges(rows):
    """Roots the declared rule is about: inside a sub-rectangle within NEAR_EDGE of one of its edges, or
    just outside an OUTER edge of the strip (an internal cut's outside is its neighbour's inside)."""
    found = []
    last = len(rows) - 1
    for k, row in enumerate(rows):
        lo, hi, bot, top = row['rect']
        for r in row['roots']:
            if r['distance_to_edge'] >= NEAR_EDGE:
                continue
            s = r['s']
            d = dict(bottom=abs(s.imag - bot), top=abs(s.imag - top), left=abs(s.real - lo), right=abs(s.real - hi))
            found.append(dict(k=k, edge=min(d, key=d.get), where='inside', root=[float(s.real), float(s.imag)]))
        for re_, im_ in row['discarded_outside']:
            d = {}
            if bot <= im_ <= top:
                d.update(left=lo - re_ if re_ < lo else np.inf, right=re_ - hi if re_ > hi else np.inf)
            if k == 0 and im_ < bot:
                d['bottom'] = bot - im_
            if k == last and im_ > top:
                d['top'] = im_ - top
            for edge, dist in d.items():
                if 0. <= dist < NEAR_EDGE:
                    found.append(dict(k=k, edge=edge, where='outside', root=[float(re_), float(im_)]))
    return found


def strip(ring, m, thresholds, shifted=False, im_max=None, re=None, **kw):
    """Every root in the strip, rectangle by rectangle, with the DECLARED edge rule: a root within 0.001
    of an edge sends that edge 0.002 further from it, the strip is swept again with the moved edges, and
    the first sweep's counts are kept beside the second."""
    re_lo, re_hi = re if re is not None else thresholds['strip_re']
    im_max = thresholds['strip_im'] if im_max is None else im_max
    cuts = partition(ring.omega, im_max, shifted)
    n = len(cuts) - 1

    def sweep(cuts, sides):
        return [solve_rectangle(ring, m, (sides[k][0], sides[k][1], cuts[k], cuts[k + 1]), thresholds, **kw)
                for k in range(n)]

    sides = [[re_lo, re_hi] for _ in range(n)]
    rows = sweep(cuts, sides)
    near, first_pass, still = _near_edges(rows), None, []
    if near:
        first_pass = [dict(rect=r['rect'], count=r['count_with_multiplicity'], winding=r['winding']['count'])
                      for r in rows]
        cuts = list(cuts)
        for hit in near:
            k, away = hit['k'], (1. if hit['where'] == 'inside' else -1.)     # away from the root
            if hit['edge'] == 'top':
                cuts[k + 1] += away*EDGE_MOVE
            elif hit['edge'] == 'bottom':
                cuts[k] -= away*EDGE_MOVE
            elif hit['edge'] == 'left':
                sides[k][0] -= away*EDGE_MOVE
            else:
                sides[k][1] += away*EDGE_MOVE
        rows = sweep(cuts, sides)
        still = _near_edges(rows)
    roots = [dict(r, rect=row['rect']) for row in rows for r in row['roots']]
    return dict(m=int(m), rectangles=rows, roots=roots, reruns=near, first_pass=first_pass,
                still_near_an_edge=still,
                counts_agree=all(r['counts_agree'] for r in rows),
                any_saturated=any(r['saturated'] for r in rows),
                smallest_gap=float(min([r['gap'] for r in rows if r['rank']] or [float('inf')])),
                worst_residual=float(max([q['residual'] for q in roots] or [0.])),
                worst_polish_movement=float(max([q['moved'] for q in roots] or [0.])),
                worst_winding_remainder=float(max(r['winding']['remainder'] for r in rows)),
                worst_winding_phase_step=float(max(r['winding']['worst_phase_step'] for r in rows)),
                smallest_det_on_contour_over_median=float(min(r['min_det_on_contour']/max(r['median_det_on_contour'], 1e-300)
                                                              for r in rows)))


def classify(roots, omega, m, floor):
    out = []
    for r in roots:
        s = r['s']
        kind = 'unstable' if s.real > floor else ('marginal' if abs(s.real) <= floor else 'stable')
        out.append(dict(re=float(s.real), im=float(s.imag), multiplicity=int(r['multiplicity']), kind=kind,
                        residual=r['residual'], moved_by_polish=r['moved'],
                        raw=[[float(z.real), float(z.imag)] for z in r['raws']],
                        e_folding_periods=float(1./(s.real*RM.T0)) if s.real > floor else None,
                        pattern_speed_over_omega=float((omega - s.imag/m)/omega) if m else None))
    return out


def continuation_check(ring, m, certified, strip_re, strip_im, tol=1e-6):
    """Stage 5's continuation as a cross-check only: distinctness, the lost partner by deflation, and what
    the certified list holds that continuation never reaches."""
    roots, _ = ring.modes(m)
    roots = [complex(z) for z in roots]
    dup = [(i, j) for i in range(len(roots)) for j in range(i + 1, len(roots))
           if abs(roots[i] - roots[j]) < 1e-9*max(1., abs(roots[i]))]
    recovered = []
    det1 = lambda s: ring.det(m, s)
    for i, j in dup:
        known = [roots[k] for k in range(len(roots)) if k != j]
        deflated = lambda s: det1(s)/np.prod([s - q for q in known])
        guess = -roots[i]
        s = polish(deflated, guess, limit=.05)
        if abs(det1(s)) < 1e-9 and all(abs(s - q) > 1e-8 for q in known):
            recovered.append(s)
    distinct = []
    for z in roots + recovered:
        if all(abs(z - q) > 1e-9 for q in distinct):
            distinct.append(z)
    inside = [z for z in distinct if strip_re[0] < z.real < strip_re[1] and abs(z.imag) < strip_im]
    cert = [complex(r['re'], r['im']) for r in certified]
    missing = [c for c in cert if all(abs(c - z) > tol*max(1., abs(c)) for z in inside)]
    extra = [z for z in inside if all(abs(c - z) > tol*max(1., abs(c)) for c in cert)]
    return dict(duplicated_branches=len(dup), recovered_by_deflation=[[float(z.real), float(z.imag)] for z in recovered],
                certified_roots_continuation_lacks=[[float(z.real), float(z.imag)] for z in missing],
                continuation_roots_not_certified=[[float(z.real), float(z.imag)] for z in extra])

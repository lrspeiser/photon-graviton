"""Contour-integral eigensolver for the nonlinear eigenproblem T(s)v = 0 (protocol-rut6.md, R1).

Stage 5 counted roots by the winding of det T(s) around a contour. That count certifies nothing and is
itself wrong: on the no-attraction control, whose determinant is exactly s^2(s^2 + Omega^2) and has no
root in the open right half-plane, it returns one unstable root at every angular mode, because the double
root at the origin sits on the contour edge with its phase winding under-resolved.

This module locates roots instead of counting them, by the method Beyn introduced (arXiv:1003.1580). For
a probe matrix V and a contour enclosing a region where T is analytic,

    A_p = (1/2 pi i) contour_integral s^p T(s)^(-1) V ds,     p = 0 .. 2K-1,

and the eigenvalues inside follow from the block Hankel matrices of those moments: the SVD of the first
gives the numerical rank, and projecting the second onto that subspace gives a small ordinary eigenvalue
problem whose eigenvalues are exactly the eigenvalues of T inside the contour, with multiplicity. No
initial guess is needed and nothing is followed from one parameter value to the next, which is what a
statement about ALL roots in a region requires.

Two properties matter for how it is used here. The contour must lie where T is analytic -- for the ring
problem B_m(s) has poles at s = -1/tau + i k Omega, so the search region is kept in Re(s) > 0 -- and a
root sitting on or very near the contour is resolved slowly, so every located root is reported with its
distance to the contour and the caller is expected to vary the contour and the resolution.
"""
import numpy as np


def rectangle(re_lo, re_hi, im_lo, im_hi, panels=12, order=16):
    """Counter-clockwise rectangle boundary: nodes and complex weights, Gauss-Legendre per panel."""
    corners = [complex(re_lo, im_lo), complex(re_hi, im_lo), complex(re_hi, im_hi), complex(re_lo, im_hi)]
    xg, wg = np.polynomial.legendre.leggauss(order)
    z, dz = [], []
    for a, b in zip(corners, corners[1:] + corners[:1]):
        edges = np.linspace(0., 1., panels + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            t = .5*(hi - lo)*xg + .5*(lo + hi)
            z.append(a + (b - a)*t)
            dz.append(.5*(hi - lo)*wg*(b - a))
    return np.concatenate(z), np.concatenate(dz)


def circle(centre, radius, points=256):
    """Counter-clockwise circle, trapezoid rule -- spectrally accurate for an analytic integrand."""
    th = 2*np.pi*np.arange(points)/points
    z = centre + radius*np.exp(1j*th)
    return z, radius*1j*np.exp(1j*th)*(2*np.pi/points)


def _inside(z, re_lo, re_hi, im_lo, im_hi, margin=0.):
    return ((z.real > re_lo + margin) & (z.real < re_hi - margin)
            & (z.imag > im_lo + margin) & (z.imag < im_hi - margin))


def eigenvalues(T, nodes, weights, n, probes=None, moments=8, rank_tol=1e-8, seed=0):
    """Eigenvalues of T inside the contour, with the singular values that determined the rank.

    `moments` is 2K: K block rows of Hankel structure, so up to K*probes eigenvalues can be resolved.
    A returned rank equal to that capacity means the region may hold more, and the caller is told.
    """
    rng = np.random.default_rng(seed)
    ell = n if probes is None else int(probes)
    V = rng.standard_normal((n, ell)) + 1j*rng.standard_normal((n, ell))
    K = moments//2
    A = np.zeros((moments, n, ell), complex)
    for z, w in zip(nodes, weights):
        R = np.linalg.solve(T(z), V)
        zp = 1.
        for p in range(moments):
            A[p] += w*zp*R
            zp *= z
    A /= 2j*np.pi
    B0 = np.block([[A[i + j] for j in range(K)] for i in range(K)])
    B1 = np.block([[A[i + j + 1] for j in range(K)] for i in range(K)])
    U, sv, Wh = np.linalg.svd(B0)
    cutoff = rank_tol*max(sv[0], 1e-300)
    k = int(np.sum(sv > cutoff))
    if k == 0:
        return np.array([], complex), sv, False
    M = U[:, :k].conj().T @ B1 @ Wh[:k].conj().T @ np.diag(1./sv[:k])
    return np.linalg.eigvals(M), sv, bool(k >= K*ell)


def newton_polish(det, s, iters=60, tol=1e-14):
    """The moments locate a root to the accuracy of the contour quadrature; this sharpens it. It cannot
    move a root far, so it refines what the contour found rather than searching."""
    s0 = complex(s)
    for _ in range(iters):
        f = det(s)
        h = 1e-7*max(1., abs(s))
        df = (det(s + h) - det(s - h))/(2*h)
        if df == 0:
            break
        ds = f/df
        if abs(ds) > .05*max(1., abs(s0)):          # refuse to wander off the located root
            break
        s = s - ds
        if abs(ds) <= tol*max(1., abs(s)):
            break
    return s if abs(s - s0) <= .05*max(1., abs(s0)) else s0


def solve_region(T, region, n=2, panels=12, order=16, moments=8, rank_tol=1e-8,
                 residual=None, dedup=1e-7, margin=1e-9, polish=None):
    """Locate every eigenvalue of T inside a declared rectangle.

    Returns the distinct roots inside, each with its multiplicity, its distance to the contour and, if a
    residual function is supplied, its residual. `saturated` warns that the moment capacity was reached,
    which means the region may hold roots this call could not resolve.
    """
    re_lo, re_hi, im_lo, im_hi = region
    z, dz = rectangle(re_lo, re_hi, im_lo, im_hi, panels, order)
    raw, sv, saturated = eigenvalues(T, z, dz, n, moments=moments, rank_tol=rank_tol)
    keep = raw[_inside(raw, re_lo, re_hi, im_lo, im_hi, margin)] if len(raw) else raw
    roots = []
    for s in sorted(keep, key=lambda v: (-v.real, v.imag)):
        for r in roots:
            if abs(r['s'] - s) <= dedup*max(1., abs(s)):
                r['multiplicity'] += 1
                break
        else:
            roots.append(dict(s=complex(s), multiplicity=1))
    for r in roots:
        if polish is not None:
            sharp = newton_polish(polish, r['s'])
            r['moved_by_polish'] = float(abs(sharp - r['s']))
            r['s'] = sharp
        s = r['s']
        r['distance_to_contour'] = float(min(abs(s.real - re_lo), abs(s.real - re_hi),
                                             abs(s.imag - im_lo), abs(s.imag - im_hi)))
        if residual is not None:
            r['residual'] = float(abs(residual(s)))
    return dict(region=[float(x) for x in region], roots=roots, count=len(roots),
                total_with_multiplicity=int(sum(r['multiplicity'] for r in roots)),
                singular_values=[float(x) for x in sv[:12]], saturated=saturated,
                nodes=int(len(z)), discarded_outside=int(len(raw) - len(keep)))


def in_region(z, region, margin=0.):
    """Is z inside the declared rectangle? Comparisons must filter both lists the same way, or a root the
    solver correctly did not seek is counted as one it missed."""
    re_lo, re_hi, im_lo, im_hi = region
    z = complex(z)
    return (re_lo + margin < z.real < re_hi - margin) and (im_lo + margin < z.imag < im_hi - margin)


def agree(a, b, tol=1e-6):
    """Match two root lists; returns the worst pairing distance and anything unmatched either way."""
    a = [complex(x) for x in a]
    b = [complex(x) for x in b]
    unmatched_a, worst, used = [], 0., set()
    for x in a:
        best, which = None, None
        for i, y in enumerate(b):
            if i in used:
                continue
            d = abs(x - y)/max(1., abs(x))
            if best is None or d < best:
                best, which = d, i
        if best is not None and best <= tol:
            used.add(which)
            worst = max(worst, best)
        else:
            unmatched_a.append([x.real, x.imag])
    unmatched_b = [[b[i].real, b[i].imag] for i in range(len(b)) if i not in used]
    return dict(worst_relative=float(worst), unmatched_first=unmatched_a, unmatched_second=unmatched_b,
                matched=len(used), agreed=bool(not unmatched_a and not unmatched_b))

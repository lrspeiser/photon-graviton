"""CC-2 stage 2A engine (see protocol.md): tracer companions on exact orbits in an evolving spherical potential, an
incoming bath in the transparent (Liouville) limit with first-order depletion, and the three collision classes of
CF-1's law (elastic, equal-mass, isotropic in the centre-of-mass frame, constant sigma/m).

Units: kpc, km/s, Msun; time in kpc/(km/s) (0.9778 Gyr). Potentials are referenced to Phi(R_b) = 0 as in CF-1. A
companion is "confined" when it cannot reach R_b in the current potential: E < 0, or E >= 0 behind the centrifugal
barrier.

Tracers carry individual masses. Incoming-incoming production is volume-dominated, so equal masses would leave the
inner region unresolved; tracers born in production shell j get a mass proportional to the square root of that
shell's initial production rate. Every collision still pairs equal masses: in bound-bound collisions a heavier
tracer is split and only the matching part collides, and the rest carries on to the next radial neighbour, so each
unit of mass is tried exactly once per step and energy and momentum are conserved exactly.

Kernels are compiled with numba without an on-disk cache, so nothing is written next to this file.
"""
import math
import time
import numpy as np
import numba as nb
from scipy.integrate import cumulative_trapezoid
from scipy.ndimage import maximum_filter1d
from scipy.optimize import brentq
from scipy.special import erf

G = 4.30091727003628e-6
PER_GYR = 1.0227121650537077            # 1 (km/s)/kpc in 1/Gyr
CM2_PER_G = .1*1.98847e30/3.0856775814913673e19**2


# ----------------------------------------------------------------------------------------------------------------
# incident distributions at R_b
def single_speed(u):
    return np.array([float(u)]), np.array([1.])


def maxwell_components(s, K=64):
    """Isotropic Maxwellian of 1-D dispersion s as K equal-probability speed bins, each at its conditional mean."""
    cdf = lambda U: erf(U/(math.sqrt(2)*s)) - math.sqrt(2/math.pi)*(U/s)*math.exp(-U*U/(2*s*s))
    m1 = lambda U: 2*math.sqrt(2/math.pi)*s*(1 - math.exp(-U*U/(2*s*s))*(1 + U*U/(2*s*s)))
    edges = [0.]
    for k in range(1, K):
        edges.append(brentq(lambda U: cdf(U) - k/K, 0., 20*s, xtol=1e-12*s))
    edges.append(60*s)
    u = np.array([(m1(edges[k + 1]) - m1(edges[k]))*K for k in range(K)])
    return u, np.full(K, 1./K)


def mixture(parts):
    """Weighted mixture of component lists [(weight, (u, f)), ...]."""
    u = np.concatenate([p[1][0] for p in parts])
    f = np.concatenate([p[0]*p[1][1] for p in parts])
    return u, f/f.sum()


# ----------------------------------------------------------------------------------------------------------------
# compiled kernels
@nb.njit(cache=False)
def _pos(lr0, dl, n, r):
    """Grid cell and fraction for radius r. Checked before any integer conversion: int() of an infinite or NaN
    value is undefined in compiled code and would index outside the tables."""
    x = (math.log(r) - lr0)/dl
    if not (x > 0.):
        return 0, 0.
    if x >= n - 1:
        return n - 2, 1.
    i = int(x)
    return i, x - i


@nb.njit(cache=False)
def _phi(lr0, dl, n, phitab, gtab, r_lo, r):
    if r < r_lo:
        return phitab[0] - .5*gtab[0]*r_lo*(1. - (r/r_lo)**2)
    i, f = _pos(lr0, dl, n, r)
    return phitab[i]*(1. - f) + phitab[i + 1]*f


@nb.njit(cache=False)
def _gacc(lr0, dl, n, gtab, r_lo, r):
    if r < r_lo:
        return gtab[0]*r/r_lo
    i, f = _pos(lr0, dl, n, r)
    return gtab[i]*(1. - f) + gtab[i + 1]*f


@nb.njit(cache=False)
def _confined(lr0, dl, n, phitab, r_lo, r, E, J2):
    """True if a companion at radius r with specific energy E and J^2 cannot reach R_b (the grid's last node)."""
    if E < 0.:
        return True
    x = (math.log(max(r, r_lo)) - lr0)/dl
    if not (x < n - 1):                        # at or beyond R_b (or not finite): nothing left to cross
        return False
    i0 = int(math.ceil(x))
    if i0 < 0:
        i0 = 0
    rr2 = math.exp(2.*(lr0 + i0*dl))          # r^2 at node i0, stepped geometrically along the log grid
    q = math.exp(2.*dl)
    for i in range(i0, n):
        if 2.*(E - phitab[i])*rr2 < J2:
            return True
        rr2 *= q
    return False


@nb.njit(cache=False)
def confined_many(r, E, J2, lr0, dl, n, phitab, r_lo):
    out = np.empty(r.shape[0], np.bool_)
    for j in range(r.shape[0]):
        out[j] = _confined(lr0, dl, n, phitab, r_lo, r[j], E[j], J2[j])
    return out


@nb.njit(cache=False)
def _unit():
    z = 2.*np.random.random() - 1.
    p = 2.*math.pi*np.random.random()
    s = math.sqrt(max(0., 1. - z*z))
    return s*math.cos(p), s*math.sin(p), z


@nb.njit(cache=False)
def _frame(x0, x1, x2):
    """Radial unit vector and two tangential unit vectors at position x."""
    r = math.sqrt(x0*x0 + x1*x1 + x2*x2)
    er0, er1, er2 = x0/r, x1/r, x2/r
    if abs(er2) < .9:
        a0, a1, a2 = 0., 0., 1.
    else:
        a0, a1, a2 = 1., 0., 0.
    e10, e11, e12 = er1*a2 - er2*a1, er2*a0 - er0*a2, er0*a1 - er1*a0
    nn = math.sqrt(e10*e10 + e11*e11 + e12*e12)
    e10, e11, e12 = e10/nn, e11/nn, e12/nn
    e20, e21, e22 = er1*e12 - er2*e11, er2*e10 - er0*e12, er0*e11 - er1*e10
    return er0, er1, er2, e10, e11, e12, e20, e21, e22


@nb.njit(cache=False)
def _bath_density(lr0, dl, n, dtab, r):
    """Sum over components of the local density factor (per unit incident density), and the grid position."""
    i, f = _pos(lr0, dl, n, max(r, math.exp(lr0)))
    tot = 0.
    for k in range(dtab.shape[0]):
        tot += dtab[k, i]*(1. - f) + dtab[k, i + 1]*f
    return tot, i, f


@nb.njit(cache=False)
def _sample_bath(dtab, stab, uk, phi, tot, i, f):
    """Draw one incoming companion at a radius with potential phi: component by local density, speed w, direction
    uniform over the allowed cones |cos theta| >= sqrt(1 - s_max^2). Returns (k, v_r, v_t1, v_t2) locally."""
    y = np.random.random()*tot
    acc = 0.
    k = dtab.shape[0] - 1
    for kk in range(dtab.shape[0]):
        acc += dtab[kk, i]*(1. - f) + dtab[kk, i + 1]*f
        if y < acc:
            k = kk
            break
    w = math.sqrt(max(uk[k]*uk[k] - 2.*phi, 0.))
    s = min(1., stab[k, i]*(1. - f) + stab[k, i + 1]*f)
    cm = math.sqrt(max(0., 1. - s*s))
    c = cm + (1. - cm)*np.random.random()
    if np.random.random() < .5:
        c = -c
    st = math.sqrt(max(0., 1. - c*c))
    p = 2.*math.pi*np.random.random()
    return k, w*c, w*st*math.cos(p), w*st*math.sin(p)


@nb.njit(cache=False)
def sample_bath_many(m, r, lr0, dl, n, dtab, stab, uk, phitab, gtab, r_lo):
    """V2: m draws at one radius; returns local-frame velocities."""
    tot, i, f = _bath_density(lr0, dl, n, dtab, r)
    phi = _phi(lr0, dl, n, phitab, gtab, r_lo, r)
    out = np.empty((m, 3))
    for j in range(m):
        k, a, b, c = _sample_bath(dtab, stab, uk, phi, tot, i, f)
        out[j, 0], out[j, 1], out[j, 2] = a, b, c
    return out


@nb.njit(cache=False)
def _peri(lr0, dl, n, phitab, gtab, r_lo, E, J2, r):
    """Pericentre by bisection in log r below the current radius (precision about 0.1%)."""
    a = math.log(r_lo*1e-4)
    b = math.log(r)
    for _ in range(22):
        m = .5*(a + b)
        rm = math.exp(m)
        if 2.*(E - _phi(lr0, dl, n, phitab, gtab, r_lo, rm))*rm*rm - J2 < 0.:
            a = m
        else:
            b = m
    return math.exp(b)


@nb.njit(cache=False)
def advance(x, v, alive, Delta, lr0, dl, n, gtab, phitab, r_lo, R_b, eta, r_soft, max_sub, e_out):
    """Leapfrog (KDK) for every live tracer through Delta with its own fixed substep, set from its pericentre
    speed. Tracers reaching R_b are marked dead with their specific energy at crossing in e_out.
    Returns the substep count and the number of tracers whose substeps were capped."""
    tot = 0
    capped = 0
    for i in range(x.shape[0]):
        if not alive[i]:
            continue
        x0, x1, x2 = x[i, 0], x[i, 1], x[i, 2]
        v0, v1, v2 = v[i, 0], v[i, 1], v[i, 2]
        r = math.sqrt(x0*x0 + x1*x1 + x2*x2)
        E = .5*(v0*v0 + v1*v1 + v2*v2) + _phi(lr0, dl, n, phitab, gtab, r_lo, r)
        j0, j1, j2 = x1*v2 - x2*v1, x2*v0 - x0*v2, x0*v1 - x1*v0
        J2 = j0*j0 + j1*j1 + j2*j2
        rp = _peri(lr0, dl, n, phitab, gtab, r_lo, E, J2, r)
        vp = math.sqrt(max(2.*(E - _phi(lr0, dl, n, phitab, gtab, r_lo, rp)), 1e-30))
        h = eta*(rp + r_soft)/vp
        k = int(math.ceil(Delta/h))
        if k > max_sub:
            k = max_sub
            capped += 1
        if k < 1:
            k = 1
        h = Delta/k
        ga = _gacc(lr0, dl, n, gtab, r_lo, r)/r
        a0, a1, a2 = -ga*x0, -ga*x1, -ga*x2
        for s in range(k):
            v0 += .5*h*a0
            v1 += .5*h*a1
            v2 += .5*h*a2
            x0 += h*v0
            x1 += h*v1
            x2 += h*v2
            r = math.sqrt(x0*x0 + x1*x1 + x2*x2)
            if r >= R_b:
                alive[i] = False
                e_out[i] = .5*(v0*v0 + v1*v1 + v2*v2)
                break
            ga = _gacc(lr0, dl, n, gtab, r_lo, r)/r
            a0, a1, a2 = -ga*x0, -ga*x1, -ga*x2
            v0 += .5*h*a0
            v1 += .5*h*a1
            v2 += .5*h*a2
        tot += k
        x[i, 0], x[i, 1], x[i, 2] = x0, x1, x2
        v[i, 0], v[i, 1], v[i, 2] = v0, v1, v2
    return tot, capped


@nb.njit(cache=False)
def energies(x, v, lr0, dl, n, phitab, gtab, r_lo):
    out = np.empty(x.shape[0])
    for i in range(x.shape[0]):
        r = math.sqrt(x[i, 0]**2 + x[i, 1]**2 + x[i, 2]**2)
        out[i] = .5*(v[i, 0]**2 + v[i, 1]**2 + v[i, 2]**2) + _phi(lr0, dl, n, phitab, gtab, r_lo, r)
    return out


@nb.njit(cache=False)
def _collide(a0, a1, a2, b0, b1, b2):
    c0, c1, c2 = .5*(a0 + b0), .5*(a1 + b1), .5*(a2 + b2)
    g = math.sqrt((a0 - b0)**2 + (a1 - b1)**2 + (a2 - b2)**2)
    n0, n1, n2 = _unit()
    h = .5*g
    return c0 + h*n0, c1 + h*n1, c2 + h*n2, c0 - h*n0, c1 - h*n1, c2 - h*n2


@nb.njit(cache=False)
def bath_bound(x, v, m, alive, nsub, Delta, sm, rho_inf, lr0, dl, n, dtab, stab, uk, phitab, gtab, r_lo, bound_only,
               spawn_x, spawn_v, spawn_m, e_rem, led, pstat):
    """Incoming-bound collisions. Tracer i takes nsub[i] sub-steps of Delta/nsub[i] at its fixed position; in each it
    meets an incoming companion drawn from the local distribution with probability sm*rho_bath*g*(Delta/nsub[i]).
    The incoming partner carries the tracer's mass. A confined partner becomes a new tracer (spawn_*); a tracer that
    ends unconfined is marked dead with its energy in e_rem. bound_only uses E < 0 in place of confinement (CF-1).
    led: [collisions, captures, ejections, dE_confined, E_bath_in, E_exported, max rel. kinetic error,
          max rel. momentum error, captured mass, ejected mass, net count]; energies are mass-weighted.
    pstat: [max P, count of P > 1]. Returns the number of spawns."""
    ns = 0
    for i in range(x.shape[0]):
        if not alive[i]:
            continue
        x0, x1, x2 = x[i, 0], x[i, 1], x[i, 2]
        r = math.sqrt(x0*x0 + x1*x1 + x2*x2)
        tot, gi, gf = _bath_density(lr0, dl, n, dtab, r)
        if tot <= 0.:
            continue
        phi = _phi(lr0, dl, n, phitab, gtab, r_lo, r)
        er0, er1, er2, e10, e11, e12, e20, e21, e22 = _frame(x0, x1, x2)
        d = Delta/nsub[i]
        mi = m[i]
        for s in range(nsub[i]):
            k, br, b1, b2 = _sample_bath(dtab, stab, uk, phi, tot, gi, gf)
            w0 = br*er0 + b1*e10 + b2*e20
            w1 = br*er1 + b1*e11 + b2*e21
            w2 = br*er2 + b1*e12 + b2*e22
            v0, v1, v2 = v[i, 0], v[i, 1], v[i, 2]
            g = math.sqrt((v0 - w0)**2 + (v1 - w1)**2 + (v2 - w2)**2)
            P = sm*rho_inf*tot*g*d
            if P > pstat[0]:
                pstat[0] = P
            if P > 1.:
                pstat[1] += 1.
            if np.random.random() >= P:
                continue
            K0 = v0*v0 + v1*v1 + v2*v2 + w0*w0 + w1*w1 + w2*w2
            p0, p1, p2, q0, q1, q2 = _collide(v0, v1, v2, w0, w1, w2)
            K1 = p0*p0 + p1*p1 + p2*p2 + q0*q0 + q1*q1 + q2*q2
            led[6] = max(led[6], abs(K1 - K0)/K0)
            dp = math.sqrt((p0 + q0 - v0 - w0)**2 + (p1 + q1 - v1 - w1)**2 + (p2 + q2 - v2 - w2)**2)
            led[7] = max(led[7], dp/math.sqrt(K0))
            led[0] += 1.
            E_old = .5*(v0*v0 + v1*v1 + v2*v2) + phi
            E_t = .5*(p0*p0 + p1*p1 + p2*p2) + phi
            E_b = .5*(q0*q0 + q1*q1 + q2*q2) + phi
            led[4] += mi*.5*uk[k]*uk[k]
            if bound_only:
                ct = E_t < 0.
                cb = E_b < 0.
            else:
                jt0, jt1, jt2 = x1*p2 - x2*p1, x2*p0 - x0*p2, x0*p1 - x1*p0
                jb0, jb1, jb2 = x1*q2 - x2*q1, x2*q0 - x0*q2, x0*q1 - x1*q0
                ct = _confined(lr0, dl, n, phitab, r_lo, r, E_t, jt0*jt0 + jt1*jt1 + jt2*jt2)
                cb = _confined(lr0, dl, n, phitab, r_lo, r, E_b, jb0*jb0 + jb1*jb1 + jb2*jb2)
            v[i, 0], v[i, 1], v[i, 2] = p0, p1, p2
            if cb and ns >= spawn_x.shape[0]:
                led[11] += 1.                   # buffer full: capture faster than the caller allowed (fail fast)
                return ns
            if cb:
                spawn_x[ns, 0], spawn_x[ns, 1], spawn_x[ns, 2] = x0, x1, x2
                spawn_v[ns, 0], spawn_v[ns, 1], spawn_v[ns, 2] = q0, q1, q2
                spawn_m[ns] = mi
                ns += 1
                led[1] += 1.
                led[3] += mi*E_b
                led[8] += mi
                led[10] += 1.
            else:
                led[5] += mi*E_b
            if ct:
                led[3] += mi*(E_t - E_old)
            else:
                alive[i] = False
                e_rem[i] = E_t
                led[2] += 1.
                led[3] -= mi*E_old
                led[5] += mi*E_t
                led[9] += mi
                led[10] -= 1.
                break
    return ns


@nb.njit(cache=False)
def bound_bound(x, v, m, alive, order, rho_c, Delta, sm, lr0, dl, n, phitab, gtab, r_lo,
                new_x, new_v, new_m, ex_E, ex_m, led, pstat):
    """Bound-bound collisions along the radial order of live tracers. A carried mass (the untried part of the
    current tracer) is paired with the next tracer; the matching mass mu = min(carry, next) is tried once with
    probability sm*rho_c*g*Delta. On a collision a heavier tracer is split: its remainder keeps the old velocity and
    the colliding part becomes a new tracer (new_*). The partner's velocity is carried into the first tracer's local
    frame with a random rotation about the radial direction, which spherical symmetry allows.
    Evaporated parts are written to ex_E/ex_m. led: [collisions, evaporations, dE_confined, E_exported,
    max rel. kinetic error, evaporated mass, splits, longest carry chain, exits written]."""
    N = order.shape[0]
    nn = 0
    if N < 2:
        return 0
    ci = order[0]
    cmass = m[ci]
    chain = 0
    j = 1
    tol = 1e-12
    while j < N:
        b = order[j]
        mb = m[b]
        mu = min(cmass, mb)
        xa0, xa1, xa2 = x[ci, 0], x[ci, 1], x[ci, 2]
        xb0, xb1, xb2 = x[b, 0], x[b, 1], x[b, 2]
        ar0, ar1, ar2, a10, a11, a12, a20, a21, a22 = _frame(xa0, xa1, xa2)
        br0, br1, br2, b10, b11, b12, b20, b21, b22 = _frame(xb0, xb1, xb2)
        va0, va1, va2 = v[ci, 0], v[ci, 1], v[ci, 2]
        vb0, vb1, vb2 = v[b, 0], v[b, 1], v[b, 2]
        u0 = va0*ar0 + va1*ar1 + va2*ar2
        u1 = va0*a10 + va1*a11 + va2*a12
        u2 = va0*a20 + va1*a21 + va2*a22
        s0 = vb0*br0 + vb1*br1 + vb2*br2
        t1 = vb0*b10 + vb1*b11 + vb2*b12
        t2 = vb0*b20 + vb1*b21 + vb2*b22
        ps = 2.*math.pi*np.random.random()
        cp, sp = math.cos(ps), math.sin(ps)
        s1, s2 = t1*cp - t2*sp, t1*sp + t2*cp
        g = math.sqrt((u0 - s0)**2 + (u1 - s1)**2 + (u2 - s2)**2)
        P = sm*.5*(rho_c[ci] + rho_c[b])*g*Delta
        if P > pstat[0]:
            pstat[0] = P
        if P > 1.:
            pstat[1] += 1.
        if np.random.random() < P:
            ra = math.sqrt(xa0*xa0 + xa1*xa1 + xa2*xa2)
            rb = math.sqrt(xb0*xb0 + xb1*xb1 + xb2*xb2)
            pa = _phi(lr0, dl, n, phitab, gtab, r_lo, ra)
            pb = _phi(lr0, dl, n, phitab, gtab, r_lo, rb)
            Ea = .5*(u0*u0 + u1*u1 + u2*u2) + pa
            Eb = .5*(s0*s0 + s1*s1 + s2*s2) + pb
            K0 = u0*u0 + u1*u1 + u2*u2 + s0*s0 + s1*s1 + s2*s2
            p0, p1, p2, q0, q1, q2 = _collide(u0, u1, u2, s0, s1, s2)
            K1 = p0*p0 + p1*p1 + p2*p2 + q0*q0 + q1*q1 + q2*q2
            led[4] = max(led[4], abs(K1 - K0)/K0)
            led[0] += 1.
            na0 = p0*ar0 + p1*a10 + p2*a20
            na1 = p0*ar1 + p1*a11 + p2*a21
            na2 = p0*ar2 + p1*a12 + p2*a22
            c1, c2 = q1*cp + q2*sp, -q1*sp + q2*cp
            nb0 = q0*br0 + c1*b10 + c2*b20
            nb1 = q0*br1 + c1*b11 + c2*b21
            nb2 = q0*br2 + c1*b12 + c2*b22
            Ea2 = .5*(p0*p0 + p1*p1 + p2*p2) + pa
            Eb2 = .5*(q0*q0 + q1*q1 + q2*q2) + pb
            ja0, ja1, ja2 = xa1*na2 - xa2*na1, xa2*na0 - xa0*na2, xa0*na1 - xa1*na0
            jb0, jb1, jb2 = xb1*nb2 - xb2*nb1, xb2*nb0 - xb0*nb2, xb0*nb1 - xb1*nb0
            conf_a = _confined(lr0, dl, n, phitab, r_lo, ra, Ea2, ja0*ja0 + ja1*ja1 + ja2*ja2)
            conf_b = _confined(lr0, dl, n, phitab, r_lo, rb, Eb2, jb0*jb0 + jb1*jb1 + jb2*jb2)
            # the colliding part of each side
            for side in range(2):
                if side == 0:
                    t, E1, E2, cf = ci, Ea, Ea2, conf_a
                    y0, y1, y2, w0, w1, w2 = xa0, xa1, xa2, na0, na1, na2
                else:
                    t, E1, E2, cf = b, Eb, Eb2, conf_b
                    y0, y1, y2, w0, w1, w2 = xb0, xb1, xb2, nb0, nb1, nb2
                split = m[t] > mu*(1. + tol)
                if split:
                    m[t] -= mu
                    led[6] += 1.
                    if cf:
                        new_x[nn, 0], new_x[nn, 1], new_x[nn, 2] = y0, y1, y2
                        new_v[nn, 0], new_v[nn, 1], new_v[nn, 2] = w0, w1, w2
                        new_m[nn] = mu
                        nn += 1
                else:
                    v[t, 0], v[t, 1], v[t, 2] = w0, w1, w2
                    if not cf:
                        alive[t] = False
                if cf:
                    led[2] += mu*(E2 - E1)
                else:
                    led[1] += 1.
                    led[2] -= mu*E1
                    led[3] += mu*E2
                    led[5] += mu
                    ne = int(led[8])
                    ex_E[ne] = E2
                    ex_m[ne] = mu
                    led[8] += 1.
        # carry
        if cmass < mb*(1. - tol):
            ci = b
            cmass = mb - mu
            j += 1
            chain = 0
        elif cmass > mb*(1. + tol):
            cmass = cmass - mu
            j += 1
            chain += 1
            if chain > led[7]:
                led[7] = chain
        else:
            if j + 1 < N:
                ci = order[j + 1]
                cmass = m[ci]
            j += 2
            chain = 0
    return nn


@nb.njit(cache=False)
def seedless_events(r_s, lr0, dl, n, dtab, stab, uk, phitab, gtab, r_lo, out):
    """Sample incoming-incoming collisions at the radii r_s (drawn with the desired weights). For each:
    out[j] = [r, g, n_conf, E_in, conf1, v1r, v1t1, v1t2, E1, conf2, v2r, v2t1, v2t2, E2, n_bound(E<0)]."""
    for j in range(r_s.shape[0]):
        r = r_s[j]
        tot, gi, gf = _bath_density(lr0, dl, n, dtab, r)
        phi = _phi(lr0, dl, n, phitab, gtab, r_lo, r)
        k1, a0, a1, a2 = _sample_bath(dtab, stab, uk, phi, tot, gi, gf)
        k2, b0, b1, b2 = _sample_bath(dtab, stab, uk, phi, tot, gi, gf)
        g = math.sqrt((a0 - b0)**2 + (a1 - b1)**2 + (a2 - b2)**2)
        p0, p1, p2, q0, q1, q2 = _collide(a0, a1, a2, b0, b1, b2)
        E1 = .5*(p0*p0 + p1*p1 + p2*p2) + phi
        E2 = .5*(q0*q0 + q1*q1 + q2*q2) + phi
        c1 = _confined(lr0, dl, n, phitab, r_lo, r, E1, r*r*(p1*p1 + p2*p2))
        c2 = _confined(lr0, dl, n, phitab, r_lo, r, E2, r*r*(q1*q1 + q2*q2))
        f1 = 1. if c1 else 0.
        f2 = 1. if c2 else 0.
        out[j, 0] = r
        out[j, 1] = g
        out[j, 2] = f1 + f2
        out[j, 3] = .5*(uk[k1]*uk[k1] + uk[k2]*uk[k2])
        out[j, 4], out[j, 5], out[j, 6], out[j, 7], out[j, 8] = f1, p0, p1, p2, E1
        out[j, 9], out[j, 10], out[j, 11], out[j, 12], out[j, 13] = f2, q0, q1, q2, E2
        out[j, 14] = (1. if E1 < 0. else 0.) + (1. if E2 < 0. else 0.)
    return 0


@nb.njit(cache=False)
def seed_rng(s):
    np.random.seed(s)


@nb.njit(cache=False)
def place(r, vloc, out_x, out_v):
    """Put tracers at radius r in random directions, with local-frame velocities vloc = (v_r, v_t1, v_t2)."""
    for j in range(r.shape[0]):
        u0, u1, u2 = _unit()
        x0, x1, x2 = r[j]*u0, r[j]*u1, r[j]*u2
        er0, er1, er2, e10, e11, e12, e20, e21, e22 = _frame(x0, x1, x2)
        ps = 2.*math.pi*np.random.random()
        t1 = vloc[j, 1]*math.cos(ps) - vloc[j, 2]*math.sin(ps)
        t2 = vloc[j, 1]*math.sin(ps) + vloc[j, 2]*math.cos(ps)
        out_x[j, 0], out_x[j, 1], out_x[j, 2] = x0, x1, x2
        out_v[j, 0] = vloc[j, 0]*er0 + t1*e10 + t2*e20
        out_v[j, 1] = vloc[j, 0]*er1 + t1*e11 + t2*e21
        out_v[j, 2] = vloc[j, 0]*er2 + t1*e12 + t2*e22
    return 0


# ----------------------------------------------------------------------------------------------------------------
class RunawayError(RuntimeError):
    """Growth faster than one doubling of the tracer population per step: the run leaves the model's regime."""


class StaticBathLost(RuntimeError):
    """The bath's deficit (depletion or barrier shielding) outweighs the enclosed mass somewhere during the run:
    the net enclosed mass turns non-positive and the declared static-bath model no longer holds."""


class Model:
    """One system, one incident distribution, one sigma/m and one incident density."""

    def __init__(self, r_b, M_b, R_b, r_half, u_k, f_k, sm_cm2g, rho_inf, n_grid=1024, self_gravity=True,
                 bath_gravity=True, freeze=False, channels=('ii', 'ib', 'bb'), seed=1, eta=.08, soft=3e-3,
                 max_sub=1 << 15, n_shell=40):
        self.R_b, self.r_half = float(R_b), float(r_half)
        self.r_lo = 1e-3*r_half
        self.lr = np.linspace(math.log(self.r_lo), math.log(R_b), n_grid)
        self.r = np.exp(self.lr)
        self.lr0, self.dl, self.n = float(self.lr[0]), float(self.lr[1] - self.lr[0]), n_grid
        keep = r_b < R_b
        rb_, Mb_ = np.r_[r_b[keep], R_b], np.r_[M_b[keep], np.interp(R_b, r_b, M_b)]
        self.M_b = np.interp(self.lr, np.log(rb_), Mb_)
        self.uk, self.fk = np.asarray(u_k, float), np.asarray(f_k, float)
        self.sm_cm2g, self.sm = float(sm_cm2g), sm_cm2g*CM2_PER_G
        self.rho_inf = float(rho_inf)
        self.self_gravity, self.bath_gravity, self.freeze = self_gravity, bath_gravity, freeze
        self.channels = set(channels)
        self.eta, self.r_soft, self.max_sub, self.n_shell = eta, soft*r_half, max_sub, n_shell
        self.lmax = 10                           # at most 2^10 collision sub-steps per step
        self.tau = np.zeros(n_grid)
        self.rho_c = np.zeros(n_grid)
        self.M_bx = np.zeros(n_grid)
        self.M_c = np.zeros(n_grid)
        self.x = np.zeros((0, 3))
        self.v = np.zeros((0, 3))
        self.m = np.zeros(0)
        self.shell_m = None
        self.pool_frac = None                    # confined fraction per shell in the last pool (sets pool sizes)
        self.T_run = None                        # set by run(): zones expecting a collision within it are collisional
        seed_rng(seed)
        self.rng = np.random.default_rng(seed)
        self._set_potential(self.M_b)
        self.bath_tables()
        (self.bath_static, self.relax_iterations, self.bath_gain, self.static_fail_rho,
         self.static_fail_reason) = self.relax_bath()

    # --- potential and bath tables
    def _set_potential(self, M):
        self.M = M
        self.g = G*M/self.r**2
        # Phi(r) = -int_r^{R_b} g dr', signed: a background-subtracted bath deficit can make g negative somewhere
        self.phi = cumulative_trapezoid((self.g*self.r)[::-1], self.lr[::-1], initial=0)[::-1]

    def _phi1(self, q):
        return _phi(self.lr0, self.dl, self.n, self.phi, self.g, self.r_lo, float(q))

    def phi_at(self, rr):
        return np.array([self._phi1(q) for q in np.atleast_1d(rr)])

    def phi_grid_at(self, rr):
        return np.interp(np.log(np.maximum(rr, self.r_lo)), self.lr, self.phi)

    def bath_tables(self):
        v2 = -2*self.phi
        K = len(self.uk)
        self.dtab = np.empty((K, self.n))
        self.stab = np.empty((K, self.n))
        dep = np.exp(-self.tau)
        for k, (u, f) in enumerate(zip(self.uk, self.fk)):
            w = np.sqrt(np.maximum(u*u + v2, 0.))
            rw = self.r*w
            Lmax = np.minimum.accumulate(rw[::-1])[::-1]
            s = np.where(rw > 0, np.minimum(1., Lmax/np.where(rw > 0, rw, 1.)), 0.)
            self.stab[k] = s
            self.dtab[k] = f*(w/u)*(1 - np.sqrt(np.maximum(0., 1 - s*s)))*dep
        self.rho_bath = self.rho_inf*self.dtab.sum(axis=0)

    def bath_excess_mass(self):
        r, ex = self.r, self.rho_bath - self.rho_inf
        return cumulative_trapezoid(4*np.pi*r**3*ex, self.lr, initial=0) + 4/3*np.pi*r[0]**3*ex[0]

    def total_mass(self):
        M = self.M_b.copy()
        if self.self_gravity and len(self.x):
            M = M + self.M_c
        if self.bath_gravity:
            M = M + self.M_bx
        return M

    def _deplete_only(self, tol=1e-10, itmax=200):
        """Depletion consistent with the bath at fixed potential (no bath gravity)."""
        for _ in range(itmax):
            old = self.tau.copy()
            self.bath_tables()
            self.update_tau()
            if np.max(np.abs(self.tau - old)) < tol*max(1., float(np.max(self.tau))):
                break
        self.bath_tables()

    def _bath_map(self, x):
        """F(x): the bath's excess enclosed mass in the potential of baryons, tracers and bath excess x, with the
        depletion made consistent (one update of the optical depth)."""
        self.M_bx = x
        self._set_potential(self.total_mass())
        self.bath_tables()
        self.update_tau()
        self.bath_tables()
        return self.bath_excess_mass()

    def _relax_at(self, tol, itmax, depth=6, beta=.5):
        """Fixed point x = F(x) by Anderson mixing, which also converges where the bath's response is a strong
        negative feedback (the centrifugal barrier shields the interior as the enclosed mass grows)."""
        scale = float(self.M_b[-1])
        x = self.M_bx.copy()
        X, Rs = [], []
        for it in range(1, itmax + 1):
            Fx = self._bath_map(x)
            if not np.all(np.isfinite(Fx)) or np.max(np.abs(Fx)) > 1e3*scale:
                return False, it
            res = Fx - x
            if float(np.max(np.abs(res)))/(scale + float(np.max(np.abs(Fx)))) < tol:
                self._bath_map(Fx)
                self.M_bx = Fx
                return True, it
            X.append(x.copy())
            Rs.append(res.copy())
            if len(X) > depth + 1:
                X.pop(0)
                Rs.pop(0)
            if len(X) == 1:
                x = x + beta*res
            else:
                dR = np.array([Rs[i + 1] - Rs[i] for i in range(len(Rs) - 1)]).T
                dX = np.array([X[i + 1] - X[i] for i in range(len(X) - 1)]).T
                gam = np.linalg.lstsq(dR, res, rcond=None)[0]
                x = x + beta*res - (dX + beta*dR) @ gam
        return False, itmax

    def _gain(self, eps=1e-4, iters=20):
        """Eigenvalues of the linearised map M_bx -> F(M_bx) at the converged state, with the depletion held fixed:
        the one of largest magnitude, and the largest positive one (by a shifted power iteration when the first is
        negative). Negative feedback comes from the centrifugal barrier; a static bath runs away gravitationally only
        if the largest positive eigenvalue reaches one."""
        saved = (self.M.copy(), self.phi.copy(), self.g.copy(), self.M_bx.copy(), self.dtab.copy(), self.stab.copy(),
                 self.rho_bath.copy())
        base_bx = self.M_bx.copy()
        amp = eps*max(float(np.max(np.abs(base_bx))), 1e-3*float(self.M_b[-1]))

        def apply(vec):
            self.M_bx = base_bx + amp*vec
            self._set_potential(self.total_mass())
            self.bath_tables()
            return (self.bath_excess_mass() - F0)/amp
        self._set_potential(self.total_mass())
        self.bath_tables()
        F0 = self.bath_excess_mass()

        def power(shift):
            vec = np.abs(base_bx) + 1e-3*float(self.M_b[-1])*self.r/self.r[-1]
            vec /= np.max(np.abs(vec))
            lam = 0.
            for _ in range(iters):
                w = apply(vec) - shift*vec
                nrm = float(np.max(np.abs(w)))
                if nrm == 0.:
                    return shift
                lam = float(w @ vec/(vec @ vec))
                vec = w/nrm
            return lam + shift
        lam1 = power(0.)
        lam_pos = lam1 if lam1 >= 0 else power(lam1)
        self.M, self.phi, self.g, self.M_bx, self.dtab, self.stab, self.rho_bath = saved
        return lam1, lam_pos

    def relax_bath(self, tol=1e-10, itmax=300, n_cont=12):
        """The bath with its own excess gravity and depletion, made consistent at fixed tracers by continuation from
        the weightless limit (1e-3 of the density) in n_cont geometric steps. At each step the fixed point must
        converge, the enclosed mass must stay positive, and the largest positive eigenvalue of the focusing response
        must stay below one; otherwise the static transparent bath has no stable solution at this density.
        Returns (static, iterations, (dominant, largest positive) eigenvalues at the last density reached,
        density at which it failed, reason)."""
        if not self.bath_gravity:
            self._deplete_only()
            return True, 0, None, None, None
        target, its, gain = self.rho_inf, 0, None
        for rho in np.geomspace(1e-3*target, target, n_cont):
            self.rho_inf = float(rho)
            ok, k = self._relax_at(tol, itmax)
            its += k
            reason = None if ok else 'fixed point not found'
            if ok and not np.all(self.M > 0):
                ok, reason = False, 'enclosed mass not positive'
            if ok:
                gain = self._gain()
                if gain[1] >= 1.:
                    ok, reason = False, 'gravitational runaway (positive gain >= 1)'
            if not ok:
                self.rho_inf = target
                return False, its, gain, float(rho), reason
        self.rho_inf = target
        return True, its, gain, None, None

    # --- tracers
    def _sorted(self):
        rr = np.linalg.norm(self.x, axis=1)
        order = np.argsort(rr)
        return rr, order

    def enclosed_tracer_mass(self):
        if not len(self.x):
            return np.zeros(self.n)
        rr, order = self._sorted()
        cm = np.cumsum(self.m[order])
        k = np.searchsorted(rr[order], self.r, side='right')
        return np.where(k > 0, cm[np.maximum(k - 1, 0)], 0.)

    def tracer_density(self, k=16):
        """Local confined density at each tracer (mass within k radial neighbours on each side) and on the grid."""
        N = len(self.x)
        rr, order = self._sorted()
        if N < 2:
            return np.zeros(N), order, np.zeros(self.n)
        rs = rr[order]
        cm = np.r_[0., np.cumsum(self.m[order])]
        idx = np.arange(N)
        lo, hi = np.clip(idx - k, 0, N - 1), np.clip(idx + k, 0, N - 1)
        # the radial window is at least 2% of r: resampled clones and repeated pool births share a radius exactly
        vol = np.maximum(4/3*np.pi*(rs[hi]**3 - rs[lo]**3), 4*np.pi*rs**3*.02)
        mass = cm[hi + 1] - cm[lo] - .5*(self.m[order][hi] + self.m[order][lo])
        rho_sorted = np.maximum(mass, 0.)/vol
        rho = np.empty(N)
        rho[order] = rho_sorted
        grid = np.interp(self.lr, np.log(rs), rho_sorted, left=rho_sorted[0], right=0.)
        return rho, order, grid

    def update_tau(self):
        """First-order depletion: optical depth along the radial inflow path, counting the bath (mean relative speed
        4w/3 for equal speeds) and the confined companions (relative speed about w)."""
        kappa = self.sm*(4/3*self.rho_bath + self.rho_c)
        self.tau = np.abs(cumulative_trapezoid((kappa*self.r)[::-1], self.lr[::-1], initial=0)[::-1])

    def _bb_phase(self, Delta, alive):
        """Bound-bound collisions. A tracer that expects at least one collision over the run (P*T/Delta > 1, with
        P = sm*rho_c*2|v|*Delta) lies in a collisional zone; zones are contiguous in radial order and smoothed over 33
        neighbours. In a collisional zone the tracers are resampled to the zone's mean mass (unbiased; booked), so
        equal masses collide and nothing is split, and a zone of level l (P up to 0.2*2^l) takes 2^l passes of
        Delta/2^l at fixed positions, with partners redrawn among near radial neighbours in every pass. Elsewhere
        collisions are rare, and the heavier tracer of a colliding pair is split so that equal masses collide."""
        lr0, dl, n = self.lr0, self.dl, self.n
        rho_ci, order_all, grid = self.tracer_density()
        sp = np.linalg.norm(self.v, axis=1)
        P = self.sm*rho_ci*2*sp*Delta
        lev = np.clip(np.ceil(np.log2(np.maximum(P/.2, 1.))), 0, self.lmax).astype(np.int64)
        coll = (lev > 0) if self.T_run is None else (lev > 0) | (P*self.T_run/Delta > 1.)
        lev_s = maximum_filter1d(lev[order_all], size=33, mode='nearest')
        coll_s = maximum_filter1d(coll[order_all].astype(np.int64), size=33, mode='nearest')
        zone = 2*lev_s + coll_s
        bounds = np.r_[0, np.nonzero(np.diff(zone))[0] + 1, len(zone)]
        led = np.zeros(9)
        out = dict(P_max=0., P_over_1=0., level_max=int(lev_s.max()), capped_mass_fraction=0., M_resample=0.,
                   E_resample=0., resampled_zones=0, collisional_mass_fraction=0.)
        newx, newv, newm, exE, exM = [], [], [], [], []
        capped = collisional = 0.

        def absorb(lb, ps, bE, bM):
            for q in (0, 1, 2, 3, 5, 6):
                led[q] += lb[q]
            led[4], led[7] = max(led[4], lb[4]), max(led[7], lb[7])
            out['P_max'] = max(out['P_max'], ps[0])
            out['P_over_1'] += ps[1]
            ne = int(lb[8])
            if ne:
                exE.append(bE[:ne].copy())
                exM.append(bM[:ne].copy())

        def buffers(k2):
            return np.zeros((k2, 3)), np.zeros((k2, 3)), np.zeros(k2), np.zeros(k2), np.zeros(k2), np.zeros(9), np.zeros(2)
        for a, b in zip(bounds[:-1], bounds[1:]):
            idx = order_all[a:b]
            idx = idx[alive[idx]]
            if len(idx) < 2:
                continue
            level = int(lev_s[a])
            if level >= self.lmax:
                capped += float(self.m[idx].sum())
            npass = 1 << level
            if not coll_s[a]:
                nx_, nv_, nm_, bE, bM, lb, ps = buffers(2*len(idx))
                k = bound_bound(self.x, self.v, self.m, alive, idx, rho_ci, Delta, self.sm, lr0, dl, n, self.phi,
                                self.g, self.r_lo, nx_, nv_, nm_, bE, bM, lb, ps)
                absorb(lb, ps, bE, bM)
                if k:
                    newx.append(nx_[:k].copy()); newv.append(nv_[:k].copy()); newm.append(nm_[:k].copy())
                continue
            n0 = len(idx)
            xz, vz, mz = self.x[idx].copy(), self.v[idx].copy(), self.m[idx].copy()
            collisional += float(mz.sum())
            out['resampled_zones'] += 1
            # equal masses in a collisional zone: unbiased resampling to the zone's mean tracer mass (booked). Splitting
            # there made a new tracer at every collision of unequal masses; in dense seed cores the count doubled
            # every few steps and the repeated thinning that followed swamped the confined mass.
            mzone = float(mz.sum())/n0
            cnt = np.floor(mz/mzone).astype(np.int64)
            cnt += (self.rng.random(n0) < (mz/mzone - cnt)).astype(np.int64)
            Ez = energies(xz, vz, lr0, dl, n, self.phi, self.g, self.r_lo)
            out['M_resample'] += mzone*float(cnt.sum()) - float(mz.sum())
            out['E_resample'] += mzone*float(cnt @ Ez) - float(mz @ Ez)
            rep = np.repeat(np.arange(n0), cnt)
            xz, vz, mz = xz[rep], vz[rep], np.full(len(rep), mzone)
            alive[idx] = False                  # replaced by the resampled set, returned as new tracers
            alz = np.ones(len(rep), np.bool_)
            rz, rhoz = np.linalg.norm(xz, axis=1), rho_ci[idx][rep]
            for _ in range(npass):
                live = np.nonzero(alz)[0]
                if len(live) < 2:
                    break
                o = live[np.argsort(rz[live], kind='stable')]
                # fresh partners in every pass: shuffle within windows of 8 radial neighbours (random offset). With
                # positions frozen, the radial order alone would pair the same two tracers in all 2^l passes.
                win = (np.arange(len(o)) + int(self.rng.integers(8)))//8
                o = o[np.lexsort((self.rng.random(len(o)), win))]
                nx_, nv_, nm_, bE, bM, lb, ps = buffers(2*len(xz))
                k = bound_bound(xz, vz, mz, alz, o, rhoz, Delta/npass, self.sm, lr0, dl, n, self.phi, self.g,
                                self.r_lo, nx_, nv_, nm_, bE, bM, lb, ps)
                absorb(lb, ps, bE, bM)
                if k:
                    rk = np.linalg.norm(nx_[:k], axis=1)
                    xz, vz, mz = np.concatenate([xz, nx_[:k]]), np.concatenate([vz, nv_[:k]]), np.concatenate([mz, nm_[:k]])
                    alz = np.concatenate([alz, np.ones(k, np.bool_)])
                    rz, rhoz = np.concatenate([rz, rk]), np.concatenate([rhoz, np.interp(np.log(rk), self.lr, grid)])
            keep = np.nonzero(alz)[0]
            if len(keep):
                newx.append(xz[keep]); newv.append(vz[keep]); newm.append(mz[keep])
        out['capped_mass_fraction'] = capped/max(float(self.m.sum()), 1e-300)
        out['collisional_mass_fraction'] = collisional/max(float(self.m.sum()), 1e-300)
        out['led'] = led
        cat = lambda z, shape: np.concatenate(z) if z else np.zeros(shape)
        out['new'] = (cat(newx, (0, 3)), cat(newv, (0, 3)), cat(newm, 0))
        out['exE'], out['exM'] = cat(exE, 0), cat(exM, 0)
        return out

    def add_tracers(self, x, v, m):
        self.x = np.concatenate([self.x, x])
        self.v = np.concatenate([self.v, v])
        self.m = np.concatenate([self.m, m])

    def energy_sum(self):
        if not len(self.x):
            return 0.
        return float(self.m @ energies(self.x, self.v, self.lr0, self.dl, self.n, self.phi, self.g, self.r_lo))

    def sample_seed(self, M_seed, sigma_grid, n_seed, margin=1e-3):
        """CF-1's counted seed: baryon-shaped, isotropic Maxwellian with the Jeans dispersion, truncated to bound
        orbits as in CF-1 (E < 0). The cut is E < -margin*|Phi(r_half)|, so that no seed companion sits within the
        orbit integrator's energy tolerance of unbinding (one of 10,000 left J1630 that way in a 10 Gyr test)."""
        N = int(n_seed)
        Mb = self.M_b
        u = self.rng.uniform(Mb[0], Mb[-1], N)
        r = np.exp(np.interp(u, Mb, self.lr))
        s = np.interp(np.log(r), self.lr, sigma_grid)
        v = self.rng.normal(size=(N, 3))*s[:, None]
        Ecut = -margin*abs(self._phi1(self.r_half))
        phi_r = self.phi_at(r)
        for it in range(1000):
            E = .5*np.einsum('ij,ij->i', v, v) + phi_r
            ok = E < Ecut
            if ok.all():
                break
            bad = ~ok
            # a position with Phi above the cut cannot hold such an orbit, and one just below it rarely does:
            # redraw those positions (the seed's outermost tail, a negligible mass) along with the velocities
            redraw = bad & ((phi_r >= Ecut) | (it % 25 == 24))
            if redraw.any():
                rr_ = np.exp(np.interp(self.rng.uniform(Mb[0], Mb[-1], int(redraw.sum())), Mb, self.lr))
                r[redraw], s[redraw] = rr_, np.interp(np.log(rr_), self.lr, sigma_grid)
                phi_r[redraw] = self.phi_at(rr_)
            v[bad] = self.rng.normal(size=(int(bad.sum()), 3))*s[bad, None]
        else:
            raise RuntimeError('seed sampling did not converge')
        xs, vs = np.empty((N, 3)), np.empty((N, 3))
        place(r, v, xs, vs)
        self.add_tracers(xs, vs, np.full(N, M_seed/N))

    # --- incoming-incoming production pools, one per logarithmic shell
    def make_pools(self, n_per=1000, n_min=100, n_cap=20000):
        """One pool of seedless events per logarithmic shell. A shell's number of draws is set from the previous pool's
        fraction of draws that left a companion confined (from a pilot of n_per draws the first time), aiming at n_min
        confined events within n_per to n_cap draws. Most of the Milky Way's production comes from outer shells where
        only 3-50 of 1,000 draws confine a companion, and births drawn from so few events are clones. The size comes
        from an earlier, independent sample, so each pool remains an unbiased estimate."""
        w = 4*np.pi*self.r**3*(self.rho_bath/self.rho_inf)**2
        cum = cumulative_trapezoid(w, self.lr, initial=0)
        edges = np.linspace(self.lr[0], self.lr[-1], self.n_shell + 1)
        ce = np.interp(edges, self.lr, cum)
        live = [bool(ce[j + 1] - ce[j] > 1e-12*max(float(cum[-1]), 1e-300)) for j in range(self.n_shell)]

        def draw(j, k):
            q = self.rng.uniform(ce[j], ce[j + 1], k)
            r_s = np.clip(np.exp(np.interp(q, cum, self.lr)), math.exp(edges[j]), math.exp(edges[j + 1]))
            out = np.empty((k, 15))
            seedless_events(r_s, self.lr0, self.dl, self.n, self.dtab, self.stab, self.uk, self.phi, self.g, self.r_lo, out)
            return out
        if self.pool_frac is None:
            self.pool_frac = np.array([float(np.mean(draw(j, n_per)[:, 2] >= 1)) if live[j] else 1.
                                       for j in range(self.n_shell)])
        sizes = np.clip(np.ceil(n_min/np.maximum(self.pool_frac, n_min/n_cap)), n_per, n_cap).astype(np.int64)
        frac = np.ones(self.n_shell)
        self.pools = []
        for j in range(self.n_shell):
            W = ce[j + 1] - ce[j]
            if not live[j]:
                # no bath left in this shell (fully depleted): nothing is produced there
                self.pools.append(dict(W=0., r=(math.exp(edges[j]), math.exp(edges[j + 1])), g_conf=0., g_n=0.,
                                       g_bound=0., ev=np.zeros((0, 15)), p=None, draws=0))
                continue
            out = draw(j, int(sizes[j]))
            g, nc = out[:, 1], out[:, 2]
            frac[j] = float(np.mean(nc >= 1))
            ev = out[nc >= 1]
            self.pools.append(dict(W=W, r=(math.exp(edges[j]), math.exp(edges[j + 1])), g_conf=float(np.mean(g*(nc >= 1))),
                                   g_n=float(np.mean(g*nc)), g_bound=float(np.mean(g*out[:, 14])),
                                   ev=ev, p=(ev[:, 1]/ev[:, 1].sum() if len(ev) else None), draws=int(sizes[j])))
        self.pool_frac = frac
        bad = [j for j, p in enumerate(self.pools) if not (np.isfinite(p['W']) and np.isfinite(p['g_n']) and np.isfinite(p['g_conf']))]
        if bad:
            raise RuntimeError(f'non-finite production pools in shells {bad[:5]}: finite phi={np.isfinite(self.phi).all()} '
                               f'g={np.isfinite(self.g).all()} rho_bath={np.isfinite(self.rho_bath).all()} '
                               f'tau={np.isfinite(self.tau).all()} M={np.isfinite(self.M).all()} tracers={len(self.x)}')
        self.pool_phi = self._phi1(self.r_half)
        return self.pools

    def shell_rates(self):
        """Seedless production of confined mass per unit time in each shell."""
        return np.array([.5*self.sm*self.rho_inf**2*p['W']*p['g_n'] for p in self.pools])

    def production_rate(self):
        return float(self.shell_rates().sum())

    def set_shell_masses(self, T, n_target, weight_range=10.):
        """Tracer mass per birth shell proportional to the square root of its production rate, normalised so that
        the initial rates would give n_target births in T, and kept within a factor weight_range of the heaviest:
        wider ranges make bound-bound splitting multiply light tracers in dense cores. A shell with no confined event
        in its first pool takes the mass of the nearest shell that has one, so it can produce later. (An infinite mass
        had switched such a shell off for the whole run: the Milky Way's outermost shell holds a quarter of the
        production but only 0-3 confined events per 1,000 draws, so the outcome depended on the first pool's luck.)"""
        rates = self.shell_rates()
        root = np.sqrt(np.maximum(rates, 0.))
        if not root.max() > 0:
            self.m_floor = self.m_top = 1.
            self.shell_m = np.ones(len(root))
            return
        c = T*root.sum()/n_target
        top = c*float(root.max())
        self.m_floor, self.m_top = top/weight_range, top
        on = np.nonzero(root > 0)[0]
        near = on[np.abs(np.arange(len(root))[:, None] - on[None, :]).argmin(axis=1)]
        self.shell_m = np.clip(c*root[near], self.m_floor, top)

    def expected_births(self, Delta):
        return sum(.5*self.sm*self.rho_inf**2*p['W']*p['g_conf']*Delta/self.shell_m[j]
                   for j, p in enumerate(self.pools) if p['p'] is not None and np.isfinite(self.shell_m[j]))

    def _thin(self, L):
        """Halve the tracer count by merging radial neighbours in pairs. A merged tracer carries both masses and
        takes the position and velocity of one of the two, chosen with probability proportional to its mass: the
        mass is exact and the energy is conserved in expectation (booked). Random removal of half the tracers was
        unbiased too, but with unequal masses one removal could delete a large share of the total; repeated in the
        dense seed cores it swung the confined mass by factors of tens. Birth masses double with the merge."""
        Et, Mt = self.energy_sum(), float(self.m.sum())
        rr, order = self._sorted()
        n2 = (len(order)//2)*2
        a, b = order[0:n2:2], order[1:n2:2]
        ma, mb = self.m[a], self.m[b]
        keep = np.where(self.rng.random(len(a)) < mb/(ma + mb), b, a)
        rest = order[n2:]
        idx = np.r_[keep, rest]
        self.x, self.v = self.x[idx], self.v[idx]
        self.m = np.r_[ma + mb, self.m[rest]]
        self.shell_m = 2*self.shell_m
        self.m_floor, self.m_top = 2*self.m_floor, 2*self.m_top
        L['E_thinning'] += self.energy_sum() - Et
        L['M_thinning'] += float(self.m.sum()) - Mt
        L['thinnings'] += 1

    def _heavier_births(self, L):
        """Birth-rate control: double the mass of future births (their number per step halves). Existing tracers
        are left alone; thinning them here, repeatedly when production rose sharply, destroyed the statistics."""
        self.shell_m = 2*self.shell_m
        self.m_floor, self.m_top = 2*self.m_floor, 2*self.m_top
        L['birth_mass_doublings'] += 1

    def births(self, Delta):
        """Incoming-incoming births in this step. Returns the new tracers, the incoming energy consumed, and the
        energies and masses of the unconfined partners (exported)."""
        xs, vs, ms, E_in, pE, pm = [], [], [], 0., [], []
        for j, p in enumerate(self.pools):
            if p['p'] is None or not np.isfinite(self.shell_m[j]):
                continue
            mj = self.shell_m[j]
            lam = .5*self.sm*self.rho_inf**2*p['W']*p['g_conf']*Delta/mj
            k = self.rng.poisson(lam)
            if not k:
                continue
            ev = p['ev'][self.rng.choice(len(p['ev']), size=k, p=p['p'])]
            E_in += mj*float(ev[:, 3].sum())
            for c, a, e in ((4, 5, 8), (9, 10, 13)):
                sel = ev[:, c] > .5
                if (~sel).any():
                    pE.append(ev[~sel, e].copy())
                    pm.append(np.full(int((~sel).sum()), mj))
                if sel.any():
                    x = np.empty((int(sel.sum()), 3))
                    v = np.empty_like(x)
                    place(ev[sel, 0].copy(), ev[sel, a:a + 3].copy(), x, v)
                    xs.append(x)
                    vs.append(v)
                    ms.append(np.full(len(x), mj))
        cat = lambda z, shape: np.concatenate(z) if z else np.zeros(shape)
        return (cat(xs, (0, 3)), cat(vs, (0, 3)), cat(ms, 0), E_in, cat(pE, 0), cat(pm, 0))

    # --- evolution
    def run(self, T_gyr, Delta_max_gyr, snaps_gyr, n_target=16000, n_max=40000, budget_s=1800., log=None, regen=10):
        """Evolve for T. Returns the history, ledgers and exits; raises RuntimeError on a non-finite state and
        TimeoutError when the wall-clock budget is exhausted."""
        t0 = time.time()
        T, Dmax = T_gyr*PER_GYR, Delta_max_gyr*PER_GYR
        self.T_run = T
        snaps = [s*PER_GYR for s in snaps_gyr]
        L = dict(collisions_ib=0., captures=0., ejections=0., collisions_bb=0., evaporations=0., births=0.,
                 orbit_escapes=0., captured_mass=0., ejected_mass=0., evaporated_mass=0., escaped_mass=0., born_mass=0.,
                 E_births=0., E_bath_in=0., E_exported=0., E_potential_work=0., E_integration=0., E_thinning=0.,
                 E_pool_lag=0., M_thinning=0., M_resample=0., E_resample=0., birth_mass_doublings=0,
                 E_orbit_escapes=0., dE_ib=0., dE_bb=0., closure_max=0., max_dK=0., max_dP=0., P_max=0., P_over_1=0.,
                 substeps=0., capped=0., steps=0, regenerations=0, thinnings=0, splits=0., longest_chain=0.,
                 ib_substeps_max=0, bb_level_max=0, capped_mass_fraction_max=0., collisional_mass_fraction_max=0.,
                 bath_max_change=0., smallest_Delta_Gyr=Delta_max_gyr)
        hist, exits_E, exits_m, exits_c = [], [], [], []
        lr0, dl, n = self.lr0, self.dl, self.n
        self.M_c = self.enclosed_tracer_mass()
        self.make_pools()
        if self.shell_m is None:
            self.set_shell_masses(T, n_target)
        L['regenerations'] += 1
        t, si, step, last_regen = 0., 0, 0, 0
        if snaps and snaps[0] <= 0:
            hist.append(self.snapshot(0., L))
            si = 1
        while t < T*(1 - 1e-12):
            N = len(self.x)
            Delta = min(Dmax, T - t)          # collisions faster than 0.2 per step are sub-cycled below
            L['smallest_Delta_Gyr'] = min(L['smallest_Delta_Gyr'], Delta/PER_GYR)
            # 1. orbits in the frozen potential of this step
            if N:
                alive = np.ones(N, np.bool_)
                e_out = np.zeros(N)
                E0 = energies(self.x, self.v, lr0, dl, n, self.phi, self.g, self.r_lo)
                sub, cap = advance(self.x, self.v, alive, Delta, lr0, dl, n, self.g, self.phi, self.r_lo, self.R_b,
                                   self.eta, self.r_soft, self.max_sub, e_out)
                L['substeps'] += sub
                L['capped'] += cap
                E1 = energies(self.x, self.v, lr0, dl, n, self.phi, self.g, self.r_lo)
                gone = ~alive
                L['E_integration'] += float(self.m[alive] @ E1[alive] + self.m[gone] @ e_out[gone] - self.m @ E0)
                if gone.any():
                    L['orbit_escapes'] += float(gone.sum())
                    L['escaped_mass'] += float(self.m[gone].sum())
                    L['E_orbit_escapes'] += float(self.m[gone] @ e_out[gone])
                    exits_E.append(e_out[gone]); exits_m.append(self.m[gone]); exits_c.append(np.full(int(gone.sum()), 0))
                    self.x, self.v, self.m = self.x[alive], self.v[alive], self.m[alive]
                if not (np.isfinite(self.x).all() and np.isfinite(self.v).all()):
                    raise RuntimeError('non-finite tracer state')
            # 2. collisions at fixed positions. Birth control first: at most n_max/20 births per step, otherwise thin
            # and double the birth masses (outside the closure check below, since thinning is booked on its own)
            if 'ii' in self.channels:
                while self.expected_births(Delta) > n_max/20:
                    self._heavier_births(L)
            N = len(self.x)
            Ebefore = self.energy_sum()
            booked = 0.
            alive = np.ones(N, np.bool_)
            e_rem = np.zeros(N)
            ns = nn = 0
            if N and 'ib' in self.channels:
                rr = np.linalg.norm(self.x, axis=1)
                sp = np.linalg.norm(self.v, axis=1)
                rb = np.interp(np.log(np.maximum(rr, self.r_lo)), self.lr, self.rho_bath)
                wmax = np.sqrt(np.maximum(self.uk.max()**2 - 2*self.phi_grid_at(rr), 0.))
                cap_n = 1 << self.lmax
                nsub = np.clip(np.ceil(self.sm*rb*(wmax + sp)*Delta/.2), 1, cap_n).astype(np.int64)
                L['ib_substeps_max'] = max(L['ib_substeps_max'], int(nsub.max()))
                L['capped_mass_fraction_max'] = max(L['capped_mass_fraction_max'],
                                                    float(self.m[nsub >= cap_n].sum())/float(self.m.sum()))
                K_ = int(min(nsub.sum(), 2*N + 1000))     # more captures than that in one step is a runaway
                sx, sv, smass = np.zeros((K_, 3)), np.zeros((K_, 3)), np.zeros(K_)
                lb, ps = np.zeros(12), np.zeros(2)
                ns = bath_bound(self.x, self.v, self.m, alive, nsub, Delta, self.sm, self.rho_inf, lr0, dl, n, self.dtab,
                                self.stab, self.uk, self.phi, self.g, self.r_lo, False, sx, sv, smass, e_rem, lb, ps)
                if lb[11]:
                    raise RunawayError(f'captures exceed twice the tracer count in one step at t={t/PER_GYR:.3f} Gyr')
                L['collisions_ib'] += lb[0]; L['captures'] += lb[1]; L['ejections'] += lb[2]
                L['dE_ib'] += lb[3]; L['E_bath_in'] += lb[4]; L['E_exported'] += lb[5]
                L['max_dK'] = max(L['max_dK'], lb[6]); L['max_dP'] = max(L['max_dP'], lb[7])
                L['captured_mass'] += lb[8]; L['ejected_mass'] += lb[9]
                L['P_max'] = max(L['P_max'], ps[0]); L['P_over_1'] += ps[1]
                booked += lb[3]
                gone = ~alive
                if gone.any():
                    exits_E.append(e_rem[gone]); exits_m.append(self.m[gone]); exits_c.append(np.full(int(gone.sum()), 1))
            bbnew = None
            if N > 1 and 'bb' in self.channels:
                bb = self._bb_phase(Delta, alive)
                lb = bb['led']
                L['collisions_bb'] += lb[0]; L['evaporations'] += lb[1]; L['dE_bb'] += lb[2]
                L['E_exported'] += lb[3]; L['max_dK'] = max(L['max_dK'], lb[4]); L['evaporated_mass'] += lb[5]
                L['splits'] += lb[6]; L['longest_chain'] = max(L['longest_chain'], lb[7])
                L['P_max'] = max(L['P_max'], bb['P_max']); L['P_over_1'] += bb['P_over_1']
                L['bb_level_max'] = max(L['bb_level_max'], bb['level_max'])
                L['capped_mass_fraction_max'] = max(L['capped_mass_fraction_max'], bb['capped_mass_fraction'])
                L['collisional_mass_fraction_max'] = max(L['collisional_mass_fraction_max'], bb['collisional_mass_fraction'])
                L['M_resample'] += bb['M_resample']
                L['E_resample'] += bb['E_resample']
                booked += lb[2] + bb['E_resample']
                if len(bb['exE']):
                    exits_E.append(bb['exE']); exits_m.append(bb['exM']); exits_c.append(np.full(len(bb['exE']), 2))
                bbnew = bb['new']
            if N:
                self.x, self.v, self.m = self.x[alive], self.v[alive], self.m[alive]
            if ns:
                self.add_tracers(sx[:ns].copy(), sv[:ns].copy(), smass[:ns].copy())
            if bbnew is not None and len(bbnew[2]):
                self.add_tracers(*bbnew)
            if 'ii' in self.channels:
                bx, bv, bm, E_in, pE, pm = self.births(Delta)
                if len(bx) or len(pE):
                    Eb = float(bm @ energies(bx, bv, lr0, dl, n, self.phi, self.g, self.r_lo)) if len(bx) else 0.
                    L['births'] += len(bx); L['born_mass'] += float(bm.sum())
                    L['E_births'] += Eb; L['E_bath_in'] += E_in; L['E_exported'] += E_in - Eb
                    L['E_pool_lag'] += (E_in - Eb) - float(pm @ pE)
                    booked += Eb
                    if len(bx):
                        self.add_tracers(bx, bv, bm)
                    if len(pE):
                        exits_E.append(pE); exits_m.append(pm); exits_c.append(np.full(len(pE), 3))
            Eafter = self.energy_sum()
            scale = abs(Ebefore) + abs(Eafter) + abs(booked) + 1e-300
            L['closure_max'] = max(L['closure_max'], abs(Eafter - Ebefore - booked)/scale)
            # 3. thinning: merge radial neighbours in pairs (and double future birth masses)
            while len(self.x) > n_max:
                self._thin(L)
            # 4. potential, bath and depletion
            Eb4 = self.energy_sum()
            self.M_c = self.enclosed_tracer_mass()
            self.rho_c = self.tracer_density()[2] if len(self.x) > 1 else np.zeros(self.n)
            self.update_tau()
            if not self.freeze:
                old_bx = self.M_bx.copy()
                self.bath_tables()
                if self.bath_gravity:
                    self.M_bx = self.bath_excess_mass()
                    L['bath_max_change'] = max(L['bath_max_change'], float(np.max(np.abs(self.M_bx - old_bx)))/float(self.M_b[-1]))
                self._set_potential(self.total_mass())
                if self.bath_gravity and not np.all(self.M > 0):
                    i_bad = int(np.argmax(self.M <= 0))
                    raise StaticBathLost(f'net enclosed mass non-positive at r={self.r[i_bad]:.4g} kpc, '
                                         f't={(t + Delta)/PER_GYR:.3f} Gyr')
            self.bath_tables()
            if not (np.isfinite(self.phi).all() and np.isfinite(self.rho_bath).all() and np.isfinite(self.tau).all()):
                raise RuntimeError(f'non-finite tables at t={(t + Delta)/PER_GYR:.3f} Gyr: finite M={np.isfinite(self.M).all()} '
                                   f'M_c={np.isfinite(self.M_c).all()} M_bx={np.isfinite(self.M_bx).all()} '
                                   f'rho_c={np.isfinite(self.rho_c).all()} tau={np.isfinite(self.tau).all()} '
                                   f'masses={np.isfinite(self.m).all()} x={np.isfinite(self.x).all()}')
            L['E_potential_work'] += self.energy_sum() - Eb4
            t += Delta
            step += 1
            L['steps'] = step
            if (step - last_regen >= regen) or (step - last_regen >= 2 and abs(self._phi1(self.r_half)/self.pool_phi - 1) > 1e-2):
                self.make_pools()
                last_regen = step
                L['regenerations'] += 1
            while si < len(snaps) and t >= snaps[si]*(1 - 1e-9):
                hist.append(self.snapshot(t, L))
                if log:
                    h = hist[-1]
                    log(f"    t={t/PER_GYR:6.2f} Gyr  M={h['M_confined']:.3e}  M(<r_half)={h['M_within_r_half']:.3e}  "
                        f"N={h['tracers']}  tau(r_half)={h['tau_r_half']:.3g}  P_max={L['P_max']:.2g}  "
                        f"{time.time() - t0:.0f}s")
                si += 1
            if time.time() - t0 > budget_s:
                raise TimeoutError(f'wall-clock budget exhausted at t={t/PER_GYR:.3f} Gyr')
        ex = dict(E=np.concatenate(exits_E) if exits_E else np.zeros(0), m=np.concatenate(exits_m) if exits_m else np.zeros(0),
                  channel=np.concatenate(exits_c) if exits_c else np.zeros(0, int))
        return dict(history=hist, ledger=L, exits=ex, runtime_s=time.time() - t0)

    def snapshot(self, t, L):
        N = len(self.x)
        rr = np.linalg.norm(self.x, axis=1) if N else np.zeros(0)
        return dict(t_Gyr=t/PER_GYR, tracers=N, M_confined=float(self.m.sum()),
                    M_within_r_half=float(self.m[rr < self.r_half].sum()),
                    tau_r_half=float(np.interp(math.log(self.r_half), self.lr, self.tau)), tau_centre=float(self.tau[0]),
                    bath_excess_within_r_half=float(np.interp(math.log(self.r_half), self.lr, self.M_bx)),
                    production_rate_Msun_per_Gyr=self.production_rate()*PER_GYR,
                    events=dict(births=L['births'], captures=L['captures'], ejections=L['ejections'],
                                evaporations=L['evaporations'], orbit_escapes=L['orbit_escapes']))

"""RUT-1 stage 5: linear modes of a rigidly rotating ring of writers (protocol-rut5.md).

Stage 4's code is frozen -- `rut4_checks.py` pins the hashes of formation.py, longrun.py, rut3.py, rut4.py
and rut1.py to what ran -- so this stage adds new files and imports the old ones unchanged.

The force on body i is the exact history integral

    a_i(t) = sum_j q_j int_0^inf h(u) grad K(X_i(t) - X_j(t-u)) du

with K the kernel and h the response kernel, the inverse Laplace transform of H(s):

    instantaneous   H = tau_keep                                   h(u) = tau_keep delta(u)
    one-stage       H = tau_keep/(1 + s tau_keep)                  h(u) = exp(-u/tau_keep)
    two-stage       H = tau_keep/[(1+s tau_keep)(1+s tau_form)]    h(u) = tk/(tk-tf) [exp(-u/tk) - exp(-u/tf)]

All three have int h = tau_keep, so they share the static gain and differ only in the lag: that is the
matched control the stage 4 attribution was missing. A fifth rung uses pairwise Newtonian attraction with
no lag, to ask whether the mode is just the classical instability of a discrete ring of attracting bodies.

Linearising about the rigidly rotating ring, in each body's own radial/tangential frame, with the discrete
Bloch decomposition u_j = u_hat exp(2 pi i m j/N) exp(s t):

    det[ s^2 I - 2 Omega s J - Omega^2 I - D - A + B_m(s) ] = 0,   J = [[0,1],[-1,0]]
    D = diag(2GM/R^3, -GM/R^3)
    A = q sum_j int h(u) HessK(d_j(u)) du                                   (the base field's Hessian)
    B_m(s) = q sum_j int h(u) e^{-su} e^{2 pi i m j/N} HessK(d_j(u)) Rot(phi_j - Omega u) du
    d_j(u) = P_0 - Rot(-Omega u) P_j

CORRECTION to the protocol's displayed formula: it wrote the last factor as Rot(-Omega u), omitting the
per-body basis rotation Lambda_j = Rot(phi_j). The displacement of body j is expressed in body j's own
radial/tangential frame -- that is what makes the problem invariant under j -> j+1 and so Bloch-separable --
so converting it to body 0's frame carries Rot(phi_j - Omega u). The implementation uses the complete
expression; the check that catches this is G1, which compares A and B against finite differences of the
full history force.

**The history sums in closed form.** On the rigidly rotating base state every geometric factor is exactly
periodic in the lag with period T = 2 pi/Omega, because Rot(phi_j - Omega(u + T)) = Rot(phi_j - Omega u).
Only the exponentials distinguish one revolution from the next, so for each decay rate lambda

    int_0^inf e^{-(lambda+s)u} f(u) du = [1/(1 - e^{-(lambda+s)T})] int_0^T e^{-(lambda+s)u} f(u) du,

which is the same geometric collapse stage 2 used for the self-force. The quadrature therefore covers ONE
revolution and the tail is exact rather than truncated. Parametrising by phase phi = Omega u makes the
geometry independent of Omega as well, so the equilibrium iteration costs nothing.

The base state is only quasi-stationary for the lagged rungs: a body's own wake pulls backward, so the ring
loses angular momentum slowly. That tangential force is measured and reported rather than assumed away.
"""
import numpy as np

import formation as FM

GM = FM.GM
T0 = 2*np.pi
JMAT = np.array([[0., 1.], [-1., 0.]])
MODELS = ('none', 'instantaneous', 'one_stage', 'two_stage', 'newtonian')
LAGGED = ('one_stage', 'two_stage')


class Gaussian:
    """K(d) = exp(-|d|^2/2w^2); a body writes q K and feels +grad of the accumulated field."""

    kind = 'gaussian'

    def __init__(self, w):
        self.w = float(w)

    def grad(self, d):
        w2 = self.w*self.w
        g = np.exp(-(d*d).sum(-1)/(2*w2))
        return -g[..., None]*d/w2

    def hess(self, d):
        w2 = self.w*self.w
        g = np.exp(-(d*d).sum(-1)/(2*w2))
        outer = d[..., :, None]*d[..., None, :]
        return g[..., None, None]*(outer/(w2*w2) - np.eye(2)/w2)


class NewtonianPair:
    """K(d) = 1/|d|: ordinary pairwise attraction, for the identification rung. No self term."""

    kind = 'newtonian'

    def grad(self, d):
        r2 = (d*d).sum(-1)
        r = np.sqrt(np.where(r2 > 0, r2, 1.))              # the self term is zeroed by the caller
        return np.where((r2 > 0)[..., None], -d/r[..., None]**3, 0.)

    def hess(self, d):
        r2 = (d*d).sum(-1)
        r = np.sqrt(np.where(r2 > 0, r2, 1.))
        outer = d[..., :, None]*d[..., None, :]
        h = 3*outer/r[..., None, None]**5 - np.eye(2)/r[..., None, None]**3
        return np.where((r2 > 0)[..., None, None], h, 0.)


def phase_quadrature(panels=256, order=8):
    """Gauss-Legendre panels over one revolution of phase. Uniform panels, because the passage of body i
    over body j's earlier position happens at phase phi_j, anywhere in the revolution -- clustering at the
    ends, as the single-writer self-force could afford to do, would miss every other writer's peak."""
    e = np.linspace(0., 2*np.pi, panels + 1)
    xg, wg = np.polynomial.legendre.leggauss(order)
    a, b = e[:-1, None], e[1:, None]
    return ((.5*(b - a)*xg[None, :] + .5*(a + b)).ravel(), (.5*(b - a)*wg[None, :]).ravel())


def decay_terms(model, tau_keep, tau_form):
    """[(coefficient, decay rate)] of h(u) = sum_i c_i exp(-lambda_i u); int h = tau_keep in every case."""
    if model == 'one_stage':
        return [(1., 1./tau_keep)]
    if model == 'two_stage':
        c = tau_keep/(tau_keep - tau_form)
        return [(c, 1./tau_keep), (-c, 1./tau_form)]
    raise ValueError(model)


def transfer(model, s, tau_keep, tau_form):
    """H(s), the field response to a source perturbation."""
    if model == 'none':
        return 0.*s
    if model in ('instantaneous', 'newtonian'):
        return tau_keep + 0.*s
    if model == 'one_stage':
        return tau_keep/(1 + s*tau_keep)
    return tau_keep/((1 + s*tau_keep)*(1 + s*tau_form))


class Ring:
    """N equal writers on a circle of radius R, rigidly rotating, with one response rung."""

    def __init__(self, model, n, strength, R=1., w=.2, tau_keep=10*T0, tau_form=3*T0,
                 kernel=None, panels=128, order=8, omega_iters=200, omega_tol=1e-15):
        if model not in MODELS:
            raise ValueError(model)
        self.model, self.n, self.R, self.w = model, int(n), float(R), float(w)
        self.tau_keep, self.tau_form = float(tau_keep), float(tau_form)
        self.strength = 0. if model == 'none' else float(strength)
        self.kernel = kernel if kernel is not None else (NewtonianPair() if model == 'newtonian'
                                                         else Gaussian(w))
        self.lagged = model in LAGGED
        self.terms = decay_terms(model, self.tau_keep, self.tau_form) if self.lagged else [(1., 0.)]
        self.gain = self.tau_keep if model in ('instantaneous', 'one_stage', 'two_stage') else 1.
        self.phi = 2*np.pi*np.arange(self.n)/self.n
        self.P = self.R*np.stack([np.cos(self.phi), np.sin(self.phi)], axis=1)
        self.panels, self.order = int(panels), int(order)
        self._geometry()
        self._solve_equilibrium(omega_iters, omega_tol)

    # ---------------------------------------------------------------- Omega-independent geometry
    def _geometry(self):
        if self.lagged:
            self.p, self.pw = phase_quadrature(self.panels, self.order)
        else:
            self.p, self.pw = np.zeros(1), np.array([1.])       # zero lag: one node at phase 0
        ang = self.phi[:, None] - self.p[None, :]                # phi_j - Omega u
        d = self.P[0][None, None, :] - self.R*np.stack([np.cos(ang), np.sin(ang)], axis=-1)
        hess, grad = self.kernel.hess(d), self.kernel.grad(d)
        if self.kernel.kind == 'newtonian':                       # no self interaction
            hess[0], grad[0] = 0., 0.
        c, s = np.cos(ang), np.sin(ang)
        L = np.stack([np.stack([c, -s], -1), np.stack([s, c], -1)], -2)
        HL = np.einsum('jpab,jpbc->jpac', hess, L)
        # stored per unit strength, so the coupling can be rescaled without rebuilding the geometry
        self.Z1 = self.n*np.fft.ifft(HL, axis=0)                  # Z[m, p] = sum_j e^{2 pi i m j/n} HL
        self.Hsum1 = hess.sum(0)                                  # (p, 2, 2)
        self.Gsum1 = grad.sum(0)                                  # (p, 2)
        self.D = np.diag([2*GM/self.R**3, -GM/self.R**3])

    @property
    def Z(self):
        return self.strength*self.Z1

    @property
    def Hsum(self):
        return self.strength*self.Hsum1

    @property
    def Gsum(self):
        return self.strength*self.Gsum1

    def set_strength(self, q, iters=200, tol=1e-15):
        """Rescale the coupling and re-solve the equilibrium: the geometry does not depend on it."""
        self.strength = float(q)
        self._solve_equilibrium(iters, tol)
        return self

    def _weights(self, s, omega):
        """w_p for int_0^inf h(u) e^{-su} f(u) du with f periodic, summed exactly over revolutions."""
        if not self.lagged:
            return self.gain*self.pw*np.exp(-s*0.)
        T = 2*np.pi/omega
        out = 0.
        for c, lam in self.terms:
            z = lam + s
            geo = 1. if np.real(z*T) > 700 else 1./(1. - np.exp(-z*T))
            out = out + c*geo*self.pw*np.exp(-z*self.p/omega)/omega
        return out

    # ---------------------------------------------------------------- base state
    def _memory_acceleration(self, omega):
        return np.einsum('p,pd->d', np.real(self._weights(0., omega)), self.Gsum)

    def _solve_equilibrium(self, iters, tol):
        omega = np.sqrt(GM/self.R**3)
        k = 0
        for k in range(iters):
            a = self._memory_acceleration(omega)
            new = np.sqrt(max(GM/self.R**2 - a[0], 1e-300)/self.R)   # +x is radially outward at body 0
            done = abs(new - omega) <= tol*max(1., abs(new))
            omega = new
            if done:
                break
        self.omega = float(omega)
        a = self._memory_acceleration(self.omega)
        self.a_radial_inward = float(-a[0])
        self.a_tangential = float(a[1])
        self.support_fraction = self.a_radial_inward/(GM/self.R**2)
        self.drag_timescale_periods = float(abs(self.omega*self.R/(self.a_tangential*T0))
                                            if self.a_tangential else np.inf)
        self.equilibrium_iterations = k + 1
        self.A = np.einsum('p,pab->ab', np.real(self._weights(0., self.omega)), self.Hsum)
        self.free = self.omega**2*np.eye(2) + self.D

    # ---------------------------------------------------------------- mode problem
    def B(self, m, s):
        if self.strength == 0.:
            return np.zeros((2, 2), complex)
        return np.einsum('p,pab->ab', self._weights(s, self.omega), self.Z[m % self.n])

    def M(self, m, s):
        return s*s*np.eye(2) - 2*self.omega*s*JMAT - self.free - self.A + self.B(m, s)

    def det(self, m, s):
        M = self.M(m, s)
        return M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]

    def _quartic_roots(self, B):
        """Exact roots of det[s^2 I - 2 Omega s J - free - A + B] = 0 for a frozen B."""
        C1 = -2*self.omega*JMAT
        C0 = -self.free - self.A + B
        comp = np.block([[np.zeros((2, 2)), np.eye(2)], [-C0, -C1]]).astype(complex)
        return np.linalg.eigvals(comp)

    def _newton(self, m, s, iters=40, tol=1e-12):
        for _ in range(iters):
            f = self.det(m, s)
            h = 1e-7*max(1., abs(s))
            df = (self.det(m, s + h) - self.det(m, s - h))/(2*h)
            if df == 0:
                break
            ds = f/df
            s = s - (ds if abs(ds) < .5*max(1., abs(s)) else .5*ds*max(1., abs(s))/abs(ds))
            if abs(ds) <= tol*max(1., abs(s)):
                break
        return s

    def modes(self, m, steps=12, weak=.02):
        """The four branches, continued in the PHYSICAL coupling from a weak ring to this one.

        Freezing B at the current s and re-solving the quartic -- the obvious inner iteration -- converges
        to whichever root is nearest, and at the stage 4 coupling that is the slow branch: it reported the
        two-stage m = 2 mode as neutral (e-folding 730 periods) when the true root of the same determinant
        gives 8.6 periods. Continuation in strength, with the equilibrium re-solved at every step, follows
        each branch instead of hunting for one, and the strength scan is monotonic as a result."""
        if not self.lagged:
            roots = self._quartic_roots(self.B(m, 0.))
            return roots, float(max(abs(self.det(m, s)) for s in roots))
        target = self.strength
        try:
            self.set_strength(weak*target)
            roots = [complex(z) for z in self._quartic_roots(self.B(m, 0.))]
            roots = [self._newton(m, z) for z in roots]
            for frac in np.linspace(weak, 1., steps + 1)[1:]:
                self.set_strength(frac*target)
                roots = [self._newton(m, z) for z in roots]
        finally:
            self.set_strength(target)
        roots = [self._newton(m, z) for z in roots]
        return np.array(roots), float(max(abs(self.det(m, z)) for z in roots))

    def unstable_count(self, m, s_max=None, y_max=None, points=500, eps_frac=1e-4):
        """Zeros of det in the right half strip, by the argument principle: a check that continuing the
        four branches did not miss an unstable root on another branch. The contour is offset from the
        imaginary axis by eps_frac*Omega, and the tracked roots must be counted with the SAME threshold --
        a marginal root sits on the boundary and would otherwise be counted on one side only."""
        s_max = s_max if s_max is not None else 10*self.omega
        y_max = y_max if y_max is not None else max(10*self.omega, 3*self.omega*max(m, 1))
        eps = eps_frac*self.omega
        corners = [complex(eps, -y_max), complex(s_max, -y_max), complex(s_max, y_max), complex(eps, y_max)]
        total = 0.
        for a, b in zip(corners, corners[1:] + corners[:1]):
            ts = np.linspace(0., 1., points)
            vals = np.array([self.det(m, a + (b - a)*t) for t in ts])
            total += np.sum(np.diff(np.unwrap(np.angle(vals))))
        return int(round(total/(2*np.pi)))

    def eigenvector(self, m, s):
        _, _, vh = np.linalg.svd(self.M(m, s))
        return vh[-1].conj()

    # ------------------------------------------------- direct history, for the finite-difference gate
    def direct_nodes(self, spans=25., panels_per_rev=64, order=8):
        """Lag nodes and h(u) weights WITHOUT the geometric collapse, over [0, spans*tau_keep]. Slow, and
        used only to check the collapsed formulas and to finite-difference the full history force."""
        if not self.lagged:
            return np.zeros(1), np.array([self.gain])
        T = 2*np.pi/self.omega
        revs = max(int(np.ceil(spans*self.tau_keep/T)), 1)
        e = np.linspace(0., revs*T, revs*panels_per_rev + 1)
        xg, wg = np.polynomial.legendre.leggauss(order)
        a, b = e[:-1, None], e[1:, None]
        u = (.5*(b - a)*xg[None, :] + .5*(a + b)).ravel()
        w = (.5*(b - a)*wg[None, :]).ravel()
        h = sum(c*np.exp(-lam*u) for c, lam in self.terms)
        return u, h*w

    def force_direct(self, nodes=None, mode=None, chunk=8192):
        """Exact history force on body 0, in its own radial/tangential frame. `mode` is (m, s, z, eps):
        body j's local-frame displacement at lag u is eps*Re[z exp(2 pi i m j/N) exp(-s u)]. The separation

            d_j(u) = (P_0 + u_0(0)) - Rot(-Omega u) (P_j + Lambda_j u_j(-u))

        is the same expression the linearisation expands, evaluated here without expanding it."""
        u, hw = nodes if nodes is not None else self.direct_nodes()
        lam_j = np.stack([np.stack([np.cos(self.phi), -np.sin(self.phi)], -1),
                          np.stack([np.sin(self.phi), np.cos(self.phi)], -1)], -2)       # Lambda_j
        here = self.P[0].copy()
        if mode is not None:
            m, s, z, eps = mode
            here = here + eps*np.real(np.asarray(z))
        total = np.zeros(2)
        for k in range(0, len(u), chunk):
            uu, ww = u[k:k + chunk], hw[k:k + chunk]
            ang = self.phi[:, None] - self.omega*uu[None, :]
            pos = self.R*np.stack([np.cos(ang), np.sin(ang)], axis=-1)                   # (n, U, 2)
            if mode is not None:
                bloch = np.exp(2j*np.pi*m*np.arange(self.n)/self.n)
                loc = eps*np.real(np.asarray(z)[None, None, :]*bloch[:, None, None]
                                  * np.exp(-s*uu)[None, :, None])
                world = np.einsum('jab,jub->jua', lam_j, loc)
                c_, s_ = np.cos(self.omega*uu), np.sin(self.omega*uu)
                back = np.stack([np.stack([c_, s_], -1), np.stack([-s_, c_], -1)], -2)   # Rot(-Omega u)
                pos = pos + np.einsum('uab,jub->jua', back, world)
            g = self.kernel.grad(here[None, None, :] - pos)
            if self.kernel.kind == 'newtonian':
                g = g.copy()
                g[0] = 0.
            total = total + self.strength*np.einsum('u,jud->d', ww, g)
        return total

    def base_field_gradient(self, x, nodes=None):
        """grad C of the equilibrium field at an arbitrary rotating-frame point, for the Hessian check."""
        u, hw = nodes if nodes is not None else self.direct_nodes()
        ang = self.phi[:, None] - self.omega*u[None, :]
        base = self.R*np.stack([np.cos(ang), np.sin(ang)], axis=-1)
        g = self.kernel.grad(np.asarray(x, float)[None, None, :] - base)
        if self.kernel.kind == 'newtonian':
            g = g.copy()
            g[0] = 0.
        return self.strength*np.einsum('u,jud->d', hw, g)


def growth_and_pattern(s, omega, m):
    """Amplitude e-folding time in reference periods, and inertial pattern speed."""
    gr = float(np.real(s))
    return dict(growth_rate=gr,
                e_folding_periods=float(1./(gr*T0)) if gr > 1e-14 else float('inf'),
                # None, not NaN: the archive is regression-compared and NaN never equals itself
                pattern_speed=float(omega - np.imag(s)/m) if m else None,
                pattern_speed_over_omega=float((omega - np.imag(s)/m)/omega) if m else None,
                frequency_rotating_frame=float(np.imag(s)))


def spectrum(ring, m_values=None, count_check=True):
    """The most unstable root for each angular mode."""
    m_values = list(m_values if m_values is not None else range(ring.n//2 + 1))
    rows, worst_residual = {}, 0.
    for m in m_values:
        roots, residual = ring.modes(m)
        worst_residual = max(worst_residual, residual)
        k = int(np.argmax(roots.real))
        row = dict(m=int(m), roots=[[float(z.real), float(z.imag)] for z in roots],
                   residual=residual, **growth_and_pattern(roots[k], ring.omega, m))
        if count_check:
            tracked = int(np.sum(roots.real > 1e-4*ring.omega))
            row.update(unstable_roots_tracked=tracked,
                       unstable_roots_in_half_plane=ring.unstable_count(m))
            row['count_matches'] = bool(row['unstable_roots_in_half_plane'] == tracked)
        rows[int(m)] = row
    return dict(model=ring.model, writers=ring.n, omega=ring.omega,
                support_fraction=ring.support_fraction,
                tangential_acceleration=ring.a_tangential,
                drag_timescale_periods=ring.drag_timescale_periods,
                worst_residual=worst_residual,
                fastest=max(rows.values(), key=lambda r: r['growth_rate']),
                all_counts_match=all(r.get('count_matches', True) for r in rows.values()),
                modes=rows)


# ---------------------------------------------------------------- nonlinear verification
class InstantaneousField(FM.MemoryField):
    """C = tau_keep * S at every instant: the owner's matched zero-delay control. Same kernel, same
    static gain, no lag. It is NOT the memory-off control, which removes the attraction altogether."""

    def __init__(self, half_width, spacing, w, tau_keep, tau_form=0., reach=5.):
        super().__init__(half_width, spacing, w, tau_keep, 0., reach)

    def advance(self, positions, rates, h):
        self.C = self.tau_keep*self._source(positions, rates)


def run_ring(positions, velocities, rates, w, tau_keep, t_end, h_max, tau_form=0., field_cls=None,
             half_width=2.5, spacing=None, prime_ring_R=None, r_min=.25, eta=.005, samples=400,
             snapshot_every=None):
    """formation.run's step, with the field class swappable. G2 checks it against formation.run itself."""
    spacing = spacing if spacing is not None else w/5
    cls = field_cls if field_cls is not None else FM.MemoryField
    field = cls(half_width, spacing, w, tau_keep, tau_form)
    instant = isinstance(field, InstantaneousField)
    if prime_ring_R is not None:
        field.prime_with_ring(prime_ring_R, A=float(np.sum(rates))*tau_keep)
    x = np.array(positions, float)
    v = np.array(velocities, float)
    rates = np.asarray(rates, float)
    if instant:
        field.advance(x, rates, 0.)
    rec = dict(t=[], r=[], theta=[], L=[])
    snaps, next_snap = [], 0.
    accel = lambda xx: (-GM*xx/np.maximum(np.linalg.norm(xx, axis=1, keepdims=True), 1e-12)**3
                        + field.sample_gradient(xx))
    a = accel(x)
    t, status, next_sample = 0., FM.COMPLETED, 0.
    while t < t_end:
        r = np.linalg.norm(x, axis=1)
        if r.min() < r_min:
            status = FM.UNRESOLVED_CENTRE
            break
        h = min(h_max, eta*float(np.min(np.sqrt(r**3/GM))), t_end - t)
        v_half = v + .5*h*a
        x_mid = x + .5*h*v_half
        x = x + h*v_half
        if np.max(np.abs(x)) > .9*half_width:
            status = FM.LEFT_DOMAIN
            break
        field.advance(x if instant else x_mid, rates, h)
        a = accel(x)
        v = v_half + .5*h*a
        t += h
        if snapshot_every and (t >= next_snap or t >= t_end):
            next_snap = t + snapshot_every
            snaps.append((t, x.copy(), v.copy()))
        if t >= next_sample or t >= t_end:
            next_sample = t + t_end/samples
            rec['t'].append(t)
            rec['r'].append(np.linalg.norm(x, axis=1))
            rec['theta'].append(np.arctan2(x[:, 1], x[:, 0]))
            rec['L'].append(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])
    out = {k: np.array(val) for k, val in rec.items()}
    out.update(status=status, t_final=t, positions=x, velocities=v, snapshots=snaps)
    return out


def seed_mode(ring, m, amplitude, s=None, vector=None):
    """Base ring plus one linear eigenmode, in positions AND velocities."""
    if s is None:
        roots, _ = ring.modes(m)
        s = roots[int(np.argmax(roots.real))]
    vec = ring.eigenvector(m, s) if vector is None else vector
    vec = vec/np.linalg.norm(vec)
    phase = np.exp(2j*np.pi*m*np.arange(ring.n)/ring.n)
    local = amplitude*np.real(vec[None, :]*phase[:, None])
    dlocal = amplitude*np.real(s*vec[None, :]*phase[:, None])
    er = np.stack([np.cos(ring.phi), np.sin(ring.phi)], axis=1)
    et = np.stack([-np.sin(ring.phi), np.cos(ring.phi)], axis=1)
    dx = local[:, :1]*er + local[:, 1:]*et
    x = ring.P + dx
    v = (ring.omega*ring.R*et + dlocal[:, :1]*er + dlocal[:, 1:]*et
         + ring.omega*np.stack([-dx[:, 1], dx[:, 0]], axis=1))
    return x, v, complex(s), vec


def measure_mode(out, m, R, t_lo, t_hi, period=T0, prony=True):
    """Growth rate and rotating-frame frequency of mode m from the recorded per-body radii.

    A seeded eigenmode never arrives alone: the same initial condition projects onto the stable epicyclic
    branch too, and a straight-line fit to log|c| of a beating signal reports neither rate. So the fit is
    also done with a two-term Prony estimate, which solves for both complex frequencies from the data with
    no knowledge of the prediction, and reports the faster-growing one."""
    t, r = out['t'], out['r']
    sel = (t >= t_lo*period) & (t <= t_hi*period)
    if sel.sum() < 16:
        return None
    n = r.shape[1]
    c = np.mean((r[sel] - R)*np.exp(-2j*np.pi*m*np.arange(n)[None, :]/n), axis=1)
    ts, good = t[sel], np.abs(c) > 0
    if good.sum() < 16:
        return None
    ts, c = ts[good], c[good]
    row = dict(growth_rate=float(np.polyfit(ts, np.log(np.abs(c)), 1)[0]),
               frequency_rotating_frame=float(np.polyfit(ts, np.unwrap(np.angle(c)), 1)[0]),
               amplitude_start=float(np.abs(c[0])), amplitude_end=float(np.abs(c[-1])),
               samples=int(len(ts)))
    if prony:
        dt = float(np.median(np.diff(ts)))
        u = np.arange(ts[0], ts[-1] + .5*dt, dt)
        z = np.interp(u, ts, c.real) + 1j*np.interp(u, ts, c.imag)
        A = np.stack([z[1:-1], z[:-2]], axis=1)
        coef, *_ = np.linalg.lstsq(A, z[2:], rcond=None)
        roots = np.roots([1., -coef[0], -coef[1]])
        ss = np.log(roots)/dt
        k = int(np.argmax(ss.real))
        row.update(prony_growth_rate=float(ss[k].real), prony_frequency=float(ss[k].imag),
                   prony_second=[float(ss[1 - k].real), float(ss[1 - k].imag)])
    return row

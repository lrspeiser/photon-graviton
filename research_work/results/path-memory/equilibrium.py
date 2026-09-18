"""RUT-1 stage 6: a supported orbital population constructed in self-consistent equilibrium (R3).

Stage 5 settled that a cold ring is collectively unstable with or without memory, so the question is no
longer how to keep that ring intact. It is which populations this response can support. Adding random
speeds to a ring is not an answer: that is another out-of-equilibrium state whose initial adjustment
hides the question. What follows constructs a stationary state instead.

In an axisymmetric potential the orbital energy and the angular momentum are constant along every orbit,
so any non-negative function of them is a steady state of the collisionless dynamics. The declared family
is a warm annulus,

    f(E, L) = exp[-(L - L0)^2 / 2 dL^2] exp[-(E - E_circ(L))^2 / 2 dE^2],

where E_circ(L) is the energy of the circular orbit of angular momentum L IN THE CURRENT POTENTIAL. dL
sets the width of the annulus, dE the spread of epicyclic amplitudes -- the "temperature". Both factors
are Gaussian, so the mass is finite and the tapering is explicit, with no need for a separate cut.

The population and the field it writes are solved together. For writing proportional to mass, at rate
alpha per unit mass, an axisymmetric surface density sigma(R) writes

    C0(r) = tau_keep alpha integral 2 pi R' sigma(R') exp[-(r-R')^2/2w^2] I0e(r R'/w^2) dR',

which is the mature-ring field of stage 1 integrated over radius: the Gaussian kernel's angular average.
The total potential is Phi_ext - C0, and E_circ(L) and the density follow from it, so the three close on
each other and are iterated to a fixed point.

Existence and stability are different questions. This module answers only the first, and supplies the
sampler that lets the second be tested.
"""
import numpy as np
from scipy.special import i0e, i1e

GM = 1.
T0 = 2*np.pi


class WarmAnnulus:
    """f(E, L) = exp[-(L-L0)^2/2 dL^2] exp[-(E - E_circ(L))^2/2 dE^2], with its self-consistent field."""

    def __init__(self, L0=1., dL=.08, dE=.02, w=.2, tau_keep=10*T0, alpha=1.,
                 r_lo=.3, r_hi=3., n_r=280, n_v=121, v_max=3.):
        self.L0, self.dL, self.dE = float(L0), float(dL), float(dE)
        self.w, self.tau_keep, self.alpha = float(w), float(tau_keep), float(alpha)
        self.r = np.linspace(r_lo, r_hi, n_r)
        self.n_v, self.v_max = int(n_v), float(v_max)
        self.C = np.zeros_like(self.r)                  # the written field on the radial grid
        self.sigma = np.zeros_like(self.r)
        self.history = []

    # ---------------------------------------------------------------- potential and its circular orbits
    def phi_total(self, r=None, C=None):
        r = self.r if r is None else r
        C = self.C if C is None else C
        return -GM/np.maximum(r, 1e-12) - C

    def _dphi(self, r, C=None, sigma=None):
        """dPhi/dr on the grid: the exact Kepler part plus the ANALYTIC derivative of the written part.

        A centred difference of C is second-order accurate, which leaves about 1e-3 of relative noise on
        this grid -- and a cold distribution amplifies it, because that noise enters E_circ(L) and dE is
        only 0.01, so a 1e-3 error moves the energy Gaussian by ten percent of its own width. The fixed
        point then stalls at 1e-3 however hard it is under-relaxed. Differentiating the kernel instead
        removes the floor: only I0e and I1e are needed, since d/du [exp(-u) I0(u)] = exp(-u)[I1(u)-I0(u)].
        """
        sigma = self.sigma if sigma is None else sigma
        return GM/r**2 - self.write_field_gradient(sigma, r)

    def circular_energy(self, L, r, C, sigma=None):
        """E_circ(L): energy of the circular orbit with angular momentum L in the current potential.

        The circular radius solves L^2/r^3 = dPhi/dr. Solved by interpolating the monotone map r -> L_c(r)
        rather than root-finding per L, which is both faster and single-valued where the map is monotone.
        """
        Lc = np.sqrt(np.maximum(r**3*self._dphi(r, C, sigma), 0.))
        good = np.isfinite(Lc) & (Lc > 0)
        rr, Ll = r[good], Lc[good]
        order = np.argsort(Ll)
        r_of_L = np.interp(L, Ll[order], rr[order])
        phi = np.interp(r_of_L, r, self.phi_total(r, C))
        return .5*(L/r_of_L)**2 + phi, r_of_L

    # ---------------------------------------------------------------- the density of the declared f
    def _grids(self, C, n_L=160, n_vr=96, reach=5., sigma=None):
        """Integration nodes in the conserved variables, sized by the distribution's own widths.

        A uniform (v_r, v_t) grid cannot do this: with dE = 0.01 the energy Gaussian is narrower than one
        velocity cell, so the density it returns is quantisation noise and the fixed point it converges to
        is the discretisation's, not the model's. Integrating over L and v_r instead puts both Gaussians
        on the integration variables, because E = E_min(L,R) + v_r^2/2 makes dE/|v_r| exactly dv_r -- the
        turning-point singularity of the energy integral is removed rather than resolved.
        """
        L = np.linspace(self.L0 - reach*self.dL, self.L0 + reach*self.dL, n_L)
        Ec, _ = self.circular_energy(L, self.r, C, sigma)
        xg, wg = np.polynomial.legendre.leggauss(n_vr)
        return L, Ec, xg, wg

    def surface_density(self, C, n_L=160, n_vr=96, reach=5., sigma=None):
        L, Ec, xg, wg = self._grids(C, n_L, n_vr, reach, sigma)
        phi = self.phi_total(self.r, C)
        dL_w = L[1] - L[0]
        out = np.empty_like(self.r)
        for i, R in enumerate(self.r):
            E_min = phi[i] + .5*(L/R)**2
            hi = 2*(Ec + reach*self.dE - E_min)
            lo = np.maximum(2*(Ec - reach*self.dE - E_min), 0.)
            live = hi > lo
            if not live.any():
                out[i] = 0.
                continue
            a, b = np.sqrt(lo[live]), np.sqrt(hi[live])
            vr = .5*(b - a)[:, None]*xg[None, :] + .5*(a + b)[:, None]
            w = .5*(b - a)[:, None]*wg[None, :]
            E = E_min[live][:, None] + .5*vr*vr
            f = np.exp(-(L[live][:, None] - self.L0)**2/(2*self.dL**2)
                       - (E - Ec[live][:, None])**2/(2*self.dE**2))
            out[i] = 2.*float(np.sum(w*f))*dL_w/R
        return out

    def write_field_gradient(self, sigma, r=None, alpha=None):
        """dC0/dr, differentiating the ring kernel rather than the result."""
        r = self.r if r is None else r
        R = self.r
        weight = 2*np.pi*R*sigma*(self.alpha if alpha is None else alpha)*self.tau_keep
        d = r[:, None] - R[None, :]
        u = r[:, None]*R[None, :]/self.w**2
        g = np.exp(-d*d/(2*self.w*self.w))
        dk = g*(-d/self.w**2*i0e(u) + (R[None, :]/self.w**2)*(i1e(u) - i0e(u)))
        return np.trapezoid(dk*weight[None, :], R, axis=1)

    def write_field(self, sigma):
        """C0(r) from the axisymmetric kernel: the stage 1 ring field integrated over radius."""
        R = self.r
        weight = 2*np.pi*R*sigma*self.alpha*self.tau_keep
        kern = np.exp(-(R[:, None] - R[None, :])**2/(2*self.w**2))*i0e(R[:, None]*R[None, :]/self.w**2)
        return np.trapezoid(kern*weight[None, :], R, axis=1)

    # ---------------------------------------------------------------- the fixed point
    def _map(self, sigma):
        """One pass of the coupled equations as a map on the density alone: sigma -> field -> sigma."""
        C = self.write_field(sigma)
        return self.surface_density(C, sigma=sigma), C

    def fixed_point(self, iters=300, tol=1e-11, damping=.3, memory=6):
        """Anderson-accelerated solve, retried with other mixing strengths if one stalls. Which strength
        converges depends on the annulus -- wide cold ones want strong mixing, wide warm ones weak -- so a
        short ladder is tried in a fixed order and the first converged result is kept."""
        start = self.sigma.copy()
        for beta in (damping, .5, .15):
            self.sigma = start.copy()
            self._anderson(iters, tol, beta, memory)
            if self.residual < tol:
                break
        self.mixing_used = beta
        return self

    def _anderson(self, iters, tol, damping, memory):
        """Solve sigma = G(sigma) at the current alpha, with Anderson acceleration.

        Plain under-relaxed iteration converges to 1e-11 for narrow annuli but is not reliably contractive
        for wide ones: depending on the damping and the starting point it either converges or stalls near
        1e-3, which is the behaviour of a map with an eigenvalue close to the unit circle. Anderson mixing
        uses the last few residuals to cancel those slow directions, and needs no derivative of the map.
        """
        R = self.r
        x = self.sigma.copy()
        if not x.any():                                   # a cold start: the Kepler-only population
            x = self.surface_density(np.zeros_like(R), sigma=np.zeros_like(R))
        X, F = [], []
        change = np.inf
        for k in range(iters):
            g, C = self._map(x)
            f = g - x
            scale = max(float(np.max(np.abs(g))), 1e-300)
            change = float(np.max(np.abs(f))/scale)
            self.history.append(change)
            if change < tol:
                x = g
                break
            X.append(x.copy())
            F.append(f.copy())
            X, F = X[-(memory + 1):], F[-(memory + 1):]
            if len(F) > 1:
                dF = np.stack([F[i + 1] - F[i] for i in range(len(F) - 1)], axis=1)
                dX = np.stack([X[i + 1] - X[i] for i in range(len(X) - 1)], axis=1)
                gamma, *_ = np.linalg.lstsq(dF, f, rcond=1e-10)
                x_new = x + damping*f - (dX + damping*dF)@gamma
            else:
                x_new = x + damping*f
            x = np.maximum(x_new, 0.)                     # a density cannot go negative
        self.sigma = x
        self.C = self.write_field(x)
        self.iterations = k + 1
        self.residual = change
        return self

    def solve(self, iters=300, tol=1e-11, damping=.3, target_support=None, outer=8, outer_tol=1e-7):
        """Solve at fixed alpha, then, if a support is targeted, secant on alpha around that solve.

        Rescaling alpha INSIDE the fixed point mixes two different couplings and stalls: the population
        chases a field whose strength is being changed underneath it. Separating them lets the inner map
        converge to machine precision and leaves the outer problem one-dimensional.
        """
        if target_support is None:
            return self.fixed_point(iters, tol, damping)
        def support_at(alpha):
            # warm start: the field scales with alpha to first order, so carry the last solution across
            self.alpha = float(alpha)                     # sigma carries over: a warm start
            self.fixed_point(iters, tol, damping)
            i = int(np.argmax(self.sigma))
            return float(-self.write_field_gradient(self.sigma)[i]/(GM/self.r[i]**2))
        a0, a1 = self.alpha, self.alpha*1.5
        s0, s1 = support_at(a0), support_at(a1)
        self.outer_history = [(a0, s0), (a1, s1)]
        for _ in range(outer):
            if abs(s1 - target_support) < outer_tol or s1 == s0:
                break
            a2 = a1 + (target_support - s1)*(a1 - a0)/(s1 - s0)
            a0, s0 = a1, s1
            a1 = float(np.clip(a2, 1e-6, 1e6))
            s1 = support_at(a1)
            self.outer_history.append((a1, s1))
        self.support_error = abs(s1 - target_support)
        return self

    # ---------------------------------------------------------------- what the state is
    def summary(self):
        R, s = self.r, self.sigma
        mass = float(np.trapezoid(2*np.pi*R*s, R))
        mean_r = float(np.trapezoid(2*np.pi*R*s*R, R)/mass) if mass else float('nan')
        var_r = float(np.trapezoid(2*np.pi*R*s*(R - mean_r)**2, R)/mass) if mass else float('nan')
        i = int(np.argmax(s))
        inward = -self.write_field_gradient(self.sigma)
        half = s.max()/2
        above = R[s >= half]
        return dict(mass=mass, mean_radius=mean_r, rms_width=float(np.sqrt(max(var_r, 0.))),
                    peak_radius=float(R[i]), peak_density=float(s[i]),
                    fwhm=float(above.max() - above.min()) if len(above) > 1 else 0.,
                    support_fraction_at_peak=float(inward[i]/(GM/R[i]**2)),
                    support_fraction_max=float(np.max(inward/(GM/R**2))),
                    alpha=self.alpha, C_peak=float(self.C.max()),
                    iterations=self.iterations, residual=self.residual,
                    edge_density_ratio=float(max(s[0], s[-1])/max(s.max(), 1e-300)))

    def dispersions(self, R=None, n_L=160, n_vr=96, reach=5.):
        """Radial and azimuthal velocity dispersion of the declared f at a radius, on the same nodes."""
        R = float(self.r[int(np.argmax(self.sigma))]) if R is None else float(R)
        L, Ec, xg, wg = self._grids(self.C, n_L, n_vr, reach)
        phi = float(np.interp(R, self.r, self.phi_total(self.r, self.C)))
        E_min = phi + .5*(L/R)**2
        hi = 2*(Ec + reach*self.dE - E_min)
        lo = np.maximum(2*(Ec - reach*self.dE - E_min), 0.)
        live = hi > lo
        if not live.any():
            return dict(radius=R, sigma_r=float('nan'), sigma_t=float('nan'), v_mean=float('nan'))
        a, b = np.sqrt(lo[live]), np.sqrt(hi[live])
        vr = .5*(b - a)[:, None]*xg[None, :] + .5*(a + b)[:, None]
        w = .5*(b - a)[:, None]*wg[None, :]
        E = E_min[live][:, None] + .5*vr*vr
        f = w*np.exp(-(L[live][:, None] - self.L0)**2/(2*self.dL**2)
                     - (E - Ec[live][:, None])**2/(2*self.dE**2))
        vt = (L[live]/R)[:, None]*np.ones_like(vr)
        tot = float(f.sum())
        if tot <= 0:
            return dict(radius=R, sigma_r=float('nan'), sigma_t=float('nan'), v_mean=float('nan'))
        vt_mean = float((f*vt).sum()/tot)
        return dict(radius=R, sigma_r=float(np.sqrt((f*vr*vr).sum()/tot)),
                    sigma_t=float(np.sqrt((f*(vt - vt_mean)**2).sum()/tot)), v_mean=vt_mean)

    def rate_per_body(self, n):
        """Writing rate of each of n equal-mass bodies carrying the population's total mass."""
        return self.alpha*self.summary()['mass']/n

    def prime(self, field):
        """Put the constructed axisymmetric field on a MemoryField's Cartesian grid, value and gradient,
        with the excitation in equilibrium with it -- the general-profile version of prime_with_ring."""
        X, Y = np.meshgrid(field.axis, field.axis, indexing='ij')
        rr = np.maximum(np.hypot(X, Y), 1e-12)
        dC = self.write_field_gradient(self.sigma)
        val = np.interp(rr, self.r, self.C, left=0., right=0.)
        der = np.interp(rr, self.r, dC, left=0., right=0.)
        field.C[0] = val
        field.C[1], field.C[2] = der*X/rr, der*Y/rr
        if field.E is not None:
            field.E[:] = field.C/field.tau_keep
        return field

    def consistency_residual(self):
        """H5: recompute both sides of the coupled equations from the converged state."""
        sigma = self.surface_density(self.C, sigma=self.sigma)
        C = self.write_field(sigma)
        return dict(density_relative=float(np.max(np.abs(sigma - self.sigma))
                                           / max(np.max(np.abs(sigma)), 1e-300)),
                    field_relative=float(np.max(np.abs(C - self.C))/max(np.max(np.abs(C)), 1e-300)))

    # ---------------------------------------------------------------- sampling, for the nonlinear tests
    def sample(self, n, seed=0, n_L=160, n_vr=96, reach=5.):
        """Draw n bodies from f by inverse transform on the same adapted nodes, weighted by 2 pi R.

        Rejection sampling in a uniform velocity box has the problem the density integral had: almost
        every proposal falls where f is negligible, and the acceptance depends on a running maximum.
        """
        rng = np.random.default_rng(seed)
        L, Ec, xg, wg = self._grids(self.C, n_L, n_vr, reach)
        phi = self.phi_total(self.r, self.C)
        cells, weights = [], []
        for i, R in enumerate(self.r):
            E_min = phi[i] + .5*(L/R)**2
            hi = 2*(Ec + reach*self.dE - E_min)
            lo = np.maximum(2*(Ec - reach*self.dE - E_min), 0.)
            live = np.nonzero(hi > lo)[0]
            if not len(live):
                continue
            a, b = np.sqrt(lo[live]), np.sqrt(hi[live])
            vr = .5*(b - a)[:, None]*xg[None, :] + .5*(a + b)[:, None]
            w = .5*(b - a)[:, None]*wg[None, :]
            E = E_min[live][:, None] + .5*vr*vr
            f = w*np.exp(-(L[live][:, None] - self.L0)**2/(2*self.dL**2)
                         - (E - Ec[live][:, None])**2/(2*self.dE**2))
            f = f*2*np.pi*R                      # the area element of the ring
            for a_, b_ in zip(*np.nonzero(f > 0)):
                cells.append((R, L[live][a_], vr[a_, b_]))
                weights.append(f[a_, b_])
        cells = np.array(cells)
        p = np.array(weights)
        p = p/p.sum()
        pick = rng.choice(len(p), size=n, p=p)
        R = cells[pick, 0]
        Lp = cells[pick, 1]
        vr = cells[pick, 2]*rng.choice([-1., 1.], size=n)
        vt = Lp/R
        th = rng.uniform(0., 2*np.pi, n)
        c, s_ = np.cos(th), np.sin(th)
        x = np.stack([R*c, R*s_], axis=1)
        v = np.stack([vr*c - vt*s_, vr*s_ + vt*c], axis=1)
        return x, v

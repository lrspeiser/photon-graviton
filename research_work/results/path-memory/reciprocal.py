"""RUT-1 stage 4R: a reciprocal, energy-accounted form of the two-stage response (protocol-rut1.md).

With writing q_i = alpha m_i and rho = sum_i m_i delta(x - X_i), the two-stage equation

    tau_form C_tt + gamma C_t + C/tau_keep = alpha K * rho,        gamma = 1 + tau_form/tau_keep,

follows from a field h with C = W * h, where the Gaussian kernel is a symmetric convolution square
W * W = K:

    tau_form h_tt + gamma h_t + h/tau_keep = alpha W * rho,         matter acceleration = grad(W * h).

Deposition (W * rho) and force readout (grad W * h) now come from ONE interaction term, and for a
time-independent external potential

    H = sum_i m_i [|v_i|^2/2 + Phi_ext(X_i) - C(X_i)] + tau_form/(2 alpha) int h_t^2 + 1/(2 alpha tau_keep) int h^2

obeys dH/dt = -(gamma/alpha) int h_t^2 <= 0 exactly: a derived dissipation rate, not a variable added
afterwards to absorb whatever energy went missing.

The implementation is spectral on a periodic box, so deposition and readout use the SAME finite mode set
and are exactly reciprocal in the truncated system; the only discretization left is time stepping. Each
mode's damped oscillator is advanced exactly for a source held at the step midpoint, and matter by the same
velocity-Verlet step as the grid model.

It assumes writing proportional to mass and the symmetric Gaussian kernel. It does not identify the
physical reservoir that receives the dissipated power, it provides no spatial causality, and it does not
finish the momentum accounting for that reservoir and the fixed centre.
"""
import numpy as np

GM = 1.


class ReciprocalField:
    """h and h_t as Fourier coefficients on a periodic square box of side L centred on the origin."""

    def __init__(self, L, n_modes, w, tau_keep, tau_form, alpha):
        self.L, self.n, self.w = float(L), int(n_modes), float(w)
        self.tau_keep, self.tau_form, self.alpha = float(tau_keep), float(tau_form), float(alpha)
        self.gamma = 1 + tau_form/tau_keep
        k1 = 2*np.pi*np.fft.fftfreq(self.n, d=L/self.n)
        self.kx, self.ky = np.meshgrid(k1, k1, indexing='ij')
        k2 = self.kx**2 + self.ky**2
        # continuous 2-D transforms: K(x) = exp(-|x|^2/2w^2) -> 2 pi w^2 exp(-w^2 k^2/2); W = sqrt of it
        self.What = np.sqrt(2*np.pi)*self.w*np.exp(-self.w**2*k2/4)
        self.area = L*L
        self.h = np.zeros((self.n, self.n), complex)
        self.hd = np.zeros((self.n, self.n), complex)
        disc = np.sqrt(self.gamma**2 - 4*self.tau_form/self.tau_keep + 0j)
        self.lp = (-self.gamma + disc)/(2*self.tau_form)       # -1/tau_keep
        self.lm = (-self.gamma - disc)/(2*self.tau_form)       # -1/tau_form

    def _phases(self, X):
        """exp(-i k . X_i) for every body and mode: shape (N, n, n)."""
        return np.exp(-1j*(self.kx[None]*X[:, 0, None, None] + self.ky[None]*X[:, 1, None, None]))

    def source(self, X, m):
        """Fourier-series coefficients of alpha (W * rho): (alpha/A) sum_i m_i W^(k) exp(-i k.X_i)."""
        ph = self._phases(X)
        return (self.alpha/self.area)*self.What*np.tensordot(m, ph, axes=1)

    def advance(self, X_mid, m, h):
        """Each mode's damped oscillator, exactly, for the source held at the midpoint positions."""
        s = self.source(X_mid, m)
        y0 = self.h - self.tau_keep*s
        yd0 = self.hd
        lp, lm = self.lp, self.lm
        b = (yd0 - lp*y0)/(lm - lp)
        a = y0 - b
        ep, em = np.exp(lp*h), np.exp(lm*h)
        self.h = a*ep + b*em + self.tau_keep*s
        self.hd = a*lp*ep + b*lm*em

    def _readout(self, X, coeff):
        """sum_k W^(k) coeff_k exp(+i k.X) and its gradient, at body positions."""
        ph = np.conj(self._phases(X))
        g = self.What*coeff
        val = np.real(np.tensordot(ph, g, axes=([1, 2], [0, 1])))
        gx = np.real(np.tensordot(ph, 1j*self.kx*g, axes=([1, 2], [0, 1])))
        gy = np.real(np.tensordot(ph, 1j*self.ky*g, axes=([1, 2], [0, 1])))
        return val, np.stack([gx, gy], axis=1)

    def C_and_grad(self, X):
        return self._readout(X, self.h)

    def dC_dt(self, X):
        return self._readout(X, self.hd)[0]

    def field_energy(self):
        """tau_form/(2 alpha) int h_t^2 + 1/(2 alpha tau_keep) int h^2, by Parseval on the box."""
        return (self.area*(self.tau_form/(2*self.alpha)*np.sum(np.abs(self.hd)**2)
                           + 1/(2*self.alpha*self.tau_keep)*np.sum(np.abs(self.h)**2)))

    def dissipation_rate(self):
        """(gamma/alpha) int h_t^2 >= 0."""
        return float(self.area*self.gamma/self.alpha*np.sum(np.abs(self.hd)**2))


def run_free(positions, velocities, masses, field, t_end, h, samples=200, checkpoints=(), period=2*np.pi):
    """Matter by velocity-Verlet in the external point-mass potential plus grad C; field advanced exactly per
    mode for the midpoint source. Steps land exactly on checkpoint times, as the grid integrator's do, so the
    two can be compared state for state. Returns the balance H(t) - H(0) + int dissipation, which must vanish."""
    x = np.array(positions, float)
    v = np.array(velocities, float)
    m = np.asarray(masses, float)

    def accel(xx):
        r = np.linalg.norm(xx, axis=1, keepdims=True)
        return -GM*xx/r**3 + field.C_and_grad(xx)[1]

    def H(xx, vv):
        C, _ = field.C_and_grad(xx)
        return float(np.sum(m*(.5*np.sum(vv*vv, axis=1) - GM/np.linalg.norm(xx, axis=1) - C))
                     + field.field_energy())

    a = accel(x)
    H0 = H(x, v)
    D = 0.
    d_prev = field.dissipation_rate()
    t = 0.
    rows, cps = [], {}
    scale = float(np.sum(m*.5*np.sum(v*v, axis=1)))
    cp_times = [c*period for c in checkpoints]
    next_sample = t_end/samples
    while t < t_end - 1e-12:
        step = min(h, t_end - t)
        for ct in cp_times:
            if t < ct - 1e-12 < t + step:
                step = ct - t
        v_half = v + .5*step*a
        x_mid = x + .5*step*v_half
        x = x + step*v_half
        field.advance(x_mid, m, step)
        a = accel(x)
        v = v_half + .5*step*a
        d_now = field.dissipation_rate()
        D += .5*step*(d_prev + d_now)
        d_prev = d_now
        t += step
        if t >= next_sample - 1e-12 or t >= t_end - 1e-12:
            next_sample = t + t_end/samples
            Hn = H(x, v)
            rows.append(dict(t=float(t), H_change=float(Hn - H0), dissipated=float(D),
                             balance=float(Hn - H0 + D), relative_to_kinetic=float(abs(Hn - H0 + D)/scale)))
        for ct in cp_times:
            if abs(t - ct) < 1e-9:
                cps[f'{ct/period:g}'] = dict(positions=x.tolist(), velocities=v.tolist())
    return dict(rows=rows, checkpoints=cps, final_positions=x.tolist(), final_velocities=v.tolist(),
                worst_balance_relative=float(max(r['relative_to_kinetic'] for r in rows)),
                total_dissipated=float(D), kinetic_scale=scale)

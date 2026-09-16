"""Semiclassical companion modes on a homogeneous propagation field (BRIDGE-1B, protocol part A).

Code units, hbar = c = 1. The field n(t) is classical; the companion field is a free quantum field in a
Gaussian state, started in the adiabatic vacuum and carried in the adiabatic (Bogoliubov) basis:

    omega_k^2 = k^2 + m0^2 + g^2 (n - n_star)^2,          theta_k' = omega_k
    alpha_k'  = (omega_k'/2 omega_k) e^{+2 i theta} beta_k
    beta_k'   = (omega_k'/2 omega_k) e^{-2 i theta} alpha_k
    K n''     = A/n^2 - V'(n) - g^2 (n - n_star) S,
    S         = sum_j w_j (|beta_j|^2 + Re(alpha_j beta_j^* e^{-2 i theta_j}))/omega_j

with isotropic quadrature weights w_j = k_j^2 dk / (2 pi^2). The adiabatic particle number is |beta_k|^2
and the companions' energy density is sum_j w_j omega_j |beta_j|^2, so the zero-point energy sum omega_k/2
and its force are subtracted: V(n) is the renormalized potential. Then

    E = K n'^2/2 + V(n) + A/n + rho_C

is conserved exactly by these equations, which is the numerical check V2.
"""
import numpy as np
from scipy.integrate import solve_ivp


def grid(kmax, J):
    """Midpoint grid in |k| with isotropic weights."""
    k = (np.arange(J) + .5)*kmax/J
    return k, k**2*(kmax/J)/(2*np.pi**2)


def landau_zener(k, g, ndot, m0=0.):
    """Occupation after one linear crossing of the mass minimum, exact for a constant crossing speed."""
    return np.exp(-np.pi*(k**2 + m0**2)/(g*abs(ndot)))


def spectrum_number(g, ndot, m0=0.):
    """Number density of that spectrum, (g|ndot|)^{3/2} exp(-pi m0^2/(g|ndot|))/(8 pi^3)."""
    return (g*abs(ndot))**1.5*np.exp(-np.pi*m0**2/(g*abs(ndot)))/(8*np.pi**3)


class Field:
    """One homogeneous trial: field, radiation and companion modes."""

    def __init__(self, K=1., g=1., m0=0., n_star=0., A=0., kmax=4., J=120, dV=None, prescribed=None,
                 track_tau=False, M=None, delta_n=None):
        self.K, self.g, self.m0, self.n_star, self.A = K, g, m0, n_star, A
        # Stage 2's saturating law, m_C^2 = m0^2 + M^2 tanh^2[(n-n*)/delta_n] with g = M/delta_n. When
        # delta_n is None every expression below is stage 1's, unchanged to the last bit.
        self.M, self.delta_n = M, delta_n
        if delta_n is not None:
            self.M = M if M is not None else g*delta_n
            self.g = self.M/delta_n
        self.dV = dV or (lambda n: 0.)
        self.prescribed = prescribed          # (t) -> (n, ndot), for trials without back-reaction
        self.track_tau = track_tau            # carry the optical time int dt/n; only for histories with n > 0
        self.k, self.w = grid(kmax, J)
        self.J = len(self.k)

    def omega(self, n):
        if self.delta_n is None:
            return np.sqrt(self.k**2 + self.m0**2 + self.g**2*(n - self.n_star)**2)
        return np.sqrt(self.k**2 + self._mass(n)**2)

    def _mass(self, n):
        """The companion mass at index n: quadratic in stage 1, saturating in stage 2."""
        if self.delta_n is None:
            return np.sqrt(self.m0**2 + self.g**2*(n - self.n_star)**2)
        return np.sqrt(self.m0**2 + self.M**2*np.tanh((n - self.n_star)/self.delta_n)**2)

    def dmass2(self, n):
        """d m_C^2 / dn, which is what the back-reaction and the adiabatic coefficient need."""
        if self.delta_n is None:
            return 2*self.g**2*(n - self.n_star)
        u = (n - self.n_star)/self.delta_n
        return 2*self.M**2*np.tanh(u)/np.cosh(u)**2/self.delta_n

    def unpack(self, y):
        J = self.J
        return (y[0], y[1], y[2:2+J],
                y[2+J:2+2*J] + 1j*y[2+2*J:2+3*J], y[2+3*J:2+4*J] + 1j*y[2+4*J:2+5*J])

    def rhs(self, t, y):
        n, nd, th, a, b = self.unpack(y)
        if self.prescribed is not None:
            n, nd = self.prescribed(t)
        x = n - self.n_star
        om = self.omega(n)
        if self.delta_n is None:
            c = self.g**2*x*nd/(2*om**2)          # omega'/(2 omega)
        else:
            c = self.dmass2(n)*nd/(4*om**2)
        e = np.exp(2j*th)
        S = float(np.sum(self.w*(np.abs(b)**2 + np.real(a*np.conj(b)*np.conj(e)))/om))
        force = self.g**2*x*S if self.delta_n is None else self.dmass2(n)*S/2
        ndd = 0. if self.prescribed is not None else (self.A/n**2 - self.dV(n) - force)/self.K
        ad, bd = c*e*b, c*np.conj(e)*a
        tail = [[1/n]] if self.track_tau else []          # the last state, where carried, is the optical time
        return np.concatenate([[nd, ndd], om, ad.real, ad.imag, bd.real, bd.imag] + tail)

    def state(self, t, y, V=None):
        """Energies, occupations and speeds at one time."""
        n, nd, th, a, b = self.unpack(y)
        if self.prescribed is not None:
            n, nd = self.prescribed(t)
        om = self.omega(n)
        nk = np.abs(b)**2
        rho = float(np.sum(self.w*om*nk))
        num = float(np.sum(self.w*nk))
        half_dm2 = self.g**2*(n - self.n_star) if self.delta_n is None else self.dmass2(n)/2
        adiab = float(np.max(np.abs(half_dm2*nd/om**3))) if self.g else 0.
        return dict(t=float(t), n=float(n), ndot=float(nd), rho_C=rho, number=num,
                    kinetic=float(self.K*nd**2/2), radiation=float(self.A/n) if self.A else 0.,
                    potential=float(V(n)) if V else 0.,
                    energy=float(self.K*nd**2/2 + (V(n) if V else 0.) + (self.A/n if self.A else 0.) + rho),
                    v_rms=float(np.sqrt(np.sum(self.w*nk*self.k**2/om**2)/num)) if num > 0 else 0.,
                    mass=float(self._mass(n)),
                    max_nonadiabaticity=adiab, nk=nk)

    def y0(self, n0, v0):
        J = self.J
        tail = [[0.]] if self.track_tau else []
        return np.concatenate([[n0, v0], np.zeros(J), np.ones(J), np.zeros(J), np.zeros(J), np.zeros(J)] + tail)

    def run(self, n0, v0, T, rtol=1e-11, atol=1e-13, n_out=201, events=None, V=None, t0=0.):
        sol = solve_ivp(self.rhs, (t0, T), self.y0(n0, v0), method='DOP853', rtol=rtol, atol=atol,
                        t_eval=np.linspace(t0, T, n_out), events=events, dense_output=True)
        states = [self.state(sol.t[i], sol.y[:, i], V) for i in range(sol.y.shape[1])]
        e0 = states[0]['energy'] or 1.
        return dict(sol=sol, states=states, nfev=int(sol.nfev),
                    energy_error=float(max(abs(s['energy'] - states[0]['energy']) for s in states)/abs(e0)),
                    events=[[self.state(t, y, V) for t, y in zip(ts, ys)]
                            for ts, ys in zip(sol.t_events or [], sol.y_events or [])])

    def kinetic(self, n0, v0, T, nk, rtol=1e-11, atol=1e-13, n_out=201, events=None, V=None, t0=0.):
        """Continuation with frozen occupations: K n'' = A/n^2 - V'(n) - g^2 (n-n*) sum w n_k/omega."""
        def f(t, y):
            n, nd = y
            om = self.omega(n)
            S = float(np.sum(self.w*nk/om))
            force = self.g**2*(n - self.n_star)*S if self.delta_n is None else self.dmass2(n)*S/2
            return [nd, (self.A/n**2 - self.dV(n) - force)/self.K]
        sol = solve_ivp(f, (t0, T), [n0, v0], method='DOP853', rtol=rtol, atol=atol,
                        t_eval=np.linspace(t0, T, n_out), events=events, dense_output=True)
        def st(t, y):
            n, nd = y
            om = self.omega(n)
            rho = float(np.sum(self.w*nk*om))
            num = float(np.sum(self.w*nk))
            return dict(t=float(t), n=float(n), ndot=float(nd), rho_C=rho, number=num,
                        kinetic=float(self.K*nd**2/2), radiation=float(self.A/n) if self.A else 0.,
                        potential=float(V(n)) if V else 0.,
                        energy=float(self.K*nd**2/2 + (V(n) if V else 0.) + (self.A/n if self.A else 0.) + rho),
                        v_rms=float(np.sqrt(np.sum(self.w*nk*self.k**2/om**2)/num)) if num > 0 else 0.,
                        mass=float(self._mass(n)))
        states = [st(sol.t[i], sol.y[:, i]) for i in range(sol.y.shape[1])]
        return dict(sol=sol, states=states,
                    energy_error=float(max(abs(s['energy'] - states[0]['energy']) for s in states)/abs(states[0]['energy'] or 1.)),
                    events=[[st(t, y) for t, y in zip(ts, ys)] for ts, ys in zip(sol.t_events or [], sol.y_events or [])])


def stop_at_rest():
    """Event: the field's speed reaches zero from above (the trapping point)."""
    def ev(t, y):
        return y[1]
    ev.terminal, ev.direction = True, -1
    return ev


def trapping_distance(K, g, ndot, m0=0.):
    """Analytic stopping distance in n: K ndot^2/(2 N g) with N the crossing's number density."""
    return K*ndot**2/(2*spectrum_number(g, ndot, m0)*g)


def ray_arrival(sol, J, t_e, distance, newton=4):
    """When a signal emitted at t_e arrives over a fixed coordinate distance.

    The optical time tau = int dt/n is the last state the field integrates, so it carries the solver's own
    accuracy; the arrival solves tau(t_o) = tau(t_e) + distance by bisection and then Newton, dtau/dt = 1/n.
    """
    i = 2 + 5*J
    tau = lambda t: float(sol.sol(t)[i])
    target = tau(t_e) + distance
    lo, hi = t_e, sol.t[-1]
    if tau(hi) < target:
        raise ValueError('the signal has not arrived within the integrated span')
    for _ in range(80):
        mid = .5*(lo + hi)
        lo, hi = (mid, hi) if tau(mid) < target else (lo, mid)
        if hi - lo <= 1e-12*max(1., abs(hi)):
            break
    t = .5*(lo + hi)
    for _ in range(newton):
        t -= (tau(t) - target)*float(sol.sol(t)[0])
    return t

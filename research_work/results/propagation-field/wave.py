"""Evolving propagation field (PF-1): exact kinematics, a driven finite-difference test and
radiation/field energy exchange. See protocol.md.

Proposed wave law (a hypothesis, not established physics), c = 1:

    d/dt (n dA/dt) - (1/n) Lap A = 0,   homogeneous index n(t) > 0,

with n defined against fixed material clocks and rulers. In tau = int dt/n this is the
ordinary wave equation, so the kinematic statements below are exact for any positive n(t),
not only in the short-wavelength limit.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import hilbert


# Exact kinematics for a linear index n(t) = 1 + ndot t (n = 1 at the present epoch t = 0).

def linear(ndot):
    n = lambda t: 1 + ndot*np.asarray(t, float)
    tau = lambda t: np.log(n(t))/ndot
    t_of_tau = lambda s: (np.exp(ndot*np.asarray(s, float)) - 1)/ndot
    return dict(n=n, tau=tau, t_of_tau=t_of_tau, ndot=ndot)


def arrival(field, t_e, D):
    """Observer time for a signal emitted at t_e over fixed separation D: tau(t_o) = tau(t_e) + D."""
    return field['t_of_tau'](field['tau'](t_e) + D)


def redshift_and_stretch(field, t_e, D):
    """1+z = omega_e/omega_o = n(t_o)/n(t_e); dt_o/dt_e from differentiating the arrival map."""
    t_o = arrival(field, t_e, D)
    return dict(t_o=t_o, one_plus_z=field['n'](t_o)/field['n'](t_e), stretch=field['n'](t_o)/field['n'](t_e))


def luminosity_distance(z, alpha):
    """Fixed Euclidean geometry, photon number conserved, linear index: D = ln(1+z)/alpha, D_L = (1+z) D."""
    z = np.asarray(z, float)
    return (1 + z)*np.log1p(z)/alpha


# Driven finite-difference test in physical time (independent of the tau mapping).

def emitter_signal(t, pulses, nu, sigma):
    """Material-clock oscillator: Gaussian-envelope carrier of frequency nu at each pulse time."""
    t = np.asarray(t, float)
    return sum(np.exp(-(t - tp)**2/(2*sigma**2))*np.sin(2*np.pi*nu*(t - tp)) for tp in pulses)


def simulate(n_of_t, t_start, t_end, x_obs, length, pulses, nu, sigma, points_per_wavelength=64, cfl=.25):
    """Integrate d/dt(n A_t) = (1/n) A_xx on [0, length] with A(0,t) driven by the emitter and a
    one-way outflow boundary at x=length. Staggered leapfrog with pi = n A_t and a fourth-order
    Laplacian. Returns the observer time series and the wave energy history."""
    # Emitted wavelength is 1/(n nu); the spatial wavevector is then preserved, so the shortest
    # wavelength present is the one emitted at the largest index.
    lam = 1/(max(float(n_of_t(tp)) for tp in pulses)*nu)
    dx = lam/points_per_wavelength
    x = np.arange(0, length + dx/2, dx)
    vmax = 1/min(float(n_of_t(t_start)), float(n_of_t(t_end)))
    dt = cfl*dx/vmax
    steps = int(np.ceil((t_end - t_start)/dt))
    A = np.zeros(len(x))
    pi = np.zeros(len(x))
    io = int(round(x_obs/dx))
    times, obs, energy = [], [], []
    t = t_start
    for s in range(steps):
        # pi at t + dt/2 from the Laplacian at t
        lap = np.zeros_like(A)
        lap[2:-2] = (-A[4:] + 16*A[3:-1] - 30*A[2:-2] + 16*A[1:-3] - A[:-4])/(12*dx*dx)
        lap[1] = (A[2] - 2*A[1] + A[0])/(dx*dx)
        lap[-2] = (A[-1] - 2*A[-2] + A[-3])/(dx*dx)
        pi += dt*lap/n_of_t(t) if s else .5*dt*lap/n_of_t(t)
        n_half = n_of_t(t + dt/2)
        A_new = A + dt*pi/n_half
        t += dt
        A_new[0] = emitter_signal(t, pulses, nu, sigma)
        # first-order outflow: A_t + (1/n) A_x = 0 at the right edge
        A_new[-1] = A[-1] - dt/n_half*(A[-1] - A[-2])/dx
        A = A_new
        times.append(t)
        obs.append(A[io])
        if s % 50 == 0:
            ax = np.gradient(A, dx)
            energy.append((t, float(np.sum(pi*pi/(2*n_half) + ax*ax/(2*n_half))*dx)))
    return dict(t=np.array(times), signal=np.array(obs), energy=np.array(energy), dx=dx, dt=dt)


def pulse_measurements(t, signal, windows):
    """Centroid time, envelope width and carrier frequency of each pulse (analytic signal)."""
    env = np.abs(hilbert(signal))
    phase = np.unwrap(np.angle(hilbert(signal)))
    freq = np.gradient(phase, t)/(2*np.pi)
    out = []
    for lo, hi in windows:
        m = (t >= lo) & (t <= hi)
        w = env[m]**2
        tc = float(np.sum(w*t[m])/np.sum(w))
        width = float(np.sqrt(np.sum(w*(t[m] - tc)**2)/np.sum(w)))
        core = m & (np.abs(t - tc) < width)
        out.append(dict(centroid=tc, width=width, carrier=float(np.sum(env[core]**2*freq[core])/np.sum(env[core]**2))))
    return out


# Coupled radiation + dynamical index, exact in Fourier modes of a periodic box.

def coupled(k, A0, P0, n0, ndot0, M, V, dV, t_end, rtol=1e-12, points=400):
    """Modes A_k, pi_k with A_k' = pi_k/n, pi_k' = -(k^2/n) A_k; index with inertia M and potential V:
    M n'' = -V'(n) + E/n, E = sum (pi_k^2 + k^2 A_k^2)/(2n). Total H = M n'^2/2 + V(n) + E is conserved."""
    k = np.asarray(k, float)
    m = len(k)

    def rhs(t, y):
        n, p = y[0], y[1]
        A, P = y[2:2 + m], y[2 + m:]
        E = np.sum(P*P + k*k*A*A)/(2*n)
        return np.concatenate([[p/M, -dV(n) + E/n], P/n, -(k*k/n)*A])

    y0 = np.concatenate([[n0, M*ndot0], A0, P0])
    sol = solve_ivp(rhs, (0, t_end), y0, method='DOP853', rtol=rtol, atol=1e-14, dense_output=True)
    ts = np.linspace(0, t_end, points)
    Y = sol.sol(ts)
    n, p = Y[0], Y[1]
    A, P = Y[2:2 + m], Y[2 + m:]
    Ek = (P*P + (k*k)[:, None]*A*A)/(2*n)
    E = Ek.sum(0)
    return dict(t=ts, n=n, ndot=p/M, wave_energy=E, field_energy=p*p/(2*M) + V(n), mode_energy=Ek,
                photon_number=(Ek*n/k[:, None]).sum(0), total=p*p/(2*M) + V(n) + E, success=sol.success)

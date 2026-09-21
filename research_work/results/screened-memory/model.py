"""SM-1: frozen screened square-root total-potential model (not a complete theory).

The dimensionless spherical law is z=h(x)=x+sqrt(x)/(1+x*x), x=g_N/a0.
mu(z)=x/z defines an AQUAL total-potential completion, NOT an independently
sourced additive scalar. All public functions accept nonnegative finite arrays.
"""
from __future__ import annotations
import numpy as np
from numpy.polynomial.legendre import leggauss

A0_SI = 6.54e-11  # frozen rounded PM-1 fitted scale; not newly fitted or derived
KPC_M = 3.0856775814913673e19
DERIVATIVE_LOWER_BOUND = 1 - (75 / 128) * (3 / 5) ** 0.75
_GLX, _GLW = leggauss(64)
_GLX, _GLW = (_GLX + 1) / 2, _GLW / 2


def positive_array(value):
    a = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError('Expected finite nonnegative values')
    return a


def excess(x):
    x = positive_array(x)
    return np.sqrt(x) / (1 + x * x)


def h(x):
    x = positive_array(x)
    return x + excess(x)


def hprime(x):
    x = positive_array(x)
    with np.errstate(divide='ignore'):
        return 1 + (1 - 3 * x * x) / (2 * np.sqrt(x) * (1 + x * x) ** 2)


def inverse(z):
    """Invert h through u=sqrt(x), using a bracketed Newton iteration.

    A converged entry stays fixed. Keeping it in the bracket iteration would
    destroy convergence when roundoff makes the residual exactly zero.
    """
    z = positive_array(z)
    lo = 2 * z / (1 + np.sqrt(1 + 4 * z))
    hi = np.sqrt(z)
    u = np.where(z < 1, lo, hi)
    done = z == 0
    for _ in range(60):
        u4 = u ** 4
        f = u * u + u / (1 + u4) - z
        done |= np.abs(f) <= 4e-15 * np.maximum(z, 1e-300)
        if np.all(done):
            break
        lo = np.where((f < 0) & ~done, u, lo)
        hi = np.where((f > 0) & ~done, u, hi)
        derivative = 2 * u + (1 - 3 * u4) / (1 + u4) ** 2
        candidate = u - f / derivative
        inside = (candidate > lo) & (candidate < hi) & np.isfinite(candidate)
        candidate = np.where(inside, candidate, (lo + hi) / 2)
        u = np.where(done, u, candidate)
    if not np.all(done):
        raise ArithmeticError('Screened constitutive inverse did not converge')
    return u * u


def mu(z):
    z = positive_array(z)
    x = inverse(z)
    return np.divide(x, z, out=np.zeros_like(z), where=z > 0)


def weight(z):
    """Log-conductance relaxation weight mu/(mu+z*mu')."""
    z = positive_array(z)
    x = inverse(z)
    with np.errstate(invalid='ignore'):
        result = mu(z) * hprime(x)
    return np.where(z == 0, 0.5, result)


def W(z):
    """Dimensionless gradient energy int_0^z mu(v)*v dv.

    Parameter u=sqrt(x) removes the singular endpoint. The integrand is
    2*u^3 + u^2*(1-3*u^4)/(1+u^4)^2 = x*dh/du > 0.
    Physical gradient energy is a0^2*W(|grad U|/a0)/(4*pi*G).
    """
    z = positive_array(z)
    top = np.sqrt(inverse(z))
    u = top[..., None] * _GLX
    u4 = u ** 4
    f = 2 * u ** 3 + u * u * (1 - 3 * u4) / (1 + u4) ** 2
    return top * (f @ _GLW)


def equation_for_repository():
    """Reuse PM-2A's solver without editing its archived equations or data."""
    import fields
    return fields.Equation('screened_root_SM1', mu, h, weight,
                           lambda u: 2 * W(np.sqrt(positive_array(u))))


def particle_weights(nodes, q, width=0.4):
    """Compact C2 particle source and its exact positional derivative.

    Normalization is differentiated too; source and reaction use the same
    finite-dimensional interaction energy, including near-boundary effects.
    """
    d = nodes - q
    t = np.maximum(1 - (d / width) ** 2, 0)
    k = t ** 3
    kp = 6 * d / width ** 2 * t ** 2
    total = k.sum()
    if total <= 0:
        raise ValueError('Particle left the computational domain')
    return k / total, kp / total - k * kp.sum() / total ** 2


def toy_run(dt=0.002, duration=2.0, gamma=0.0, n=161):
    """1D finite-domain reciprocal matter/field test, NOT a galaxy simulation.

    Units 4*pi*G=a0=1. Boundaries U=0. c_star=.7 is a declared numerical
    choice, not a physical prediction. The primary gamma=0 system is
    conservative. gamma>0 adds an explicit nonnegative heat ledger, but not
    a microscopic bath or a closed momentum ledger for that bath.
    """
    if dt <= 0 or gamma < 0 or n < 20:
        raise ValueError('Invalid integration parameters')
    nodes = np.linspace(-4, 4, n)
    dx = nodes[1] - nodes[0]
    cstar = 0.7
    field_mass = dx / cstar ** 2
    if dt / dx * cstar / np.sqrt(DERIVATIVE_LOWER_BOUND) > 0.3:
        raise ValueError('Time step violates the declared conservative CFL bound')
    U, P = np.zeros(n), np.zeros(n)
    q, p, mass = np.array([-1., 1.]), np.zeros(2), np.full(2, .2)
    heat = 0.
    steps = int(round(duration / dt))
    if abs(steps * dt - duration) > 1e-12:
        raise ValueError('Duration must be an integer multiple of dt')

    def forces():
        gradient = np.diff(U) / dx
        flux = mu(np.abs(gradient)) * gradient
        field = np.zeros(n)
        field[1:-1] = np.diff(flux)
        body = np.zeros(2)
        for i in range(2):
            w, dw = particle_weights(nodes, q[i])
            field -= mass[i] * w
            body[i] = -mass[i] * (dw @ U)
        field[[0, -1]] = 0
        return field, body

    def channels():
        kinetic = float(np.sum(p * p / (2 * mass)))
        fkin = float(np.sum(P * P) / (2 * field_mass))
        gradient = float(dx * np.sum(W(np.abs(np.diff(U) / dx))))
        interaction = float(sum(mass[i] * (particle_weights(nodes, q[i])[0] @ U)
                                for i in range(2)))
        return np.array([kinetic, fkin, gradient, interaction, heat])

    records = [channels()]
    ff, fp = forces()
    for step in range(steps):
        if gamma:
            old = float(P @ P / (2 * field_mass))
            P *= np.exp(-gamma * dt / 2)
            heat += old - float(P @ P / (2 * field_mass))
        P += dt / 2 * ff
        p += dt / 2 * fp
        U += dt * P / field_mass
        q += dt * p / mass
        ff, fp = forces()
        P += dt / 2 * ff
        p += dt / 2 * fp
        if gamma:
            old = float(P @ P / (2 * field_mass))
            P *= np.exp(-gamma * dt / 2)
            heat += old - float(P @ P / (2 * field_mass))
        if step % max(1, steps // 50) == 0 or step == steps - 1:
            records.append(channels())
    history = np.asarray(records)
    energy = history.sum(1)
    scale = max(float(np.max(np.abs(history).sum(1))), 1e-300)
    return dict(dt=dt, n=n, gamma=gamma, positions=q.tolist(), velocities=(p / mass).tolist(),
                normalized_energy_error=float(np.max(np.abs(energy - energy[0])) / scale),
                final_channels=dict(zip(('matter_kinetic', 'field_kinetic', 'field_gradient',
                                         'interaction', 'heat'), history[-1].tolist())),
                heat_nondecreasing=bool(np.all(np.diff(history[:, -1]) >= -1e-14)),
                source='two initially resting extended particles; initially empty field',
                scope='one-dimensional discretized energy/reaction control; not physical gravity validation')

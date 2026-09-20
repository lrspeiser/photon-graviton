"""RUT-1 stage 12: the angular-momentum budget of the reciprocal memory field.

Stage 4R gave the reciprocal form C = W*h, with W*W = K, an exact energy ledger:

    H = sum_i m_i [|v_i|^2/2 + Phi_ext(X_i) - C(X_i)] + tau_form/(2 alpha) int h_t^2 + 1/(2 alpha tau_keep) int h^2
    dH/dt = -(gamma/alpha) int h_t^2 <= 0

and said in as many words that it "does not finish the momentum accounting". This module finishes it.

The external potential -GM/r is axisymmetric and exerts no torque, so the whole matter torque comes from the
field. With the acceleration +grad C and rho = sum_i m_i delta(x - X_i),

    dL_m/dt = int rho (x x grad C)_z = int rho d_theta (W*h).

W is radially symmetric, so d_theta commutes with W*, and convolution by a symmetric kernel is self-adjoint:

    int rho d_theta (W*h) = int (W*rho) d_theta h = (1/alpha) int (tau_f h_tt + gamma h_t + h/tau_keep) d_theta h.

The last term is (1/2) int d_theta(h^2) = 0, and int h_tt d_theta h = d/dt int h_t d_theta h, because
int h_t d_theta h_t vanishes for the same reason. So, writing

    J = int h_t d_theta h,     L_field = -(tau_form/alpha) J,     torque = (gamma/tau_form) L_field,

the budget closes exactly, in the same shape as the energy one:

    d(L_m + L_field)/dt = -(gamma/tau_form) L_field.

Unlike the dissipation, the torque is NOT sign-definite: the field gives angular momentum back as readily as it
takes it, which is what a memory that decays should do.

Two things make this approximate in the implementation, and both are measured rather than assumed. The box is
periodic, so d_theta is only meaningful where the field has died away before the edge; and the square lattice
of modes is not rotation invariant, so the truncation itself can exert a small torque. `L_field` is a physical
quantity and must therefore not depend on the box or the mode count -- a gate, not a hope.

The moment integrals are taken on the field's own grid. `ifft2` returns the field at x_j = j L / n, on [0, L),
so the signed coordinate a moment needs is that wrapped into [-L/2, L/2). Getting this wrong by half a box
leaves the energy ledger untouched -- energy is computed by Parseval and never sees a position -- while making
L_field vary threefold with the box and destroying the closure. It is a declared negative control below.
"""
import numpy as np


def coordinates(field):
    """The signed position grid of `ifft2`'s output: x_j = j L/n wrapped into [-L/2, L/2)."""
    n, L = field.n, field.L
    ax = (np.arange(n)*L/n + L/2) % L - L/2
    X, Y = np.meshgrid(ax, ax, indexing='ij')
    return X, Y, (L/n)**2


def real_space(field, coeff):
    """h(x) = sum_k coeff_k exp(+i k.x) on that grid, matching `_readout`'s sign convention."""
    return np.real(np.fft.ifft2(coeff)*field.n*field.n)


def d_theta(field, coeff, offset=0.):
    """x d_y h - y d_x h, the derivatives spectral. `offset` displaces the coordinate origin and exists only
    so the control can be run: offset = 0.5 is the half-box mistake."""
    X, Y, _ = coordinates(field)
    if offset:
        L = field.L
        X, Y = X + offset*L, Y + offset*L
    hx = real_space(field, 1j*field.kx*coeff)
    hy = real_space(field, 1j*field.ky*coeff)
    return X*hy - Y*hx


def J(field, offset=0.):
    """int h_t d_theta h over the box."""
    _, _, cell = coordinates(field)
    return float(np.sum(real_space(field, field.hd)*d_theta(field, field.h, offset))*cell)


def field_angular_momentum(field, offset=0.):
    return -(field.tau_form/field.alpha)*J(field, offset)


def torque(field, offset=0.):
    """(gamma/tau_form) L_field: what leaves the matter-plus-field total, of either sign."""
    return (field.gamma/field.tau_form)*field_angular_momentum(field, offset)


def matter_angular_momentum(x, v, m):
    return float(np.sum(m*(x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])))


class Ledger:
    """Both budgets carried along a run by the trapezoid rule, in the same shape:

        energy:  H(t) - H(0) + int dissipation dt        must vanish
        angular: (L_m + L_field)(t) - (...)(0) + int torque dt   must vanish
    """

    def __init__(self, field, x, v, m, H_of, offset=0.):
        self.field, self.m, self.H_of, self.offset = field, m, H_of, offset
        self.H0 = H_of(x, v)
        self.L0 = matter_angular_momentum(x, v, m) + field_angular_momentum(field, offset)
        self.D = self.T = 0.
        self._d = field.dissipation_rate()
        self._t = torque(field, offset)

    def step(self, h):
        d, t = self.field.dissipation_rate(), torque(self.field, self.offset)
        self.D += .5*h*(self._d + d)
        self.T += .5*h*(self._t + t)
        self._d, self._t = d, t

    def read(self, x, v):
        Lm = matter_angular_momentum(x, v, self.m)
        Lf = field_angular_momentum(self.field, self.offset)
        return dict(energy_balance=self.H_of(x, v) - self.H0 + self.D,
                    angular_balance=Lm + Lf - self.L0 + self.T,
                    L_matter=Lm, L_field=Lf, L_total_0=self.L0, H_0=self.H0,
                    dissipated=self.D, torque_integral=self.T)

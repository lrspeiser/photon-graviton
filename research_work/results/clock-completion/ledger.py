"""CC-1 energy ledger (N5) and the history the declared field dynamics imply (O7, O8). Units c = 1.

Per comoving volume, the homogeneous toy has the Lagrangian
    L = (M/2) ndot^2 - V(n)                                     clock field (V = 0 in the toy)
      + sum_i [(n/2) qdot_i^2 - k_i^2 q_i^2/(2n)]               radiation modes (PF-1's wave law)
      - sum_a m_a sqrt(1 - n^2 xdot_a^2)                        free particles in g_m
      + (mu/2) n^2 |xdot|^2 + kc/(n |x|)                        one bound pair (nonrelativistic Coulomb)
Every matter term is that matter's action in g_m, so the field receives exactly the work matter does
on it. The Hamiltonian is
    H = p_n^2/(2M) + V + sum (pi^2 + k^2 q^2)/(2n) + sum sqrt(m^2 + p^2/n^2) + |P|^2/(2 mu n^2) - kc/(n|x|).
"""
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RTOL, ATOL = 1e-12, 1e-15

TOY = dict(M=1e7, ndot0=1e-3,
           modes=dict(k=[1., 2., 3.], W=[1., .6, .4]),
           particles=dict(m=[1., 1.], p=[3., .01]),
           pair=dict(mu=1., kc=1e-7, r=1e-3))


class Toy:
    def __init__(self, cfg=TOY):
        self.cfg = cfg
        self.M = cfg['M']
        self.k, self.W0 = np.array(cfg['modes']['k']), np.array(cfg['modes']['W'])
        self.m, self.p = np.array(cfg['particles']['m']), np.array(cfg['particles']['p'])
        self.mu, self.kc, self.r0 = cfg['pair']['mu'], cfg['pair']['kc'], cfg['pair']['r']
        self.omega_pair = math.sqrt(self.kc/(self.mu*self.r0**3))

    def initial(self):
        nd = self.cfg['ndot0']
        v = math.sqrt(self.kc/(self.mu*self.r0))
        # A circular orbit in proper terms at n = 1: the comoving velocity (-ndot r, v) cancels the Hubble flow.
        return np.concatenate([[1., self.M*nd], np.sqrt(2*self.W0)/self.k, np.zeros(3),
                               [self.r0, 0.], [-self.mu*nd*self.r0, self.mu*v]])

    def rhs(self, t, y):
        n, pn, q, pi, x, P = y[0], y[1], y[2:5], y[5:8], y[8:10], y[10:12]
        r = math.hypot(x[0], x[1])
        Ep = np.sqrt(self.m**2 + self.p**2/n**2)
        force = (np.sum(pi*pi + self.k**2*q*q)/(2*n*n) + np.sum(self.p**2/(n**3*Ep))
                 + (P @ P)/(self.mu*n**3) - self.kc/(n*n*r))          # -dH/dn
        return np.concatenate([[pn/self.M, force], pi/n, -self.k**2*q/n, P/(self.mu*n*n), -self.kc*x/(n*r**3)])

    def energies(self, Y):
        """Energy terms for states stored as columns."""
        n, pn, q, pi, x, P = Y[0], Y[1], Y[2:5], Y[5:8], Y[8:10], Y[10:12]
        W = (pi*pi + (self.k**2)[:, None]*q*q)/2
        free = np.sqrt((self.m**2)[:, None] + (self.p**2)[:, None]/n**2)
        pair = (P*P).sum(0)/(2*self.mu*n*n) - self.kc/(n*np.hypot(x[0], x[1]))
        field = pn*pn/(2*self.M)
        rad = W.sum(0)/n
        return dict(field=field, radiation=rad, free=free, pair=pair, W=W, total=field + rad + free.sum(0) + pair)


def run(toy, backward=False, span=1e5, samples_per_orbit=32):
    """Forward until n = 2, or backward until the field turns around (p_n = 0). States are sampled on
    whole bound-pair orbits; the pair's proper period, 2 pi/omega, is constant in t."""
    dt = 2*math.pi/toy.omega_pair/samples_per_orbit
    ev = (lambda t, y: y[1]) if backward else (lambda t, y: y[0] - 2.)
    ev.terminal = True
    sgn = -1 if backward else 1
    sol = solve_ivp(toy.rhs, (0, sgn*span), toy.initial(), method='DOP853', rtol=RTOL, atol=ATOL,
                    t_eval=sgn*np.arange(0, span, dt), events=ev)
    if sol.status != 1 or not np.all(np.isfinite(sol.y)):
        raise FloatingPointError(f'ledger run did not reach its event: {sol.message}')
    return sol


def ledger(toy=None):
    toy = toy or Toy()
    fwd, bwd = run(toy), run(toy, backward=True)
    e = toy.energies(fwd.y)
    e_end = toy.energies(fwd.y_events[0].T)
    H0 = e['total'][0]
    k = len(e['pair'])//32*32
    pair_orbit = e['pair'][:k].reshape(-1, 32).mean(1)
    delta = {name: float(np.sum(e_end[name][..., 0] - e[name][..., 0])) for name in ('field', 'radiation', 'pair')}
    delta['free_relativistic'], delta['free_nonrelativistic'] = (float(e_end['free'][i, 0] - e['free'][i, 0]) for i in (0, 1))
    closure = (delta['field'] + delta['radiation'] + delta['free_relativistic'] + delta['free_nonrelativistic']
               + delta['pair'])/H0
    # turnaround predicted from energy conservation with V = 0: sum W/n + sum sqrt(m^2 + p^2/n^2) + E_pair = H
    f = lambda n: toy.W0.sum()/n + np.sum(np.sqrt(toy.m**2 + toy.p**2/n**2)) + pair_orbit[0] - H0
    n_min_analytic = brentq(f, 1e-6, 1., xtol=1e-15, rtol=1e-15)
    n_min = float(bwd.y_events[0][0, 0])
    checks = dict(
        energy_drift=float(np.max(np.abs(e['total']/H0 - 1))),
        ledger_closure=float(abs(closure)),
        mode_invariant_W_drift=float(np.max(np.abs(e['W']/e['W'][:, :1] - 1))),
        pair_orbit_average_drift=float(np.max(np.abs(pair_orbit/pair_orbit[0] - 1))),
        turnaround_error=float(abs(n_min/n_min_analytic - 1)))
    passed = (checks['energy_drift'] < 1e-8 and checks['ledger_closure'] < 1e-8
              and checks['pair_orbit_average_drift'] < 1e-6 and checks['turnaround_error'] < 1e-6)
    return dict(toy=toy.cfg, H0=float(H0), t_forward=float(fwd.t_events[0][0]), deltas_n_1_to_2=delta,
                field_gain_over_radiation_loss=float(delta['field']/-delta['radiation']),
                turnaround=dict(n_min=n_min, n_min_analytic=float(n_min_analytic), t=float(bwd.t_events[0][0]),
                                one_plus_z_max=float(1/n_min)),
                checks=checks, passed=bool(passed))


def history(T0=2.7255, H_kms_Mpc=70.48):
    """O8 with V = 0 and radiation dominating the turnaround. The field's energy today equals the radiation
    energy it absorbed since the turnaround, so z_max = E_field/E_radiation and q0 = -1/(2 z_max) (coasting).
    The present redshift rate ndot/n is the repository's supernova value, used only as a unit here, not as a
    history. Whether the field's energy gravitates is open; the last column assumes it does."""
    a_rad, G, c, Mpc = 7.565723e-16, 6.67430e-11, 299792458., 3.0856775814913673e22
    u_gamma = a_rad*T0**4
    H = H_kms_Mpc*1e3/Mpc
    u_unit = 3*H*H*c*c/(8*math.pi*G)
    rows = {}
    for label, T in (('recombination, 3000 K', 3000.), ('nucleosynthesis, 8e8 K', 8e8)):
        z = T/T0 - 1
        rows[label] = dict(one_plus_z_max=T/T0, field_over_radiation_energy=z, q0=-1/(2*z),
                           field_energy_density_J_per_m3=z*u_gamma, over_3H2c2_over_8piG=z*u_gamma/u_unit)
    return dict(T0_K=T0, u_gamma_J_per_m3=u_gamma, present_rate_km_s_Mpc=H_kms_Mpc,
                unit_3H2c2_over_8piG_J_per_m3=u_unit, minimum_turnaround=rows)

"""RUT-1 stage 4: the instrumented long-run integrator (protocol-rut1.md, stage 4).

The numerical core is frozen: this module reuses formation.MemoryField and formation.ring_start unchanged
and the same velocity-Verlet step with its adaptive timestep, so what it adds is measurement only. The
question it exists to answer is whether the radial motion left in the two-stage runs is initial adjustment,
coherent oscillation, or a genuine slow instability -- which a last-window average cannot tell apart.

Measured during the integration rather than reconstructed from sparse snapshots:

* the orbital-energy consistency residual. A body's instantaneous specific orbital energy
  eps_i = |v|^2/2 - GM/r - C(X_i, t) is NOT conserved in an evolving field; along the declared equations
  d(eps_i)/dt = -dC/dt at X_i exactly. So R_i(t) = eps_i(t) - eps_i(0) + int_0^t dC/dt(X_i(s), s) ds must
  vanish, and it checks the moving body, the sampled force, the sampled potential and the field update
  together. It is not the matter-field energy budget, and a positive eps is not a permanent-escape test.
* mode-resolved memory power and torque on every body, by angular band, from rings of the force-producing
  field around each body's own radius: specific torque d(C_b)/dtheta and power v . grad(C_b). The bands
  must sum to the total, which is checked against the interpolated force.
* radial mean flow and dispersion in radial bins, so a coherent contraction or breathing is not counted as
  heating the way a bare radial-velocity RMS would count it.
* angular-mode band powers of the writing pattern S, the excitation E and the matured field C on fixed
  rings, since fine structure at the writer spacing is exactly what maturation is meant to suppress.
"""
import numpy as np
from scipy.ndimage import map_coordinates

import formation as FM

GM = FM.GM
BANDS = ((0, 0), (1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 32), (33, 64))


def band_label(b):
    return f'm{b[0]}' if b[0] == b[1] else f'm{b[0]}-{b[1]}'


def _sample(arr, field, pts):
    """Cubic-spline value of one grid array at points, the same interpolation the force uses."""
    c = ((np.atleast_2d(pts) - field.x0)/field.d).T
    return map_coordinates(arr, c, order=3, mode='nearest')


def _writing_pattern(pts, positions, rates, w):
    """S at arbitrary points, analytically: sum_j q_j exp(-|x - X_j|^2 / 2w^2)."""
    d2 = ((pts[:, None, :] - positions[None, :, :])**2).sum(-1)
    return (np.exp(-d2/(2*w*w))*rates[None, :]).sum(1)


def _dC_dt(field, x, rates, positions):
    """dC/dt at the bodies, from the model's own equation evaluated on the current state."""
    C = _sample(field.C[0], field, x)
    if field.E is None:
        return _writing_pattern(x, positions, rates, field.w) - C/field.tau_keep
    return _sample(field.E[0], field, x) - C/field.tau_keep


def _ring_modes(values, M):
    """Angular Fourier coefficients c_m = (1/M) sum_j f(phi_j) exp(-i m phi_j), m = 0..M/2."""
    return np.fft.rfft(values, axis=-1)/M


def _band_power_on_rings(field_values, M):
    """Sum over each band of |c_m|^2 (m > 0 counted twice, for +-m)."""
    c = _ring_modes(field_values, M)
    p = np.abs(c)**2
    out = {}
    for b in BANDS:
        lo, hi = b
        hi = min(hi, c.shape[-1] - 1)
        if lo > hi:
            out[band_label(b)] = 0.
            continue
        s = p[..., lo:hi + 1].sum(-1)
        out[band_label(b)] = s if lo == 0 and hi == 0 else 2*s
    return out


def _mode_power_and_torque(field, x, v, M=128, delta_frac=.5):
    """Per-body, per-band memory power v . grad(C_b) and specific torque dC_b/dtheta.

    For each body, C is sampled on five rings at its own radius r and r +- delta, r +- 2 delta,
    Fourier-decomposed in angle, and each band is reassembled at the body's angle. The radial derivative
    is the fourth-order central difference across the rings -- a second-order one left a 1% discrepancy
    against the interpolated force, (delta/w)^2 as expected -- and the angular derivative is analytic in m.
    """
    n = len(x)
    r = np.linalg.norm(x, axis=1)
    th = np.arctan2(x[:, 1], x[:, 0])
    delta = delta_frac*field.d
    phi = 2*np.pi*np.arange(M)/M
    radii = np.stack([r - 2*delta, r - delta, r, r + delta, r + 2*delta], axis=1)      # (n, 5)
    pts = np.stack([radii[..., None]*np.cos(phi), radii[..., None]*np.sin(phi)], axis=-1)
    vals = _sample(field.C[0], field, pts.reshape(-1, 2)).reshape(n, 5, M)
    c = _ring_modes(vals, M)                                               # (n, 5, M/2+1)
    m = np.arange(c.shape[-1])
    phase = np.exp(1j*m[None, :]*th[:, None])                              # (n, M/2+1)
    vr = np.sum(v*x, axis=1)/r
    vt = (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])/r
    power, torque = {}, {}
    for b in BANDS:
        lo, hi = b
        hi = min(hi, c.shape[-1] - 1)
        sel = slice(lo, hi + 1)
        weight = np.where(m[sel] == 0, 1., 2.)
        val = np.real((c[:, :, sel]*phase[:, None, sel]*weight).sum(-1))              # (n, 5)
        dth = np.real((1j*m[sel]*c[:, 2, sel]*phase[:, sel]*weight).sum(-1))           # (n,)
        dr = (val[:, 0] - 8*val[:, 1] + 8*val[:, 3] - val[:, 4])/(12*delta)
        power[band_label(b)] = vr*dr + vt/r*dth
        torque[band_label(b)] = dth
    return power, torque


def _radial_bins(x, v, edges):
    """Mean radial flow and radial/azimuthal dispersion about per-bin means, equal-mass bodies."""
    r = np.linalg.norm(x, axis=1)
    vr = np.sum(v*x, axis=1)/r
    vt = (x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0])/r
    idx = np.digitize(r, edges)
    res_r = np.zeros_like(r)
    res_t = np.zeros_like(r)
    counted = np.zeros(len(r), bool)
    for k in np.unique(idx):
        sel = idx == k
        if sel.sum() >= 2:
            res_r[sel] = vr[sel] - vr[sel].mean()
            res_t[sel] = vt[sel] - vt[sel].mean()
            counted[sel] = True
    return dict(mean_radial_flow=float(vr.mean()),
                radial_dispersion=float(np.sqrt(np.mean(res_r[counted]**2))) if counted.any() else 0.,
                azimuthal_dispersion=float(np.sqrt(np.mean(res_t[counted]**2))) if counted.any() else 0.,
                bodies_in_shared_bins=int(counted.sum()),
                radial_velocity_rms=float(np.sqrt(np.mean(vr**2))))


def run_instrumented(positions, velocities, rates, w, tau_keep, t_end, h_max, tau_form=0., eta=.01,
                     r_min=.25, half_width=2.5, spacing=None, memory=True, prime_ring_R=None,
                     period=2*np.pi, record_every=.05, band_record_every=.5, diag_every=.0025,
                     checkpoints=(5, 10, 20, 50, 100), ring_radii=(.8, .9, 1., 1.1), ring_M=512, body_ring_M=128,
                     band_diagnostics=True, rotate=0.):
    """The formation step of formation.run, unchanged, with the stage 4 instrumentation around it."""
    spacing = spacing if spacing is not None else w/5
    field = FM.MemoryField(half_width, spacing, w, tau_keep, tau_form)
    if prime_ring_R is not None:
        field.prime_with_ring(prime_ring_R, A=float(np.sum(rates))*tau_keep)
    x = np.array(positions, float)
    v = np.array(velocities, float)
    if rotate:                                   # the grid-rotation control: same bodies, turned on the grid
        c_, s_ = np.cos(rotate), np.sin(rotate)
        rot = np.array([[c_, -s_], [s_, c_]])
        x, v = x@rot.T, v@rot.T
    rates = np.asarray(rates, float)
    n = len(x)
    edges = np.arange(.2, 2.5001, .1)

    def accel(xx):
        r = np.linalg.norm(xx, axis=1, keepdims=True)
        a = -GM*xx/np.maximum(r, 1e-12)**3
        return a + field.sample_gradient(xx) if memory else a

    def eps(xx, vv):
        base = .5*np.sum(vv*vv, axis=1) - GM/np.linalg.norm(xx, axis=1)
        return base - _sample(field.C[0], field, xx) if memory else base

    a = accel(x)
    eps0 = eps(x, v)
    kin0 = .5*np.sum(v*v, axis=1) - GM/np.linalg.norm(x, axis=1)
    dCdt = _dC_dt(field, x, rates, x) if memory else np.zeros(n)
    g = field.sample_gradient(x) if memory else np.zeros_like(x)
    P = np.sum(v*g, axis=1)
    Q = x[:, 0]*g[:, 1] - x[:, 1]*g[:, 0]      # specific torque (x cross grad C)_z = dC/dtheta
    I_dC = np.zeros(n)                 # int dC/dt at the body, trapezoid every step
    W_tot = np.zeros(n)                # int v . grad C, every step
    J_tot = np.zeros(n)                # int (x x grad C)_z, every step
    W_band = {band_label(b): np.zeros(n) for b in BANDS}
    J_band = {band_label(b): np.zeros(n) for b in BANDS}
    W_tot_cadence = np.zeros(n)
    band_check = []
    last_diag = None
    L0 = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]

    rec, cps = [], {}
    t, status = 0., FM.COMPLETED
    next_rec, next_diag = 0., 0.
    band_state = {'next': 0.}
    cp_times = [c*period for c in checkpoints if c*period <= t_end + 1e-9]

    def diagnostics(tt):
        if not memory or not band_diagnostics:
            return None
        pw, tq = _mode_power_and_torque(field, x, v, M=body_ring_M)
        gg = field.sample_gradient(x)
        tot = np.sum(v*gg, axis=1)
        s = sum(pw.values())
        scale = max(float(np.max(np.abs(tot))), 1e-300)
        band_check.append(float(np.max(np.abs(s - tot))/scale))
        return tt, pw, tq, tot

    def record(tt):
        r = np.linalg.norm(x, axis=1)
        L = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
        gg = field.sample_gradient(x) if memory else np.zeros_like(x)
        support = -np.sum(gg*x/r[:, None], axis=1)/(GM/r**2)
        e_now = eps(x, v)
        kin = .5*np.sum(v*v, axis=1) - GM/r
        row = dict(t_periods=float(tt/period), mean_radius=float(r.mean()), min_radius=float(r.min()),
                   max_radius=float(r.max()), radius_spread=float(r.max() - r.min()),
                   L_mean=float(np.mean(L/L0)), L_min=float(np.min(L/L0)), L_max=float(np.max(L/L0)),
                   support_mean=float(support.mean()), support_min=float(support.min()),
                   support_max=float(support.max()),
                   eps_residual_max=float(np.max(np.abs(e_now - eps0 + I_dC))),
                   eps_residual_rel=float(np.max(np.abs(e_now - eps0 + I_dC)/np.abs(eps0))),
                   work_residual_max=float(np.max(np.abs(kin - kin0 - W_tot))),
                   torque_residual_max=float(np.max(np.abs((L - L0) - J_tot))),
                   specific_orbital_energy_max=float(e_now.max()),
                   **_radial_bins(x, v, edges))
        if memory and tt >= band_state['next'] - 1e-12:
            band_state['next'] = tt + band_record_every*period
            ring_pts = np.stack([np.outer(ring_radii, np.cos(2*np.pi*np.arange(ring_M)/ring_M)),
                                 np.outer(ring_radii, np.sin(2*np.pi*np.arange(ring_M)/ring_M))], axis=-1)
            flat = ring_pts.reshape(-1, 2)
            Cv = _sample(field.C[0], field, flat).reshape(len(ring_radii), ring_M)
            Sv = _writing_pattern(flat, x, rates, field.w).reshape(len(ring_radii), ring_M)
            row['C_band_power'] = {k: np.asarray(vv).tolist() for k, vv in _band_power_on_rings(Cv, ring_M).items()}
            row['S_band_power'] = {k: np.asarray(vv).tolist() for k, vv in _band_power_on_rings(Sv, ring_M).items()}
            if field.E is not None:
                Ev = _sample(field.E[0], field, flat).reshape(len(ring_radii), ring_M)
                row['E_band_power'] = {k: np.asarray(vv).tolist() for k, vv in _band_power_on_rings(Ev, ring_M).items()}
            row['W_band_sum_over_bodies'] = {k: float(vv.sum()) for k, vv in W_band.items()}
            row['J_band_sum_over_bodies'] = {k: float(vv.sum()) for k, vv in J_band.items()}
            row['W_total_sum'] = float(W_tot.sum())
            row['J_total_sum'] = float(J_tot.sum())
            row['W_total_cadence_sum'] = float(W_tot_cadence.sum())
        rec.append(row)

    def checkpoint(tt):
        r = np.linalg.norm(x, axis=1)
        L = x[:, 0]*v[:, 1] - x[:, 1]*v[:, 0]
        cps[f'{tt/period:g}'] = dict(
            t_periods=float(tt/period), positions=x.tolist(), velocities=v.tolist(),
            radius=r.tolist(), L_ratio=(L/L0).tolist(),
            eps_residual=(eps(x, v) - eps0 + I_dC).tolist(),
            W_total=W_tot.tolist(), J_total=J_tot.tolist(),
            W_band={k: vv.tolist() for k, vv in W_band.items()},
            J_band={k: vv.tolist() for k, vv in J_band.items()})

    record(0.)
    last_diag = diagnostics(0.)
    next_rec, next_diag = record_every*period, diag_every*period
    while t < t_end - 1e-12:
        r = np.linalg.norm(x, axis=1)
        if r.min() < r_min:
            status = FM.UNRESOLVED_CENTRE
            break
        h = min(h_max, eta*float(np.min(np.sqrt(r**3/GM))), t_end - t)
        for ct in cp_times:                                   # land exactly on checkpoint times
            if t < ct - 1e-12 < t + h:
                h = ct - t
        v_half = v + .5*h*a
        x_mid = x + .5*h*v_half
        x_new = x + h*v_half
        if np.max(np.abs(x_new)) > .9*half_width:
            status = FM.LEFT_DOMAIN
            break
        x = x_new
        if memory:
            field.advance(x_mid, rates, h)
        a = accel(x)
        v = v_half + .5*h*a
        t += h
        if memory:
            dCdt_new = _dC_dt(field, x, rates, x)
            I_dC += .5*h*(dCdt + dCdt_new)
            dCdt = dCdt_new
            g_new = field.sample_gradient(x)
            P_new = np.sum(v*g_new, axis=1)
            W_tot += .5*h*(P + P_new)
            Q_new = x[:, 0]*g_new[:, 1] - x[:, 1]*g_new[:, 0]
            J_tot += .5*h*(Q + Q_new)             # each end paired with its OWN position and gradient
            P, Q, g = P_new, Q_new, g_new
            if band_diagnostics and t >= next_diag - 1e-12:
                d = diagnostics(t)
                if last_diag is not None:
                    dt = d[0] - last_diag[0]
                    for k in W_band:
                        W_band[k] += .5*dt*(last_diag[1][k] + d[1][k])
                        J_band[k] += .5*dt*(last_diag[2][k] + d[2][k])
                    W_tot_cadence += .5*dt*(last_diag[3] + d[3])
                last_diag = d
                next_diag = t + diag_every*period
        if t >= next_rec - 1e-12:
            record(t)
            next_rec = t + record_every*period
        if any(abs(t - ct) < 1e-9 for ct in cp_times):
            checkpoint(t)
    if not rec or abs(rec[-1]['t_periods'] - t/period) > 1e-9:
        record(t)
    final = dict(t_periods=float(t/period), positions=x.tolist(), velocities=v.tolist())
    return dict(status=status, t_final_periods=float(t/period), records=rec, checkpoints=cps,
                termination_state=final,
                band_sum_check_worst=float(max(band_check)) if band_check else None)

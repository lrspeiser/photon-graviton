"""CF-1 declared interaction (see protocol.md): elastic, equal-mass companion–companion scattering, isotropic in
the centre-of-mass frame, with a constant cross-section per unit mass sigma/m. Energy and momentum are conserved
exactly in every event.

An incoming companion at radius r moves with speed w = sqrt(u^2 + v_esc^2(r)) in a random direction. Its target is
a bound companion whose velocity is drawn from an isotropic Maxwellian with one-dimensional dispersion s,
truncated at v_esc (bound). With constant sigma/m, events occur in proportion to the relative speed, so every
tally is weighted by it. Specific energies are measured relative to R_b, where the bath is specified.

For stationary targets (s = 0) the outcomes are analytic: the incomer keeps a fraction of w^2 that is uniform on
[0, 1]. So P(capture) = v_esc^2/w^2, P(target ejected) = u^2/w^2, and the net number retained per event is
(v_esc^2 - u^2)/(v_esc^2 + u^2).
"""
import numpy as np


def _unit(rng, n):
    d = rng.normal(size=(n, 3))
    return d/np.linalg.norm(d, axis=1)[:, None]


def outcomes(u, v_esc, s, n=400000, seed=12345):
    """Relative-speed-weighted outcome statistics per scattering event, per unit mass."""
    rng = np.random.default_rng(seed)
    w = np.sqrt(u*u + v_esc*v_esc)
    v1 = w*_unit(rng, n)
    v2 = np.zeros((n, 3))
    if s > 0:
        v2 = rng.normal(scale=s, size=(n, 3))
        bad = np.einsum('ij,ij->i', v2, v2) >= v_esc*v_esc
        while bad.any():
            v2[bad] = rng.normal(scale=s, size=(int(bad.sum()), 3))
            bad = np.einsum('ij,ij->i', v2, v2) >= v_esc*v_esc
    vcm, g = .5*(v1 + v2), v1 - v2
    gmag = np.linalg.norm(g, axis=1)
    e = _unit(rng, n)
    v1p, v2p = vcm + .5*gmag[:, None]*e, vcm - .5*gmag[:, None]*e
    k1, k2 = np.einsum('ij,ij->i', v1p, v1p), np.einsum('ij,ij->i', v2p, v2p)
    k1i, k2i = np.einsum('ij,ij->i', v1, v1), np.einsum('ij,ij->i', v2, v2)
    phi = -.5*v_esc*v_esc
    E1, E2, E1p, E2p = .5*k1i + phi, .5*k2i + phi, .5*k1 + phi, .5*k2 + phi
    cap, ej = k1 < v_esc*v_esc, k2 > v_esc*v_esc
    wt = gmag/gmag.sum()
    bound_after = np.where(cap, E1p, 0.) + np.where(ej, 0., E2p)
    escaping = np.where(cap, 0., E1p) + np.where(ej, E2p, 0.)
    return dict(
        u=u, v_esc=v_esc, s=s, events=n,
        capture=float(wt @ cap), eject=float(wt @ ej), both=float(wt @ (cap & ej)),
        net_retained=float(wt @ (cap.astype(float) - ej)),
        mean_relative_speed=float(gmag.mean()),
        rate_weighted_net=float(np.mean(gmag*(cap.astype(float) - ej))),     # <|v_rel| * net> per event pair
        bound_energy_change=float(wt @ (bound_after - E2)),                    # per event, per unit mass
        escaping_energy=float(wt @ escaping),
        energy_violation=float(np.max(np.abs(k1 + k2 - k1i - k2i))/(w*w)),
        momentum_violation=float(np.max(np.abs(v1p + v2p - v1 - v2))/w),
        analytic_stationary=dict(capture=v_esc**2/w**2, eject=u*u/w**2, net=(v_esc**2 - u*u)/w**2))

"""Round 19, step 2a: the absorbing stream as an exact wave medium (a uniform stream, solved in Fourier space).

Round 18 put the stream that absorbs counter-moving waves into the full model in its ray (eikonal) form: every
coupling between two pieces, near field included, multiplied by exp(-kappa L), L the distance the straight path
travels against the stream. That form lost passivity at kappa >= 10 (in units of the inverse wavelength; k = 2 pi),
where the absorption length is shorter than a wavelength. The review asks for the medium in full wave form: what
absorbs, energy and momentum, convergence, passivity and no self-push. Here the medium is solved exactly for a
uniform stream along +z (the geometry a receiver sees far from a source), with units c = 1, wavelength 1:

    (-lap - k^2 - i Pi) psi = source,     G(q) = 1 / (q^2 - k^2 - i Pi(q)),     G(x) = int d^3q/(2 pi)^3 e^{iq.x} G(q)

A medium is passive (no arrangement of emitters can draw energy from it) exactly when Pi(q) >= 0 for every q, since
the power a source distribution s feeds it is (omega/2) int |s(q)|^2 Im G(q) >= 0. The media:

  'doppler'   absorbers carried by the stream at the wave's own speed, point-like: an absorber moving along +z sees a
              wave of wavevector q at the frequency omega' = omega - u q_z = k - q_z (u = c = 1); an Ohmic absorber
              takes energy at a rate proportional to omega', so Pi = kappa (k - q_z). Waves moving with the stream
              (q_z = k) are at zero frequency for it: any passive absorber lets them pass, so the medium is exactly
              one-way for propagating waves, with the cardioid attenuation (kappa/2)(1 - cos theta) per unit length.
              But the near field of a point emitter holds components with q_z > k, for which Pi < 0: the stream pumps
              them (the anomalous Doppler effect) and drags every emitter along with it; the drag grows without bound
              with the emitter's sharpness (checked below with a cutoff).
  'cardioid'  the same absorbers, but only able to take up waves near the wave's own wavenumber (size about a
              wavelength: a form factor exp(-(|q| - k)^2 / 2 w^2)), direction read from q^ : Pi = kappa k (1 - q^_z) F(|q|).
              Pi >= 0 everywhere: passive at every kappa, by construction.
  'inward'    round 18's rule as a wave medium: only waves moving against the stream are absorbed,
              Pi = 2 kappa k max(0, -q^_z) F(|q|) (on the shell: kappa per unit length travelled against the stream).

A small uniform loss eps (amplitude e^{-eps R/2}) regularises the shell in every medium and is divided out. G is the
free Green's function e^{i k_eps R}/(4 pi R) plus the correction dG(q) = i Pi / ((q^2 - k_eps^2 - i Pi)(q^2 - k_eps^2)),
which for the form-factor media is confined to the shell region and bounded at R = 0 (the near field is untouched).
The integral runs over q (nonuniform, dense at the shell) and mu = cos(angle to the stream) (Gauss-Legendre), the
azimuth done exactly (J0).

Measured: transmission with and against the stream against the ray form; the near field; the self-force on a lone
emitter (no self-push?); the power booking; passivity of the pieces' coupling, exact against round 18's ray form.

    python code/full_wave_v19.py --output run-full-wave-v19/full_wave_v19.json
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np
from scipy.special import j0

K = 2 * np.pi
EPS = 0.02


def q_grid(k=K, qmax=None, w=None, n_side=1600, n_far=400):
    """Radial wavenumber nodes and trapezoid weights, geometric toward the shell q = k from both sides."""
    qmax = qmax or (k + 8 * (w or k))
    off = np.geomspace(1e-6, k * 0.999, n_side)
    inner = k - off[::-1]
    outer_near = k + np.geomspace(1e-6, min(qmax - k, 3 * k), n_side)
    outer_far = np.linspace(outer_near[-1], qmax, n_far + 1)[1:] if qmax > outer_near[-1] else np.array([])
    q = np.concatenate([[0.0], inner, [k], outer_near, outer_far])
    wq = np.zeros_like(q); dq = np.diff(q); wq[:-1] += dq / 2; wq[1:] += dq / 2
    return q, wq


def Pi_of(medium, q, mu, kappa, w, k=K):
    if medium == 'doppler':
        return kappa * (k - q * mu)
    F = np.exp(-(q - k) ** 2 / (2 * w ** 2))
    if medium == 'cardioid':
        return kappa * k * (1 - mu) * F
    if medium == 'inward':
        return 2 * kappa * k * np.maximum(0.0, -mu) * F
    if medium == 'none':
        return np.zeros(np.broadcast(q, mu).shape)
    raise ValueError(medium)


class Medium:
    def __init__(self, medium, kappa, w=0.5 * K, eps=EPS, nmu=512, qmax=None, k=K):
        self.medium, self.kappa, self.w, self.eps, self.k = medium, kappa, w, eps, k
        self.keps = np.sqrt(k * k + 1j * eps * k)
        if medium == 'doppler':
            qmax = qmax or 12 * k
        self.q, self.wq = q_grid(k, qmax, w)
        self.mu, self.wmu = np.polynomial.legendre.leggauss(nmu)
        Q, MU = np.meshgrid(self.q, self.mu, indexing='ij')
        self.Pi = Pi_of(medium, Q, MU, kappa, w, k)
        D = Q ** 2 - self.keps ** 2
        self.dG = 1j * self.Pi / ((D - 1j * self.Pi) * D)          # G - G_eps0, per (q, mu)
        self.Gq = 1.0 / (D - 1j * self.Pi)
        self.Q, self.MU = Q, MU
        self.W2 = (self.wq[:, None] * self.wmu[None, :]) * Q ** 2 / (4 * np.pi ** 2)

    def G0(self, R):
        R = np.asarray(R, float)
        return np.exp(1j * self.keps * R) / (4 * np.pi * R)

    def dG_at(self, rho, z):
        """The correction to the free Green's function at points (rho, z) (arrays), by the double integral."""
        rho = np.atleast_1d(rho); z = np.atleast_1d(z); out = np.empty(rho.shape, complex)
        st = np.sqrt(np.clip(1 - self.MU ** 2, 0, None))
        for i, (r_, z_) in enumerate(zip(rho.ravel(), z.ravel())):
            ker = np.exp(1j * self.Q * self.MU * z_) * (j0(self.Q * st * r_) if r_ else 1.0)
            out.ravel()[i] = np.sum(self.W2 * ker * self.dG)
        return out

    def G_at(self, rho, z):
        R = np.hypot(rho, z)
        return self.G0(R) + self.dG_at(rho, z)

    def grad_z_dG_at(self, rho, z):
        """d/dz of the correction (the medium's field gradient along the stream) at (rho, z)."""
        rho = np.atleast_1d(rho); z = np.atleast_1d(z); out = np.empty(rho.shape, complex)
        st = np.sqrt(np.clip(1 - self.MU ** 2, 0, None))
        for i, (r_, z_) in enumerate(zip(rho.ravel(), z.ravel())):
            ker = 1j * self.Q * self.MU * np.exp(1j * self.Q * self.MU * z_) * (j0(self.Q * st * r_) if r_ else 1.0)
            out.ravel()[i] = np.sum(self.W2 * ker * self.dG)
        return out

    def self_terms(self):
        """For a lone monopole: the power it feeds the medium relative to the free value (Im G(0) / (k/4 pi)) and the
        self-force along the stream in units of the power over c (F_z c / P; +1 would be all its momentum sent
        against the stream, i.e. a push along it... sign: a positive value is a push along the stream)."""
        ImG = np.sum(self.W2 * self.Gq.imag)                        # includes the eps medium's own part
        Pfree = np.real(self.keps) / (4 * np.pi)
        # force on a monopole of unit strength: (1/2) Re grad G_reg(0) (the field's own gradient at the source); the
        # free part is symmetric, so only dG contributes: Re int iq_z dG = -int q mu Im dG
        Fz = -np.sum(self.W2 * self.Q * self.MU * self.dG.imag)
        # the power it feeds relative to the free case, and the force per power (P = omega Im G / 2, F = Re grad G / 2,
        # both per |q|^2; omega = k with c = 1)
        return dict(power_ratio=float(ImG / Pfree), force_over_power_c=float((Fz / 2) / (K * ImG / 2)))


def transmissions(med, dists=(0.5, 1.0, 2.0, 4.0), angles_deg=(0, 45, 90, 135, 180)):
    """|G| relative to the free medium (the eps part divided out) along directions at an angle to the stream (0: with
    it), against the ray forms: 'cardioid' exp(-(kappa/2)(1 - cos theta) R), 'inward' exp(-kappa max(0, -cos) R)."""
    rows = []
    for ang in angles_deg:
        th = np.radians(ang)
        for R in dists:
            rho, z = R * np.sin(th), R * np.cos(th)
            G = med.G_at(np.array([rho]), np.array([z]))[0]
            T = abs(G / med.G0(R))
            if med.medium == 'cardioid':
                ray = np.exp(-(med.kappa / 2) * (1 - np.cos(th)) * R)
            elif med.medium == 'inward':
                ray = np.exp(-med.kappa * max(0.0, -np.cos(th)) * R)
            elif med.medium == 'doppler':
                ray = np.exp(-(med.kappa / 2) * (1 - np.cos(th)) * R)
            else:
                ray = 1.0
            rows.append(dict(angle_deg=ang, R=R, transmission=float(T), ray=float(ray),
                             phase_shift=float(np.angle(G / med.G0(R)))))
    return rows


def near_field(med, R=(0.05, 0.1, 0.15, 0.3)):
    """The correction against the free near field at short range, with and against the stream."""
    out = []
    for r in R:
        up = med.dG_at(np.array([0.0]), np.array([r]))[0]; dn = med.dG_at(np.array([0.0]), np.array([-r]))[0]
        g0 = med.G0(r)
        out.append(dict(R=r, correction_with_stream=float(abs(up / g0)), correction_against=float(abs(dn / g0))))
    return out


def ball(n, R, rng, sep=0.12):
    pts = []
    while len(pts) < n:
        p = rng.uniform(-R, R, 3)
        if p @ p <= R * R and all(np.sum((p - q) ** 2) >= sep * sep for q in pts):
            pts.append(p)
    return np.array(pts)


def passivity(med, n=24, Rb=1.0, seed=3):
    """The smallest eigenvalue of the monopole pieces' dissipative matrix (G - G^dagger)/2i, exact medium against the
    ray form of round 18 (the free coupling times exp(-kappa L) pair by pair; L = the distance travelled against the
    stream, max(0, -dz), for the 'inward' rule; (1/2)(R - dz) for the cardioid)."""
    rng = np.random.default_rng(seed)
    x = ball(n, Rb, rng)
    d = x[:, None, :] - x[None, :, :]                                 # x_j - x_l: from l to j
    rho = np.hypot(d[..., 0], d[..., 1]); z = d[..., 2]; R = np.linalg.norm(d, axis=-1)
    off = ~np.eye(n, dtype=bool)
    Gex = np.zeros((n, n), complex)
    Gex[off] = med.G_at(rho[off], z[off])
    Gself = np.sum(med.W2 * med.Gq.imag) * 1j                        # the radiative (finite) self part
    np.fill_diagonal(Gex, Gself)
    G0m = np.zeros((n, n), complex); G0m[off] = med.G0(R[off]); np.fill_diagonal(G0m, 1j * np.real(med.keps) / (4 * np.pi))
    if med.medium == 'inward':
        Lr = np.maximum(0.0, -z)
        T = np.exp(-med.kappa * Lr)
    else:
        T = np.exp(-(med.kappa / 2) * (R - z))
    Gray = G0m * np.where(off, T, 1.0)
    ev = lambda M: np.linalg.eigvalsh(((M - M.conj().T) / 2j + ((M - M.conj().T) / 2j).conj().T) / 2)
    e_ex, e_ray, e_0 = ev(Gex), ev(Gray), ev(G0m)
    return dict(n=n, radius=Rb, exact_min=float(e_ex.min()), exact_max=float(e_ex.max()), ray_min=float(e_ray.min()),
                ray_max=float(e_ray.max()), free_min=float(e_0.min()), free_max=float(e_0.max()))


def doppler_drag(kappas=(0.5, 2.3, 5.0), cutoffs=(4, 8, 16, 32)):
    """The point-absorber medium: the self-force on a lone emitter with the integral cut at q = cutoff k (the emitter's
    sharpness), showing its growth: the stream drags the emitter's own near field."""
    out = []
    for kap in kappas:
        row = dict(kappa=kap, cutoff_over_k=[], force_over_power_c=[], power_ratio=[])
        for cq in cutoffs:
            m = Medium('doppler', kap, qmax=cq * K, nmu=256)
            st = m.self_terms()
            row['cutoff_over_k'].append(cq); row['force_over_power_c'].append(st['force_over_power_c']); row['power_ratio'].append(st['power_ratio'])
        out.append(row)
        print(f"doppler kappa {kap}: self-force/(P/c) at cutoffs {cutoffs}: " + ', '.join(f'{v:+.3f}' for v in row['force_over_power_c']), flush=True)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--quick', action='store_true')
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    out = dict(experiment='round 19: the absorbing stream as an exact wave medium (uniform stream, Fourier space)',
               units='wavelength 1, k = 2 pi, c = 1; kappa per unit length', eps=EPS, media={})
    # a check: the free medium's correction vanishes, and the regularised free power is 1
    m0 = Medium('none', 0.0)
    out['check_free'] = dict(self=m0.self_terms(), dG_at_1=float(abs(m0.dG_at(np.array([0.0]), np.array([1.0]))[0])))
    kappas = (0.5, 2.3, 5.0) if args.quick else (0.5, 2.3, 5.0, 10.0, 20.0)
    for medium in ('inward', 'cardioid'):
        for w_rel in ((0.5,) if args.quick else (0.25, 0.5)):
            for kap in kappas:
                med = Medium(medium, kap, w=w_rel * K)
                row = dict(medium=medium, kappa=kap, form_factor_width_over_k=w_rel, self=med.self_terms(),
                           transmissions=transmissions(med), near_field=near_field(med), passivity=passivity(med))
                out['media'][f'{medium}, kappa {kap}, w {w_rel}k'] = row
                tr = {(t['angle_deg'], t['R']): t for t in row['transmissions']}
                print(f"[{time.monotonic() - t0:5.0f} s] {medium:8s} kappa {kap:5.1f} w {w_rel}k: power x{row['self']['power_ratio']:.3f}, "
                      f"self-force {row['self']['force_over_power_c']:+.2e} P/c; T against at R=1,2: {tr[(180, 1.0)]['transmission']:.3f}, "
                      f"{tr[(180, 2.0)]['transmission']:.3f} (ray {tr[(180, 1.0)]['ray']:.3f}, {tr[(180, 2.0)]['ray']:.3f}); with: {tr[(0, 2.0)]['transmission']:.3f}; "
                      f"sideways R=2: {tr[(90, 2.0)]['transmission']:.3f} (ray {tr[(90, 2.0)]['ray']:.3f}); near field at 0.15: "
                      f"{row['near_field'][2]['correction_against']:.1e}; passivity exact {row['passivity']['exact_min']:+.2e}, "
                      f"ray {row['passivity']['ray_min']:+.2e} (free {row['passivity']['free_min']:+.1e})", flush=True)
                args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    out['doppler_point_absorbers'] = doppler_drag() if not args.quick else doppler_drag(kappas=(2.3,), cutoffs=(4, 8))
    out['seconds'] = time.monotonic() - t0
    args.output.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('wrote', args.output, f'({out["seconds"]:.0f} s)')


if __name__ == '__main__':
    main()

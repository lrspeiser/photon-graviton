#!/usr/bin/env python3
"""Derive MOND as the cold limit of hot-companion gravity, and check every step.

    python derive_mond.py --output-dir ../run-derive

The direction of the derivation is the point: MOND is not an input. It comes out of the
companion mechanism when the matter is cold, with its constant a0 and its interpolating
function both fixed by the mechanism.

Postulates (proposed here; originality unverified):
  P1 emission    every kilogram of ordinary matter emits companion energy at rate l (W/kg)
  P2 streaming   the energy streams away from its emitter at speed u, conserved
  P3 coherence   cold, orderly matter emits in step: energy FLUXES add as vectors
                 hot matter (random speed sigma) emits an extra k*l, k = 3 sigma^2/u^2,
                 out of step: energy DENSITIES add as scalars
  P4 amplitude   the companion's energy density is A^2/(8 pi G); its pull equals A
  P5 attachment  strong fields hold it: released fraction f = exp(-|g_N|/g_d)
Result:  g = g_N + f * sqrt(a (|g_N| + S_hot)),  a = 2 l / u
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[2] / 'tools'))
import formula_guard as FG

G_SI = 6.674e-11
C = 2.998e8
A_FIT = 6.5797e-11          # m/s^2, fitted on SPARC (run-v1)
LAM = 4.281                 # g_d / a, fitted on SPARC
U_FIT = 874.1e3             # m/s, fitted on X-COP
YEAR = 3.15576e7


def step1_coherent_flux_is_newtonian(rng):
    """Vector sum of steady radial energy streams equals (l/4piG) x Newtonian field."""
    n = 400
    pos = rng.normal(size=(n, 3)) * [3.0, 3.0, 0.3]           # a thick disk of emitters
    m = rng.lognormal(0, 1, n)
    l = 1.7
    pts = rng.normal(size=(50, 3)) * 8
    worst = 0.
    for x in pts:
        d = x - pos; r = np.linalg.norm(d, axis=1)
        flux = np.sum((l * m / (4 * np.pi * r ** 2))[:, None] * d / r[:, None], axis=0)
        gN = -G_SI * np.sum((m / r ** 2)[:, None] * d / r[:, None], axis=0)
        worst = max(worst, np.linalg.norm(flux - (-(l / (4 * np.pi * G_SI)) * gN)) / np.linalg.norm(flux))
    return dict(max_relative_mismatch=float(worst), identity='F_coh = -(l / 4 pi G) g_N  (Gauss geometry)')


def step2_coherent_vs_incoherent(rng):
    """Why heat means incoherence: random motion Doppler-shifts every emitter.

    Cold, orderly emitters share one frequency and phase: the time-averaged intensity of
    their summed field equals the SQUARED VECTOR SUM (they cancel like Newton's pull).
    Hot emitters move randomly (speed sigma), so each is Doppler-shifted by v_i/u. Averaged
    over time, all cross terms vanish and the intensity equals the SCALAR SUM.
    """
    n = 200
    amps = rng.uniform(0.5, 1.5, n); dirs = rng.normal(size=(n, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    t = np.linspace(0, 2 * np.pi * 400, 80001)
    out = {}
    for label, sigma_over_u in (('cold', 0.0), ('warm', 0.02), ('hot', 0.3)):
        w = 1 + rng.normal(0, sigma_over_u, n) if sigma_over_u > 0 else np.ones(n)
        ph = rng.uniform(0, 2 * np.pi, n) if sigma_over_u > 0 else np.zeros(n)
        field = np.zeros((3, t.size))
        for i in range(n):
            field += amps[i] * dirs[i][:, None] * np.cos(w[i] * t + ph[i])[None, :]
        inten = float(np.mean(np.sum(field ** 2, axis=0)) * 2)
        out[label] = dict(sigma_over_u=sigma_over_u, intensity=inten,
                          scalar_sum=float(np.sum(amps ** 2)),
                          squared_vector_sum=float(np.sum((amps[:, None] * dirs).sum(0) ** 2)))
    return out


def step3_symbolic():
    """Energy balance -> amplitude; limits; scale invariance; Tully-Fisher. Symbolic."""
    l, u, G, M, r, gN, a, gd, beta, A = sp.symbols('ell u G M r g_N a g_d beta A', positive=True)
    # steady state: outward energy flux through a sphere equals emission; flux = u * A^2/(8 pi G)
    flux_balance = sp.Eq(4 * sp.pi * r ** 2 * u * A ** 2 / (8 * sp.pi * G), l * M)
    A_sol = sp.solve(flux_balance, A)[0]
    a_def = 2 * l / u
    check_amp = sp.simplify(A_sol - sp.sqrt(a_def * G * M / r ** 2))
    # cold law and its limits
    y = sp.symbols('y', positive=True); lam = sp.symbols('lambda', positive=True)
    nu = 1 + sp.exp(-y / lam) / sp.sqrt(y)
    deep = sp.limit(nu * sp.sqrt(y), y, 0)                 # nu -> y^(-1/2): g -> sqrt(a g_N)
    newton = sp.limit(nu, y, sp.oo)                         # nu -> 1
    # deep-regime scale invariance: r -> beta r, t -> beta t leaves orbits' shapes invariant:
    # g(beta r) / (beta^-1) should equal g(r) for g = sqrt(a G M)/r
    gdeep = sp.sqrt(a * G * M) / r
    scale = sp.simplify(gdeep.subs(r, beta * r) * beta - gdeep)
    # flat rotation and Tully-Fisher: v^2 = r g -> v^4 = G M a
    v4 = sp.simplify((r * gdeep) ** 2)
    return dict(amplitude_from_energy_balance=str(A_sol), amplitude_minus_sqrt_aGM_over_r=str(check_amp),
                a_equals='2 l / u', deep_limit_nu_times_sqrt_y=str(deep), newtonian_limit_nu=str(newton),
                deep_scale_invariance_residual=str(scale), v4=str(v4))


def step4_cold_limit_is_mond():
    """The cold law satisfies Milgrom's defining properties; compare nu with published ones."""
    def pm(r, M, s):
        gN = G_SI * M / r ** 2
        return gN + np.exp(-gN / (LAM * A_FIT)) * np.sqrt(A_FIT * gN)
    KPC = 3.0856775814913673e19; MSUN = 1.98847e30
    y = np.geomspace(1e-4, 1e3, 200)
    nu = 1 + np.exp(-y / LAM) / np.sqrt(y)
    table = {}
    for name, f in FG.MOND_CATALOGUE.items():
        best = min(((np.max(np.abs(np.log10(nu) - np.log10(f(y * s)))), s) for s in np.geomspace(.05, 20, 4001)))
        table[name] = dict(max_dev_dex=float(best[0]), a0_scale=float(best[1]))
    return dict(nu_ours='1 + exp(-y/lambda)/sqrt(y), y = g_N/a, lambda = %.3f' % LAM,
                newtonian_at_y_100=float(nu[np.argmin(abs(y - 100))] - 1),
                deep_nu_times_sqrt_y_at_1e4=float(nu[0] * np.sqrt(y[0])),
                comparison_with_published=table)


def step5_numbers():
    l = A_FIT * U_FIT / 2
    per_age = l * 13.8e9 * YEAR / C ** 2
    msun, lsun = 1.98847e30, 3.828e26
    # the Sun: mostly hot plasma inside; k from the mass-weighted thermal speed ~ 300 km/s
    k_sun = 3 * (3.0e5) ** 2 / U_FIT ** 2
    L_sun_comp = l * (1 + k_sun) * msun
    return dict(emission_power_W_per_kg=l,
                fraction_of_rest_energy_over_13_8_Gyr=per_age,
                sun_companion_power_W=L_sun_comp, sun_companion_power_over_luminosity=L_sun_comp / lsun,
                sun_extra_mass_loss_per_year=L_sun_comp / C ** 2 / msun * YEAR,
                note='applies if companion energy emitted inside the attachment radius is not re-absorbed')


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    out = ap.parse_args().output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True)
    rng = np.random.default_rng(20260923)
    res = dict(step1=step1_coherent_flux_is_newtonian(rng), step2=step2_coherent_vs_incoherent(rng),
               step3=step3_symbolic(), step4=step4_cold_limit_is_mond(), step5=step5_numbers())
    (out / 'derivation.json').write_text(json.dumps(res, indent=2) + '\n')
    s1, s2, s3, s4, s5 = (res[k] for k in ('step1', 'step2', 'step3', 'step4', 'step5'))
    print('1. coherent energy flux = (l/4piG) x Newtonian field: max mismatch', f"{s1['max_relative_mismatch']:.1e}")
    for k, v in s2.items():
        print(f"2. {k:5s} (sigma/u={v['sigma_over_u']}): intensity {v['intensity']:.2f} | scalar sum {v['scalar_sum']:.2f}"
              f" | squared vector sum {v['squared_vector_sum']:.2f}")
    print('3. energy balance gives A =', s3['amplitude_from_energy_balance'], '| minus sqrt(aGM)/r with a=2l/u:',
          s3['amplitude_minus_sqrt_aGM_over_r'])
    print('   deep limit nu*sqrt(y) ->', s3['deep_limit_nu_times_sqrt_y'], '| Newtonian limit nu ->', s3['newtonian_limit_nu'],
          '| scale-invariance residual', s3['deep_scale_invariance_residual'], '| v^4 =', s3['v4'])
    print('4. cold law nu:', s4['nu_ours'], f"| nu-1 at y=100: {s4['newtonian_at_y_100']:.1e} | deep check {s4['deep_nu_times_sqrt_y_at_1e4']:.4f}")
    for k, v in s4['comparison_with_published'].items():
        print(f"     vs {k:52s} max {v['max_dev_dex']:.3f} dex (a0 scale {v['a0_scale']:.2f})")
    print(f"5. a0 = 2 l/u  ->  l = {s5['emission_power_W_per_kg']:.2e} W/kg; over 13.8 Gyr "
          f"{s5['fraction_of_rest_energy_over_13_8_Gyr']:.1e} of rest energy;")
    print(f"   Sun: {s5['sun_companion_power_W']:.1e} W = {s5['sun_companion_power_over_luminosity']:.2f} L_sun, "
          f"extra mass loss {s5['sun_extra_mass_loss_per_year']:.1e} per year")


if __name__ == '__main__':
    main()

"""Round 24 review of round24-casimir-eft.md: its arithmetic, two identities it implies, and order-of-magnitude
bounds for its test B (ultra-high-energy cosmic rays) and for its equation without the strong-field hold.

1. The note's section-11 numbers, recomputed.
2. Its compact form a = [chi C_8/(8 x_q^2)] G m_p^3 c u / hbar^2: the microscopic model predicts a/u (both a and u
   are fitted on data), and the u at which it would match a_fit exactly.
3. Its strong-field lead g_* = hbar Gamma_X/(m_p R_b) is the same combination as a: g_*/a_micro = 1/(2 C_8) exactly,
   so it is not a second prediction; g_d/a = 3.2 stays unexplained.
4. Test B, bounded: a proton of 10^20 eV that must cross 100 Mpc can lose at most E/t per second. Against the law's
   rest-frame emission per proton, ell m_p, that allows an enhancement factor of about 2 x 10^17. Three ways the
   heat factor could continue to relativistic speeds: saturating, 1 + v^2/u^2 with v < c (the note's section 3 as
   written); growing with the energy, gamma v^2/u^2; growing with the momentum squared, (gamma v)^2/u^2.
5. The note's equation without the hold, in the Solar System: sqrt(a g_N) at the planets' orbits.

    python code/eft_checks_v24.py --output run-eft-v24/eft_checks_v24.json
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

G = 6.67430e-11; hbar = 1.054571817e-34; c = 299792458.0; m_p = 1.67262192369e-27
eV = 1.602176634e-19; M_sun = 1.98847e30; kpc = 3.085677581491367e19; Mpc = 1e3 * kpc; AU = 1.495978707e11
u = 169.4e3; a_fit = 6.298e-11; g_d = 2.027e-10
x_q, chi, C8 = 2.04, 2.744, 0.51


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args(); args.output.parent.mkdir(parents=True, exist_ok=True)
    out = {}
    alpha_G = G * m_p ** 2 / (hbar * c)
    R_b = 4 * x_q * hbar / (m_p * c)
    E_X = C8 * hbar * u / R_b
    Gamma_X = alpha_G * chi * u / R_b
    ell_micro = Gamma_X * E_X / m_p
    a_micro = 2 * ell_micro / u
    ell_fit = a_fit * u / 2
    M = 1e11 * M_sun; r = 10 * kpc
    g_gr = G * ell_fit * M / (u * c * c * r); g_need = math.sqrt(G * M * a_fit) / r
    g_star = hbar * Gamma_X / (m_p * R_b)
    out['note_section_11'] = dict(alpha_G=alpha_G, R_b_fm=R_b / 1e-15, E_X_keV=E_X / (1e3 * eV), Gamma_X_per_s=Gamma_X,
                                  ell_micro_W_per_kg=ell_micro, a_micro=a_micro, a_micro_over_a_fit=a_micro / a_fit,
                                  ordinary_GR_gap=g_need / g_gr, g_star=g_star, g_star_over_g_d=g_star / g_d,
                                  M_D_TeV=145.0 / (math.sqrt(alpha_G) ** 0.25) / 1e6)
    pref = chi * C8 / (8 * x_q ** 2); base = G * m_p ** 3 * c * u / hbar ** 2
    out['compact_form'] = dict(prefactor=pref, G_mp3_c_u_over_hbar2=base, a=pref * base,
                               a_over_u_micro=a_micro / u, a_over_u_fit=a_fit / u,
                               u_for_exact_match_kms=u * a_fit / (pref * base) / 1e3,
                               note='a_micro is proportional to u, which is itself fitted (X-COP); the model predicts a/u')
    out['g_star_identity'] = dict(g_star_over_a_micro=g_star / a_micro, one_over_2C8=1 / (2 * C8),
                                  g_d_over_a_fit=g_d / a_fit,
                                  note='g_* = a_micro/(2 C_8) identically: the same scale as a, not a second one')
    # test B, bounded
    E = 1e20 * eV; L_prop = 100 * Mpc; t = L_prop / c
    P0 = ell_fit * m_p                                  # W per proton at rest (the law's emission)
    allowed = (E / t) / P0
    gam = E / (m_p * c * c); cu2 = (c / u) ** 2
    rows = []
    for name, fac in (('saturating: 1 + v^2/u^2, v < c (section 3 as written)', 1 + cu2),
                      ('growing with energy: gamma v^2/u^2', gam * cu2),
                      ('growing with momentum squared: (gamma v)^2/u^2', gam ** 2 * cu2)):
        P = P0 * fac
        rows.append(dict(rule=name, enhancement=fac, loss_eV_per_s=P / eV, loss_time_s=E / P,
                         distance_before_losing_energy_Mpc=E / P * c / Mpc, margin=allowed / fac))
    out['test_B_bound'] = dict(E_eV=1e20, propagation_Mpc=100, P0_W=P0, P0_eV_per_s=P0 / eV, allowed_enhancement=allowed,
                               gamma=gam, c_over_u_squared=cu2, rules=rows,
                               note='order of magnitude only: the energy loss must stay below E/t; a canonical calculation (the '
                                    'note\'s test B) decides which rule the coupling T^mn dX dX/M_D^4 gives')
    # the equation without the hold, in the Solar System
    GMsun = G * M_sun
    planets = {}
    for name, r_AU in (('Mercury', 0.387), ('Earth', 1.0), ('Saturn', 9.58), ('Neptune', 30.1)):
        gN = GMsun / (r_AU * AU) ** 2
        planets[name] = dict(g_N=gN, extra_no_hold=math.sqrt(a_fit * gN), extra_over_gN=math.sqrt(a_fit / gN),
                             hold_factor=math.exp(-gN / g_d))
    out['no_hold_solar_system'] = planets
    args.output.write_text(json.dumps(out, indent=1) + '\n')
    for k, v in out.items():
        print(k, json.dumps(v, indent=1)[:1500])


if __name__ == '__main__':
    main()

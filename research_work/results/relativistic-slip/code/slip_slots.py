"""T0.2, relativistic part: where can a gravitational slip eta != 1 live?

Four checks, each printed and saved to results.json:

  1. Symbolic: the linearised Einstein tensor of a static, spherical metric
         ds^2 = -(1+2 Phi) dt^2 + (1-2 Psi)(dr^2 + r^2 dOmega^2)
     gives, identically for ANY theory,
         M_Phi(<R) - M_Psi(<R) = 4 pi R^3 p_r(R) / c^2
     i.e. lensing and dynamics can differ inside R only if the gravitational
     field carries a radial stress at R.  (The radial Einstein equation; used
     for dark matter by Faber & Visser 2006, MNRAS 372, 136.)
  2. Symbolic: how a matter-coupling slot sets eta.
       conformal   g~ = A(phi)^2 g                   -> eta = -1 (light blind)
       TeVeS-type  g~ = e^{-2phi} g - B U U          -> eta = 1/(2k-1)
       AeST + conformal coupling kappa to its scalar -> eta = (1-kappa)/(1+kappa)
  3. Toy model: a 1e11 Msun Hernquist galaxy under QUMOND, with the slip put in
     through a constant kappa.  The effective stress is computed from the
     metric by finite differences, and the identity in (1) is checked at every
     radius.
  4. Numbers: what every measured eta requires, and the solar-system size of
     the slip.

Usage:  python slip_slots.py --output-dir ../run-v1
"""
import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp

G = 6.674e-11
C = 2.998e8
MSUN = 1.989e30
KPC = 3.0856775814913673e19
A0 = 1.171e-10  # m/s^2, the SPARC-fitted constant (JR-13)


# --------------------------------------------------------------------------
# 1. Linearised Einstein tensor, static spherical metric
# --------------------------------------------------------------------------
def einstein_linear():
    t, r, th, ph, eps = sp.symbols('t r theta phi epsilon', positive=True)
    Phi = sp.Function('Phi')(r)
    Psi = sp.Function('Psi')(r)
    x = [t, r, th, ph]
    g = sp.diag(-(1 + 2 * eps * Phi), (1 - 2 * eps * Psi),
                (1 - 2 * eps * Psi) * r**2,
                (1 - 2 * eps * Psi) * r**2 * sp.sin(th)**2)
    gi = g.inv()
    n = 4
    Gam = [[[sum(gi[a, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b])
                             - sp.diff(g[b, c], x[d])) for d in range(n)) / 2
             for c in range(n)] for b in range(n)] for a in range(n)]

    def ricci(b, c):
        return sp.simplify(sum(
            sp.diff(Gam[a][b][c], x[a]) - sp.diff(Gam[a][b][a], x[c])
            + sum(Gam[a][a][d] * Gam[d][b][c] - Gam[a][c][d] * Gam[d][b][a]
                  for d in range(n)) for a in range(n)))

    R = sp.zeros(n)
    for b in range(n):
        for c in range(n):
            R[b, c] = ricci(b, c)
    Rs = sum(gi[a, b] * R[a, b] for a in range(n) for b in range(n))
    Gt = R - g * Rs / 2

    def first(expr):
        return sp.simplify(sp.series(expr, eps, 0, 2).removeO().coeff(eps, 1))

    # mixed components G^mu_nu at first order = source 8 pi G T^mu_nu
    G_tt = first(-Gt[0, 0] / g[0, 0])         # -G^t_t -> 8 pi G rho
    G_rr = first(Gt[1, 1] / g[1, 1])          # -> 8 pi G p_r
    G_hh = first(Gt[2, 2] / g[2, 2])          # -> 8 pi G p_t
    return r, Phi, Psi, G_tt, G_rr, G_hh


def check_identity():
    r, Phi, Psi, G_tt, G_rr, G_hh = einstein_linear()
    lap = lambda f: sp.diff(r**2 * sp.diff(f, r), r) / r**2
    out = {
        '-G^t_t == 2 lap Psi': sp.simplify(G_tt - 2 * lap(Psi)) == 0,
        'G^r_r == (2/r)(Phi - Psi)\'': sp.simplify(
            G_rr - 2 * sp.diff(Phi - Psi, r) / r) == 0,
    }
    # conservation d p_r/dr + 2 (p_r - p_t)/r = 0 must hold identically
    cons = sp.diff(G_rr, r) + 2 * (G_rr - G_hh) / r
    out['static conservation holds identically'] = sp.simplify(cons) == 0
    # the mass identity: with rho = G_tt/8piG, p_r = G_rr/8piG, and
    # M_X(<R) = R^2 X'(R)/G:  M_Phi - M_Psi = R^2 (Phi-Psi)'/G = 4 pi R^3 p_r
    Gs = sp.symbols('G', positive=True)
    lhs = r**2 * sp.diff(Phi - Psi, r) / Gs
    rhs = 4 * sp.pi * r**3 * G_rr / (8 * sp.pi * Gs)
    out['M_Phi - M_Psi == 4 pi R^3 p_r'] = sp.simplify(lhs - rhs) == 0
    out['G^r_r'] = str(G_rr)
    out['G^theta_theta'] = str(sp.simplify(G_hh))
    return out


# --------------------------------------------------------------------------
# 2. Coupling slots
# --------------------------------------------------------------------------
def coupling_slots():
    eps, Phi, Psi, phi, k, kap = sp.symbols('epsilon Phi Psi phi k kappa')
    res = {}
    # conformal, A = exp(kappa phi): read off Phi~, Psi~ at first order
    A2 = sp.exp(2 * kap * eps * phi)
    g00 = sp.expand(sp.series(-A2 * (1 + 2 * eps * Phi), eps, 0, 2).removeO())
    gij = sp.expand(sp.series(A2 * (1 - 2 * eps * Psi), eps, 0, 2).removeO())
    Phit = sp.simplify(-(g00 + 1) / 2 / eps)
    Psit = sp.simplify(-(gij - 1) / 2 / eps)
    res['conformal'] = {'Phi~': str(Phit), 'Psi~': str(Psit),
                        'eta of the conformal part':
                        str(sp.simplify((Psit - Psi) / (Phit - Phi)))}
    # TeVeS-type disformal: g~ = e^{-2 phi} g - 2 k sinh(2 phi) U U,
    # static U_0 U_0 = (1 + 2 Phi)
    g00 = -(1 + 2 * eps * Phi) * (sp.exp(-2 * eps * phi)
                                  + 2 * k * sp.sinh(2 * eps * phi))
    gij = sp.exp(-2 * eps * phi) * (1 - 2 * eps * Psi)
    g00 = sp.expand(sp.series(g00, eps, 0, 2).removeO())
    gij = sp.expand(sp.series(gij, eps, 0, 2).removeO())
    Phit = sp.simplify(-(g00 + 1) / 2 / eps)
    Psit = sp.simplify(-(gij - 1) / 2 / eps)
    res['disformal (TeVeS k=1)'] = {
        'Phi~': str(Phit), 'Psi~': str(Psit),
        'eta of the scalar part': str(sp.simplify((Psit - Psi) / (Phit - Phi)))}
    # AeST quasi-static gives Phi = Psi = PhiN + PhiA (Skordis & Zlosnik 2021).
    # Add a conformal coupling of matter to the MOND scalar, carried here as
    # kappa * PhiA.  Matter-frame extra potentials:
    PhiA = sp.symbols('Phi_A')
    Phix = PhiA + kap * PhiA
    Psix = PhiA - kap * PhiA
    res['AeST + conformal kappa'] = {
        'extra in dynamics': str(Phix), 'extra in spatial metric': str(Psix),
        'eta': str(sp.simplify(Psix / Phix)),
        'lensing weight (1+eta)/2': str(sp.simplify((1 + Psix / Phix) / 2))}
    return res


# --------------------------------------------------------------------------
# 3. Toy model: Hernquist galaxy, QUMOND, constant kappa
# --------------------------------------------------------------------------
def nu(y):
    return 0.5 + np.sqrt(0.25 + 1.0 / y)


def toy_galaxy(eta, M=1e11 * MSUN, a=3.0 * KPC):
    r = np.logspace(np.log10(0.05 * a), np.log10(300 * a), 4000)
    gN = G * M / (r + a)**2
    gx = (nu(gN / A0) - 1.0) * gN          # extra dynamical acceleration
    kap = (1 - eta) / (1 + eta)
    gA = gx / (1 + kap)                    # the metric (AeST) part
    # matter-frame potential gradients
    dPhi = gN + (1 + kap) * gA             # = gN + gx
    dPsi = gN + (1 - kap) * gA
    MPhi = r**2 * dPhi / G
    MPsi = r**2 * dPsi / G
    p_r = (dPhi - dPsi) / (4 * np.pi * G * r) * C**2       # Pa (J/m^3)
    rho = np.gradient(r**2 * dPsi, r) / (4 * np.pi * G * r**2)  # kg/m^3
    rho_N = M * a / (2 * np.pi * r * (r + a)**3)            # Hernquist baryons
    rho_x = rho - rho_N
    p_t = p_r + 0.5 * r * np.gradient(p_r, r)                 # conservation
    ident = (MPhi - MPsi) - 4 * np.pi * r**3 * p_r / C**2
    sel = {}
    for R_kpc in (1, 10, 30, 100, 300):
        i = np.argmin(abs(r - R_kpc * KPC))
        sel[f'{R_kpc} kpc'] = {
            'gN/a0': float(gN[i] / A0),
            'p_r / (rho_x c^2)': float(p_r[i] / (rho_x[i] * C**2)),
            'p_t / (rho_x c^2)': float(p_t[i] / (rho_x[i] * C**2)),
            'eta measured from masses':
                float((MPsi[i] - r[i]**2 * gN[i] / G)
                      / (MPhi[i] - r[i]**2 * gN[i] / G)),
        }
    return {'eta_input': eta, 'kappa': kap,
            'max |identity residual| / M_Phi': float(np.max(abs(ident) / MPhi)),
            'far-field expectation p_r/(rho_x c^2) = 1/eta - 1': 1 / eta - 1,
            'profile': sel}


# --------------------------------------------------------------------------
# 4. What the measured values require; solar-system size of the slip
# --------------------------------------------------------------------------
MEASURED = [
    ('three over-bent SLACS lenses', 0.82, 0.88),
    ('eleven X-COP clusters (published hydrostatic bias)', 1.27, 1.42),
    ('three under-bent SLACS lenses', 1.36, 1.54),
]


def requirements():
    rows = []
    for name, lo, hi in MEASURED:
        row = {'systems': name, 'eta': [lo, hi]}
        for key, f in [
            ('light sees the extra gravity at (1+eta)/2 of its pull on matter',
             lambda e: (1 + e) / 2),
            ('conformal share kappa = (1-eta)/(1+eta)',
             lambda e: (1 - e) / (1 + e)),
            ('radial stress p_r / (mean extra density c^2) = (1/eta-1)/3',
             lambda e: (1 / e - 1) / 3),
            ('far-field radial equation of state w_r = 1/eta - 1',
             lambda e: 1 / e - 1),
        ]:
            v = sorted([f(lo), f(hi)])
            row[key] = [round(v[0], 3), round(v[1], 3)]
        rows.append(row)
    return rows


def solar_system(eta_worst=1.54):
    """Slip's effect on the Cassini light-deflection gamma."""
    Msun = MSUN
    Rsun = 6.957e8
    b = 1.6 * Rsun                        # Cassini 2002 closest approach
    L = 1.5e12                            # ~10 AU path half-length
    gN = G * Msun / b**2
    # simple nu: extra acceleration -> a0 (a constant) at high g
    gx = (nu(gN / A0) - 1.0) * gN
    # transverse deflection integrals: Newtonian 2GM/b; constant-gx field
    # gives 2 gx b asinh(L/b)
    newt = 2 * G * Msun / b
    extra = 2 * gx * b * np.arcsinh(L / b)
    frac = extra / newt
    kap = (1 - eta_worst) / (1 + eta_worst)
    # Einstein frame is slip-free; matter frame differs by 2 kappa phi.
    # gamma - 1 ~ -2 kappa (extra/Newtonian) in the deflection
    return {'gN at 1.6 Rsun (m/s^2)': gN,
            'extra / Newtonian deflection': frac,
            '|gamma - 1| from slip (worst eta)': abs(2 * kap * frac),
            'Cassini bound on |gamma - 1|': 2.3e-5}


def cluster_kinematic_slip():
    """Our prediction for the one published direct cluster slip measurement.

    Pizzuti et al. 2016 (arXiv:1602.03385): MACS J1206.2-0847, galaxy
    kinematics vs strong+weak lensing, eta(r200) = 1.01 +0.31 -0.28 (68%) for
    the TOTAL potential.  Galaxy orbits, not hot gas, so no non-thermal
    pressure enters.  For the total potential, eta_tot ~ 1 + f (eta_x - 1)
    with f = 1 - M_b/M_dyn; f from the X-COP gas fractions (0.79-0.87).
    """
    obs, up, dn = 1.01, 0.31, 0.28
    preds = [1 + f * (ex - 1) for f in (0.79, 0.87) for ex in (1.27, 1.42)]
    lo, hi = min(preds), max(preds)
    return {'measured eta_total (MACS J1206)': [obs, up, dn],
            'predicted eta_total from the X-COP slip': [round(lo, 3), round(hi, 3)],
            'offset in sigma (upper error)':
                [round((lo - obs) / up, 2), round((hi - obs) / up, 2)],
            'AeST (eta_x = 1) prediction': 1.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output-dir', default='../run-v1')
    args = ap.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    res = {'identity': check_identity(),
           'coupling_slots': coupling_slots(),
           'toy_galaxy': [toy_galaxy(e) for e in (0.85, 1.34, 1.45)],
           'requirements': requirements(),
           'solar_system': solar_system(),
           'cluster_kinematic_slip': cluster_kinematic_slip()}
    (out / 'results.json').write_text(json.dumps(res, indent=2, default=str))

    print('1. Identity (radial Einstein equation, any metric theory)')
    for k, v in res['identity'].items():
        print(f'   {k}: {v}')
    print('\n2. Coupling slots')
    for k, v in res['coupling_slots'].items():
        print(f'   {k}: {v}')
    print('\n3. Toy galaxy (1e11 Msun Hernquist, QUMOND, constant kappa)')
    for t in res['toy_galaxy']:
        print(f"   eta={t['eta_input']:.2f} kappa={t['kappa']:+.3f} "
              f"identity residual {t['max |identity residual| / M_Phi']:.1e} "
              f"far-field w_r expected {t['far-field expectation p_r/(rho_x c^2) = 1/eta - 1']:+.3f}")
        for R, d in t['profile'].items():
            print(f"      {R:>8}: gN/a0={d['gN/a0']:8.3f}  "
                  f"p_r/rho={d['p_r / (rho_x c^2)']:+.3f}  "
                  f"p_t/rho={d['p_t / (rho_x c^2)']:+.3f}  "
                  f"eta={d['eta measured from masses']:.3f}")
    print('\n4. What each measured eta requires')
    for row in res['requirements']:
        print(f"   {row['systems']}  eta {row['eta']}")
        for k, v in row.items():
            if k not in ('systems', 'eta'):
                print(f'      {k}: {v}')
    print('\n5. Solar system')
    for k, v in res['solar_system'].items():
        print(f'   {k}: {v:.3g}')
    print('\n6. The published cluster slip measurement (galaxy kinematics)')
    for k, v in res['cluster_kinematic_slip'].items():
        print(f'   {k}: {v}')


if __name__ == '__main__':
    main()

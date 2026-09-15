"""RC-1: the radiation budget and spectrum on the source side (see protocol.md).

    python rc1.py [--canonical] [--output-dir DIR]

For the incident companion densities that CC-2 stage 2A requires, the companions' energy is compared with the
radiation present today, and three declared laws of taking that energy from the microwave background are scored
against FIRAS: (a) a frequency shift at fixed photon number, (b) removal of whole photons at fixed shape, and (c) the
Planck-preserving combination dN/dt = -3hN, dU/dt = -4hU. The FIRAS scores are diagnostic (diagonal errors, fitted
colour temperature and Galaxy coefficient, no calibration model), as in the thermal-conversion pass.
"""
import json
import math
import sys
import time
from pathlib import Path
import numpy as np
from scipy.constants import h as H_PL, k as K_B, c as C, hbar, sigma as SIGMA_SB
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.special import zeta

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
ROOT = RESULTS.parents[1]
sys.path.insert(0, str(RESULTS/'companion-extensions'))
import evidence_io  # noqa: E402

T0 = 2.72548                         # K, Fixsen (2009)
T_REF = 2.725                        # K, the blackbody behind the LAMBDA residual columns
UVMM = (45., 100., 170.)             # nW m^-2 sr^-1, Hauser & Dwek's historical UV-to-mm range
ALPHA0 = 2.488993286382367e-4        # Mpc^-1, the archived loss coefficient
U_KMS = 300.                         # km/s: <v^2> = U^2 for both of stage 2A's incident distributions
DIST_MPC = (100., 1000., 3000.)
T_STAR = 5800.                       # K, a starlight proxy for the thermalization control
MPC_M = 3.0856775814913673e22
KPC_M = MPC_M/1e3
MSUN_KG = 1.98847e30                 # as in CF-1 (cf1.py), which stage 2A's densities use
GYR_S = 3.15576e16
A_RAD = 4*SIGMA_SB/C
FIRAS = ROOT/'time_revision/data/firas.txt'
W3_REFERENCE = {1.: 45.02, .9999: 86.00, .999: 4339.94}   # thermal-conversion pass, derivation.md section 5


def log(msg):
    print(f'[{time.time() - T_START:7.1f}s] {msg}', flush=True)


def planck_mjy(nu, T):
    """B_nu in MJy/sr."""
    return 2*H_PL*nu**3/C**2/np.expm1(H_PL*nu/(K_B*T))/1e-20


class Firas:
    """The thermal-conversion pass's diagnostic: residuals in kJy/sr against 2.725 K, diagonal errors, a fitted
    colour temperature and a signed coefficient on the Galaxy template (solved linearly)."""

    def __init__(self, path):
        d = np.loadtxt(path)
        assert d.ndim == 2 and d.shape[1] == 5 and np.all(d[:, 3] > 0)
        self.nu, self.obs, self.sig, self.gal = d[:, 0]*100*C, d[:, 2], d[:, 3], d[:, 4]
        self.ref = planck_mjy(self.nu, T_REF)
        self.wg = self.gal/self.sig

    def model(self, T, A):
        return 1000*(A*planck_mjy(self.nu, T) - self.ref)

    def profile(self, T, A):
        d = (self.obs - self.model(T, A))/self.sig
        g = float(self.wg @ d/(self.wg @ self.wg))
        r = d - g*self.wg
        return float(r @ r)

    def chi2(self, A):
        """chi^2 for a fixed amplitude A, profiled over the colour temperature and the Galaxy coefficient; the
        bounded scalar fit is cross-checked with a joint nonlinear fit."""
        f = minimize_scalar(lambda T: self.profile(T, A), bounds=(.5, 4.), method='bounded', options={'xatol': 1e-12})
        assert f.success and .5001 < f.x < 3.9999
        rep = least_squares(lambda x: (self.model(x[0], A) + x[1]*self.gal - self.obs)/self.sig, [f.x, 0.],
                            bounds=([.5, -np.inf], [4., np.inf]), x_scale='jac', ftol=1e-12, xtol=1e-12, gtol=1e-12)
        c2 = float(rep.fun @ rep.fun)
        assert rep.success and abs(c2 - f.fun) < 1e-5*max(1., f.fun)
        return c2, float(rep.x[0])


def planck_moments(T):
    """Photon number and energy densities of a Planck spectrum by numerical integration (W1)."""
    occ = lambda x: 0. if x > 700 else 1/math.expm1(x)
    I2 = quad(lambda x: x*x*occ(x), 0, np.inf, epsabs=0, epsrel=1e-13, limit=200)[0]
    I3 = quad(lambda x: x**3*occ(x), 0, np.inf, epsabs=0, epsrel=1e-13, limit=200)[0]
    s = K_B*T/(hbar*C)
    return s**3*I2/math.pi**2, hbar*C*s**4*I3/math.pi**2


def law_a_moments(q):
    """W2: photon number and energy after a frequency shift keeping a fraction q of each photon's energy."""
    occ = lambda x: 0. if x > 700 else 1/math.expm1(x)
    n0 = quad(lambda x: x*x*occ(x), 0, np.inf, epsabs=0, epsrel=1e-13)[0]
    u0 = quad(lambda x: x**3*occ(x), 0, np.inf, epsabs=0, epsrel=1e-13)[0]
    n1 = quad(lambda x: q**-3*x*x*occ(x/q), 0, np.inf, epsabs=0, epsrel=1e-13)[0]
    u1 = quad(lambda x: q**-3*x**3*occ(x/q), 0, np.inf, epsabs=0, epsrel=1e-13)[0]
    return n1/n0, u1/u0


def law_c_transport(T, ht, nu):
    """Law (c) along characteristics: frequency nu at time t came from nu*e^{ht}; number density per unit frequency
    carries the Jacobian e^{ht} and the removal factor e^{-3ht}. Returns the spectral energy density ratio to a Planck
    spectrum at T e^{-ht} (W4)."""
    x0 = H_PL*nu*math.exp(ht)/(K_B*T)
    n_t = math.exp(ht)*math.exp(-3*ht)*(nu*math.exp(ht))**2/np.expm1(x0)
    n_planck = nu**2/np.expm1(H_PL*nu/(K_B*T*math.exp(-ht)))
    return n_t/n_planck


def energy_density(rho_msun_kpc3):
    """Energy per unit volume, J/m^3, carried by companions at mass density rho: rest energy and the kinetic energy of
    stage 2A's incident distributions (<v^2> = U^2 for both S and M)."""
    rho = rho_msun_kpc3*MSUN_KG/KPC_M**3
    return rho*C**2*(1 + (U_KMS*1e3)**2/(2*C**2))


def boundary(fir, law, c2_ref):
    """The largest energy fraction of the microwave background a law can transfer before the diagnostic chi^2 rises
    by 4 above no transfer. Law (a): amplitude q^-3 with q = 1 - f. Law (b): amplitude s = 1 - f."""
    amp = (lambda f: (1 - f)**-3) if law == 'a' else (lambda f: 1 - f)
    g = lambda f: fir.chi2(amp(f))[0] - c2_ref - 4.
    lo, hi = 0., 1e-6
    while g(hi) < 0:
        lo, hi = hi, hi*2
        assert hi < .5
    return brentq(g, lo, hi, xtol=1e-14, rtol=1e-12)


def main():
    args = evidence_io.parse(__doc__)
    np.seterr(over='raise', invalid='raise', divide='raise')
    out = dict(scope='RC-1: the radiation budget and spectrum on the source side, for CC-2 stage 2A\'s requirement',
               inputs=dict(T0_K=T0, UVmm_nW_m2_sr=UVMM, alpha0_per_Mpc=ALPHA0, U_kms=U_KMS, firas=str(FIRAS.relative_to(ROOT)),
                           stage_2A=str((RESULTS/'companion-formation/formation-results.json').relative_to(ROOT))),
               validation={})
    # ---- validation
    n_num, u_num = planck_moments(T0)
    n_an = 2*zeta(3)/math.pi**2*(K_B*T0/(hbar*C))**3
    u_an = A_RAD*T0**4
    w1 = max(abs(n_num/n_an - 1), abs(u_num/u_an - 1))
    out['validation']['W1'] = dict(max_relative=w1, passed=bool(w1 < 1e-10))
    w2 = max(max(abs(nr - 1), abs(ur/q - 1)) for q in (.5, .9, .99) for nr, ur in [law_a_moments(q)])
    out['validation']['W2'] = dict(max_relative=w2, passed=bool(w2 < 1e-10))
    fir = Firas(FIRAS)
    w3 = {f'{q:g}': fir.chi2(q**-3)[0] for q in W3_REFERENCE}
    out['validation']['W3'] = dict(chi2=w3, reference=W3_REFERENCE,
                                   passed=bool(all(abs(w3[f'{q:g}'] - v) < .01 for q, v in W3_REFERENCE.items())))
    # on a grid scaled to each final temperature (x = h nu/kT_t from 0.01 to 30), so no exponent overflows
    w4_shape = max(abs(law_c_transport(T0, ht, x*K_B*T0*math.exp(-ht)/H_PL) - 1)
                   for ht in (1., 5., 10.) for x in np.geomspace(.01, 30., 50))
    split = quad(lambda t: 3*math.exp(-4*t), 0, 10, epsabs=0, epsrel=1e-13)[0]/quad(lambda t: math.exp(-4*t), 0, 10, epsabs=0, epsrel=1e-13)[0]
    out['validation']['W4'] = dict(shape_max_relative=w4_shape, energy_split=split,
                                   passed=bool(w4_shape < 1e-9 and abs(split/3 - 1) < 1e-9))
    log('validation: ' + ', '.join(f"{k}={v['passed']}" for k, v in out['validation'].items()))
    # ---- radiation today
    u_cmb = u_an
    u_uv = {f'{i:g}': 4*math.pi*i*1e-9/C for i in UVMM}
    out['radiation'] = dict(u_CMB_J_m3=u_cmb, u_UVmm_J_m3=u_uv)
    # ---- laws (a) and (b) against FIRAS
    c2_ref = fir.chi2(1.)[0]
    fa, fb = boundary(fir, 'a', c2_ref), boundary(fir, 'b', c2_ref)
    out['firas'] = dict(chi2_no_transfer=c2_ref, law_a_max_energy_fraction=fa, law_b_max_energy_fraction=fb,
                        law_a_colour_temperature_at_bound=fir.chi2((1 - fa)**-3)[1],
                        law_b_colour_temperature_at_bound=fir.chi2(1 - fb)[1])
    log(f'FIRAS: chi2(no transfer)={c2_ref:.2f}; largest transferable energy fraction (a) {fa:.3e}, (b) {fb:.3e}')
    # ---- law (c): cost to every beam, and the selectivity it would need
    h_per_gyr = ALPHA0*C/MPC_M*GYR_S
    out['law_c'] = dict(rate_h_per_Gyr=h_per_gyr, temperature_efold_Gyr=1/h_per_gyr,
                        extra_dimming_mag={f'{D:g} Mpc': 2.5/math.log(10)*3*ALPHA0*D for D in DIST_MPC},
                        photon_survival={f'{D:g} Mpc': math.exp(-3*ALPHA0*D) for D in DIST_MPC},
                        optical_to_microwave_rate_ratio_for_1pct_at_1Gpc=-math.log(.99)/(3*ALPHA0*1000.))
    # ---- controls from the plan: a thermal input (preservation) and a nonthermal input (thermalization). Each law maps
    # a diluted blackbody W*B(T) to: (a) W q^-3 B(qT); (b) s W B(T); (c) W B(T e^{-ht}), the dilution unchanged.
    def after(law, T, W):
        if law == 'a':
            return lambda x: W*.99**-3*planck_mjy(x, .99*T)
        if law == 'b':
            return lambda x: .99*W*planck_mjy(x, T)
        return lambda x: W*planck_mjy(x, T*math.exp(-.01))

    def deviation(spec, grid, tb):
        """Largest relative deviation from the best-fitting unit-amplitude Planck spectrum on the grid, with the fitted
        temperature bounded to tb (bounds that keep every exponent finite on that grid)."""
        f = minimize_scalar(lambda T: float(np.sum(np.log(spec(grid)/planck_mjy(grid, T))**2)),
                            bounds=tb, method='bounded', options={'xatol': 1e-12})
        return float(np.max(np.abs(spec(grid)/planck_mjy(grid, f.x) - 1)))
    grids = {'thermal': (np.geomspace(1e9, 3e12, 400), (.5, 20.)),
             'starlight proxy': (np.geomspace(1e13, 3e15, 400), (1e3, 1e5))}
    ctl = {}
    for name, T, W in (('thermal', T0, 1.), ('starlight proxy', T_STAR, 1e-13)):
        g_, tb = grids[name]
        ctl[name] = dict(input=deviation(lambda x: W*planck_mjy(x, T), g_, tb),
                         **{f'law {l}': deviation(after(l, T, W), g_, tb) for l in 'abc'})
    out['controls'] = dict(deviation_from_unit_planck=ctl, transfer=dict(law_a_q=.99, law_b_s=.99, law_c_ht=.01),
                           starlight_dilution=1e-13, starlight_temperature_K=T_STAR)
    # ---- stage 2A's requirement
    fr = json.loads((RESULTS/'companion-formation/formation-results.json').read_text(encoding='utf-8'))
    rows, labels = {}, {}
    u_all = u_cmb + u_uv['170']
    for name, cmb in fr['combos'].items():
        lab = cmb.get('labels', {})
        row = dict(system=cmb['system'], bath=cmb['bath'], sigma_over_m_cm2g=cmb['sigma_over_m_cm2g'],
                   bath_gravity=cmb['bath_gravity'], rho_B1=cmb.get('rho_B1'), rho_B2=cmb.get('rho_B2'),
                   in_regime=bool(cmb.get('rho_B1') and lab and not lab.get('outside_model_regime', True)))
        for bench in ('B1', 'B2'):
            rho = cmb.get(f'rho_{bench}')
            if not rho:
                continue
            E = energy_density(rho)
            Ti = T0*(1 + E/u_cmb)**.25
            row[bench] = dict(E_req_J_m3=E, over_CMB=E/u_cmb, over_CMB_plus_UVmm170=E/u_all,
                              law_a_supplies_fraction=fa*u_cmb/E, law_b_supplies_fraction=fb*u_cmb/E,
                              law_c_past_temperature_K=Ti, law_c_photons_removed=1 - (T0/Ti)**3,
                              law_c_duration_Gyr_at_alpha0=math.log(Ti/T0)/h_per_gyr)
        rows[name] = row
        if row['in_regime']:
            B1 = row['B1']
            labels[name] = dict(present_radiation_can_hold_requirement=bool(u_all >= B1['E_req_J_m3']),
                                firas_permits_law_a=bool(B1['law_a_supplies_fraction'] >= 1),
                                firas_permits_law_b=bool(B1['law_b_supplies_fraction'] >= 1))
    out['requirement'] = rows
    out['labels'] = dict(per_combination=labels,
                         law_c_needs_frequency_selectivity=bool(out['law_c']['extra_dimming_mag']['1000 Mpc'] > 2.5/math.log(10)*-math.log(.99)))
    ins = [r['B1'] for r in rows.values() if r['in_regime']]
    if ins:
        out['summary'] = dict(in_regime_combinations=len(ins),
                              B1_over_CMB_range=[min(b['over_CMB'] for b in ins), max(b['over_CMB'] for b in ins)],
                              B1_law_c_past_temperature_range_K=[min(b['law_c_past_temperature_K'] for b in ins),
                                                                 max(b['law_c_past_temperature_K'] for b in ins)])
    out['all_validation_passed'] = all(v['passed'] for v in out['validation'].values())
    text = json.dumps(out, indent=1, default=float)
    status = evidence_io.finish(args, 'radiation-budget', text, HERE/'rc1-results.json', ignore={'/runtime_seconds'})
    log(f"done; validation {out['all_validation_passed']}")
    return status if out['all_validation_passed'] else 1


T_START = time.time()
if __name__ == '__main__':
    sys.exit(main())

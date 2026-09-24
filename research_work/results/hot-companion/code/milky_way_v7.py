#!/usr/bin/env python3
"""Round 7: the Milky Way under our law, against Gaia and the classic Galactic measurements.

    python milky_way_v7.py --output-dir ../run-milky-way-v7

Nothing here is fitted to the Galaxy: the law's three constants come from SPARC and X-COP
(run-v3), and the Galaxy's ordinary matter comes from published models:
* M17: McMillan (2017) best fit (bulge, thin and thick stellar disks, HI and H2 disks) plus the
  stellar halo (1.4 x 10^9 Msun, Deason et al. 2019; broken power law of Deason et al. 2011);
* local census: the same shapes rescaled to the solar-neighbourhood census of McKee et al. (2015):
  stars and remnants 33.4, gas 13.7 Msun/pc^2 at the Sun;
* BR13: the short, 'maximal' stellar disk of Bovy & Rix (2013): R_d = 2.15 kpc, 38 Msun/pc^2 at
  8 kpc (gas and bulge as M17);
* M17 + corona: M17 plus a hot gas corona of 1.3 x 10^10 Msun inside 250 kpc (beta model,
  n = 0.0135 r^-1.5 cm^-3, Miller & Bregman 2015; gas, so it adds mass only).

Measurements compared (data/mw_literature_v7.json, verified against the papers):
* circular speed 5-27 kpc: Eilers et al. 2019, Zhou et al. 2023, Ou et al. 2024, Jiao et al. 2023;
* vertical pull 1.1 kpc above the Sun: Bovy & Rix 2013 (68 +- 4 Msun/pc^2), Holmberg & Flynn 2004
  (74 +- 6);
* enclosed mass at 20-200 kpc (Posti & Helmi 2019, Watkins et al. 2019, Vasiliev 2019, Deason
  et al. 2021, Correa Magnus & Vasiliev 2022);
* escape speed at the Sun (Piffl 2014, Monari 2018, Deason 2019, Koppelman & Helmi 2021, Necib &
  Lin 2022, Prudil 2022, Roche 2024);
* the inner Galaxy: baryons supply 0.88 +- 0.07 of the peak circular speed (Wegg, Gerhard &
  Portail 2016, from bulge microlensing and star counts).

Laws: ours (field-equation form, lap Phi = -div h), MOND (QUMOND, simple function, a0 = 1.2e-10)
and Newton with the same baryons; McMillan's dark halo is shown as the dark-matter reference.
The solver is in mw_model.py.
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import mw_model as M
import run as RUN
import law as LAW

G, K_SI, PC2, PC3 = M.G, M.K_SI, M.PC2, M.PC3
R0 = 8.2
REACH = dict(young=197.41 * 1.0227 * 10.0, old=197.41 * 1.0227 * 13.0)      # u t for t = 10 and 13 Gyr, kpc


def components():
    bulge = M.Spheroid('bulge', 98.4 * PC3, 0.075, 2.1, 1.8, 0.5)
    thin = M.Disk('thin', 895.9 * PC2, 2.50, 0.300)
    thick = M.Disk('thick', 183.0 * PC2, 3.02, 0.900)
    HI = M.Disk('HI', 53.1 * PC2, 7.0, 0.085, vert='sech2', Rm=4.0, kind='gas')
    H2 = M.Disk('H2', 2179.5 * PC2, 1.5, 0.045, vert='sech2', Rm=12.0, kind='gas')
    sh = lambda r, rb=27.0: np.where(r < rb, (np.maximum(r, 1.0) / rb) ** -2.3, (r / rb) ** -4.6)
    norm = 1.4e9 / M.Sphere('tmp', sh, 300.0).mass()
    # stellar halo: sigma_r 141, sigma_theta 75, sigma_phi 85 km/s (inner halo; Bond et al. 2010)
    halo = M.Sphere('stellar halo', lambda r: norm * sh(r), 300.0, sigma=float(np.sqrt((141 ** 2 + 75 ** 2 + 85 ** 2) / 3)))
    br = M.Disk('BR13 disk', 38 * np.exp(8 / 2.15) * PC2, 2.15, 0.300)
    mp, kpc3, msun = 1.6726e-24, (3.0857e21) ** 3, 1.989e33
    corona = M.Sphere('corona', lambda r: 1.14 * mp * 1.35e-2 * np.maximum(r, 1.0) ** -1.5 * kpc3 / msun, 250.0, kind='gas')
    return dict(bulge=bulge, thin=thin, thick=thick, HI=HI, H2=H2, halo=halo, br=br, corona=corona)


def models(C):
    star0 = (C['thin'].Sigma(R0) + C['thick'].Sigma(R0)) / PC2
    gas0 = (C['HI'].Sigma(R0) + C['H2'].Sigma(R0)) / PC2
    return {
        'M17': dict(bulge=1, thin=1, thick=1, HI=1, H2=1, halo=1),
        'local census': dict(bulge=1, thin=33.4 / star0, thick=33.4 / star0, HI=13.7 / gas0, H2=13.7 / gas0, halo=1),
        'BR13': dict(bulge=1, br=1, HI=1, H2=1, halo=1),
        'M17 + corona': dict(bulge=1, thin=1, thick=1, HI=1, H2=1, halo=1, corona=1),
    }


class Scaled:
    """A component with its mass multiplied by f (for the heat sums)."""
    def __init__(self, c, f):
        self.c, self.f, self.kind, self.sigma = c, f, c.kind, c.sigma

    def m_profile(self, s):
        return self.f * self.c.m_profile(s)


def nfw_mcmillan(R, z):
    """McMillan (2017) dark halo, for reference: rho0 = 0.00854 Msun/pc^3, r_h = 19.6 kpc."""
    rho0, rh = 0.00854 * PC3, 19.6
    r = np.sqrt(R ** 2 + z ** 2) + 1e-9
    m = 4 * np.pi * rho0 * rh ** 3 * (np.log(1 + r / rh) - (r / rh) / (1 + r / rh))
    g = -G * m / r ** 2
    return g * R / r, g * z / r


class Evaluator:
    def __init__(self, C, grid, fields, consts):
        self.C, self.grid, self.F, self.consts = C, grid, fields, consts
        self.eR = lambda x: grid.Re[np.argmin(abs(grid.Re - x))]
        self.ez = lambda x: grid.ze[np.argmin(abs(grid.ze - x))]
        self.Rt = np.unique([self.eR(x) for x in np.r_[np.arange(1.0, 4.0, 0.5), np.arange(4.0, 30.01, 0.5), [32, 35, 40, 50, 60, 80, 100]]])
        self.Rfar = np.unique([self.eR(x) for x in np.geomspace(R0, 2950, 70)])
        self.R0e = self.eR(R0)

    def run(self, parts, law, reach=REACH['old'], gd_scale=1.0):
        grid, C, F = self.grid, self.C, self.F
        gR = sum(f * F[k][0] for k, f in parts.items()); gz = sum(f * F[k][1] for k, f in parts.items())
        rho_b = sum(f * C[k].rho(grid.RR, grid.ZZ) for k, f in parts.items())
        jz = grid.nzh
        vb = np.sqrt(np.maximum(-parts.get('bulge', 0) * F['bulge'][0][:, jz] * grid.R, 0))
        C['bulge'].sigma = 0.65 * vb[grid.R < 10].max()                   # the SPARC rule for bulges
        hot = [Scaled(C[k], f) for k, f in parts.items() if k in ('bulge', 'halo')]
        consts = dict(self.consts, lam=self.consts['lam'] * gd_scale)
        if law == 'ours':
            S, hR, hz = M.heat_fields(hot, grid, consts['u_kms'])
        else:
            S = hR = hz = np.zeros_like(gR)
        eR, ez = M.law_extra(gR, gz, S, hR, hz, consts, law=law if law != 'lcdm' else 'newton', reach=reach if law == 'ours' else None, grid=grid)
        rho_ph = -M.divergence(grid, eR, ez) / (4 * np.pi * G) if law in ('ours', 'mond') else np.zeros_like(rho_b)
        rho = rho_b + rho_ph
        g1, _ = M.ring_sum(self.Rt, 0 * self.Rt, grid, rho)
        if law == 'lcdm':
            g1 = g1 + nfw_mcmillan(self.Rt, 0 * self.Rt)[0]
        v = np.sqrt(-g1 * self.Rt)
        # algebraic h in the plane (the SPARC convention), for comparison
        v_alg = np.sqrt(np.interp(self.Rt, grid.R, -(gR[:, jz] + eR[:, jz])) * self.Rt)
        # vertical pull above the Sun
        zs = [0.3, 0.5, 0.7, 1.1, 1.5, 2.0]
        Sig = []
        for zz in zs:
            _, gzz = M.ring_sum(np.array([self.R0e]), np.array([self.ez(zz)]), grid, rho)
            if law == 'lcdm':
                gzz = gzz + nfw_mcmillan(np.array([self.R0e]), np.array([self.ez(zz)]))[1]
            Sig.append(float(-gzz[0] / (2 * np.pi * G) / PC2))
        # enclosed 'as if' mass by Gauss's law on spheres (flux of the pull = flux of h)
        fR = RGI((grid.R, grid.z), gR + eR, bounds_error=False, fill_value=None)
        fz = RGI((grid.R, grid.z), gz + ez, bounds_error=False, fill_value=None)
        mu, w = np.polynomial.legendre.leggauss(600)
        Meff = {}
        for r in (20.0, 21.1, 50.0, 100.0, 200.0):
            st = np.sqrt(1 - mu ** 2); pts = np.c_[r * st, r * mu]
            hr = fR(pts) * st + fz(pts) * mu
            m = r ** 2 * 0.5 * np.sum(-hr * w) / G
            if law == 'lcdm':
                m += -nfw_mcmillan(np.array([r]), np.array([0.0]))[0][0] * r ** 2 / G
            Meff[str(r)] = float(m)
        # escape speed: integrate the pull outward along the plane
        gf, _ = M.ring_sum(self.Rfar, 0 * self.Rfar, grid, rho)
        if law == 'lcdm':
            gf = gf + nfw_mcmillan(self.Rfar, 0 * self.Rfar)[0]
        phi = np.concatenate([[0], np.cumsum(0.5 * (gf[1:] + gf[:-1]) * np.diff(self.Rfar))])     # Phi(R0) - Phi(R) (negative)
        Mb = float(sum(f * C[k].mass() for k, f in parts.items()))
        vesc = {str(x): float(np.sqrt(-2 * np.interp(x, self.Rfar, phi))) for x in (200.0, 400.0, 540.0, 1000.0)}
        if law == 'ours':
            vesc['reach'] = float(np.sqrt(-2 * np.interp(reach, self.Rfar, phi)))
            vesc['infinity'] = float(np.sqrt(-2 * phi[-1] + 2 * G * Mb / self.Rfar[-1]))
        if law == 'newton':
            vesc['infinity'] = float(np.sqrt(-2 * phi[-1] + 2 * G * Mb / self.Rfar[-1]))
        # the Sun: local pulls for the Cassini and wide-binary calculations
        iR = np.argmin(abs(grid.R - R0)); sun = dict(g_N_SI=float(np.hypot(gR[iR, jz], gz[iR, jz]) * K_SI), S_SI=float(S[iR, jz] * K_SI),
                                                    g_hot_SI=float(np.hypot(hR[iR, jz], hz[iR, jz]) * K_SI), g_total_SI=float(np.interp(R0, self.Rt, v) ** 2 / R0 * K_SI))
        return dict(R=self.Rt.tolist(), v=v.tolist(), v_algebraic=v_alg.tolist(), Sigma_z=dict(zip(map(str, zs), Sig)), M_eff=Meff, v_esc=vesc,
                    M_baryons=Mb, bulge_sigma=float(C['bulge'].sigma), sun=sun)


def compare_rc(res, rc):
    out = {}
    for name, d in rc.items():
        R, v, e = np.array(d['R']), np.array(d['v']), np.array(d['err'])
        e_tot = np.hypot(e, d['sys_frac'] * v)
        m = (R >= 5.0) & (R <= 27.5)
        p = np.interp(R[m], res['R'], res['v'])
        out[name] = dict(rms_kms=float(np.sqrt(np.mean((p - v[m]) ** 2))), mean_offset_kms=float(np.mean(p - v[m])),
                         chi2=float(np.sum(((p - v[m]) / e_tot[m]) ** 2)), n=int(m.sum()),
                         inner_5_10=float(np.mean((p - v[m])[(R[m] <= 10)])) if np.any(R[m] <= 10) else None,
                         outer_15_27=float(np.mean((p - v[m])[(R[m] >= 15)])) if np.any(R[m] >= 15) else None)
    return out


def sparc_transition(consts):
    """Median SPARC residual log10(g_obs/g_pred) by Newtonian pull, for our law and MOND."""
    gals = RUN.load_sparc()
    a, lam, u = consts['a_code'], consts['lam'], consts['u_kms']
    rows = []
    for g in gals:
        gobs = g['v'] ** 2 / g['r']
        rows.append(np.c_[g['gN'] * K_SI, gobs * K_SI, RUN.galaxy_g(g, a, u, lam) * K_SI, LAW.mond_simple(g['gN'], 1.2e-10 / K_SI) * K_SI])
    A = np.vstack(rows)
    out = []
    for lo in np.arange(-12.0, -8.99, 0.25):
        m = (np.log10(A[:, 0]) >= lo) & (np.log10(A[:, 0]) < lo + 0.25)
        if m.sum() >= 20:
            out.append(dict(log_gN=[float(lo), float(lo + 0.25)], n=int(m.sum()), ours=float(np.median(np.log10(A[m, 1] / A[m, 2]))),
                            mond=float(np.median(np.log10(A[m, 1] / A[m, 3])))))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args(); out = args.output_dir
    if out.exists(): raise FileExistsError('use a fresh output directory')
    out.mkdir(parents=True); t0 = time.monotonic()
    consts = json.loads((HERE.parent / 'run-v3/results.json').read_text())['constants']
    rc = json.loads((HERE.parent / 'data/mw_rotation_curves.json').read_text())
    C = components(); grid = M.Grid()
    F = {}
    for k, c in C.items():
        F[k] = M._one_field(c, grid, 60.0, 10, 40.0, 0.005)
        print(f'  Newtonian field of {k:12s} ({c.mass():.3e} Msun) done at {time.monotonic() - t0:.0f} s', flush=True)
    ev = Evaluator(C, grid, F, consts)
    MOD = models(C)
    results = {}
    for mname, parts in MOD.items():
        for law in ('ours', 'mond', 'newton') + (('lcdm',) if mname == 'M17' else ()):
            r = ev.run(parts, law)
            r['vs_rotation_curves'] = compare_rc(r, rc)
            results[f'{mname} | {law}'] = r
            vr = lambda x: np.interp(x, r['R'], r['v'])
            print(f"{mname:13s} {law:6s} v(3,5,8.2,12,16,20,25,27.5) = " + ' '.join(f'{vr(x):5.1f}' for x in (3, 5, R0, 12, 16, 20, 25, 27.5)) +
                  f" | rms vs Eilers/Zhou/Ou/Jiao " + '/'.join(f"{r['vs_rotation_curves'][k]['rms_kms']:.1f}" for k in ('Eilers2019', 'Zhou2023', 'Ou2024', 'Jiao2023')) +
                  f" | Sigma(1.1) {r['Sigma_z']['1.1']:.1f} | M(<20,50,100)/1e11 " + ' '.join(f"{r['M_eff'][k] / 1e11:.2f}" for k in ('20.0', '50.0', '100.0')) +
                  f" | v_esc(400,540) {r['v_esc']['400.0']:.0f} {r['v_esc']['540.0']:.0f}" + (f", reach {r['v_esc']['reach']:.0f}" if 'reach' in r['v_esc'] else ''), flush=True)
    # the inner Galaxy: share of the peak circular speed supplied by the ordinary matter (Wegg et al. 2016: 0.88 +- 0.07)
    for mname in MOD:
        vN = np.array(results[f'{mname} | newton']['v']); Rr = np.array(results[f'{mname} | newton']['R'])
        for law in ('ours', 'mond'):
            v = np.array(results[f'{mname} | {law}']['v']); m = (Rr >= 1.0) & (Rr <= 6.0)
            i = np.argmax(np.where(m, v, 0))
            results[f'{mname} | {law}']['inner_peak'] = dict(R=float(Rr[i]), v=float(v[i]), v_baryons=float(vN[i]), f_v=float(vN[i] / v[i]))
        print(f"{mname}: inner peak baryon share f_v = ours {results[f'{mname} | ours']['inner_peak']['f_v']:.2f} (R {results[f'{mname} | ours']['inner_peak']['R']:.1f}), "
              f"MOND {results[f'{mname} | mond']['inner_peak']['f_v']:.2f}; Wegg et al. 2016: 0.88 +- 0.07", flush=True)
    # what would have to change at the Sun: (a) the stellar disk our law needs; (b) a later switch-off
    need = {}
    for law in ('ours', 'mond'):
        scan = []
        for f in (0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5):
            parts = dict(bulge=1, thin=f, thick=f, HI=1, H2=1, halo=1)
            r = ev.run(parts, law)
            c = compare_rc(r, rc)
            scan.append(dict(disk_scale=f, Sigma_star_R0=float(f * 45.8), v_R0=float(np.interp(R0, r['R'], r['v'])), rms_Eilers=c['Eilers2019']['rms_kms'],
                             rms_all=float(np.mean([c[k]['rms_kms'] for k in c])), Sigma_1p1=r['Sigma_z']['1.1']))
        need[law] = scan
        print(f'{law}: disk scale / Sigma*(R0) / v(R0) / rms Eilers / Sigma(1.1): ' + '; '.join(f"{s['disk_scale']:.1f} {s['Sigma_star_R0']:.0f} {s['v_R0']:.0f} {s['rms_Eilers']:.1f} {s['Sigma_1p1']:.0f}" for s in scan), flush=True)
    gd_scan = []
    for gs in (1.0, 1.25, 1.5, 2.0, 3.0):
        r = ev.run(MOD['M17'], 'ours', gd_scale=gs)
        gd_scan.append(dict(g_d_scale=gs, g_d_SI=float(consts['g_d_SI'] * gs), v_R0=float(np.interp(R0, r['R'], r['v'])), v_20=float(np.interp(20, r['R'], r['v'])),
                            rms_Eilers=compare_rc(r, rc)['Eilers2019']['rms_kms'], Sigma_1p1=r['Sigma_z']['1.1']))
        print(f"  g_d x{gs}: v(R0) {gd_scan[-1]['v_R0']:.1f}, v(20) {gd_scan[-1]['v_20']:.1f}, rms Eilers {gd_scan[-1]['rms_Eilers']:.1f}, Sigma(1.1) {gd_scan[-1]['Sigma_1p1']:.1f}", flush=True)
    sparc = sparc_transition(consts)
    print('SPARC median residual log10(g_obs/g_pred) by log g_N: ' + '; '.join(f"{b['log_gN'][0]:.2f}: ours {b['ours']:+.3f}, MOND {b['mond']:+.3f}" for b in sparc), flush=True)
    # the Galaxy's mass profile for the dwarf-galaxy script (spherical shells of baryons and heat-weighted mass, M17)
    s = np.geomspace(1e-3, 3000.0, 1500)
    dm = np.zeros(s.size); dmk = np.zeros(s.size)
    for k, f in MOD['M17'].items():
        c = C[k]
        if isinstance(c, M.Disk):
            Rg = np.linspace(0, 200, 40001); Mc = np.concatenate([[0], np.cumsum(0.5 * np.diff(Rg) * (2 * np.pi * Rg * c.Sigma(Rg))[1:] + 0.5 * np.diff(Rg) * (2 * np.pi * Rg * c.Sigma(Rg))[:-1])])
            Ms = f * np.interp(s, Rg, Mc)
        else:
            Ms = f * c.m_profile(s)
        d = np.diff(np.concatenate([[0], Ms]))
        dm += d
        if k in ('bulge', 'halo'):
            dmk += LAW.k_from_sig2(c.sigma ** 2, consts['u_kms']) * d
    (out / 'milky_way_v7.json').write_text(json.dumps(dict(
        experiment='The Milky Way under our law vs Gaia rotation curves, vertical pull, enclosed mass, escape speed (round 7)',
        constants=consts, R0_kpc=R0, reach_kpc=REACH, models={k: v for k, v in MOD.items()},
        component_masses={k: float(c.mass()) for k, c in C.items()}, results=results, disk_needed=need, g_d_scan=gd_scan,
        sparc_transition=sparc, sun=results['M17 | ours']['sun'],
        galaxy_mass_profile=dict(s_kpc=s.tolist(), dm_Msun=dm.tolist(), dmk_Msun=dmk.tolist()),
        observed=dict(Sigma_1p1=dict(BovyRix2013=[68, 4], HolmbergFlynn2004=[74, 6]), wegg2016_f_v=[0.88, 0.07],
                      M_enclosed={'PostiHelmi2019 <20': [1.91e11, 0.18e11], 'Watkins2019 <21.1': [2.1e11, 0.35e11], 'Vasiliev2019 <50': [5.4e11, 0.95e11],
                                  'CorreaMagnusVasiliev2022 <50': [4.5e11, 0.4e11], 'Deason2021 <100': [6.07e11, 1.24e11], 'CorreaMagnusVasiliev2022 <100': [7.3e11, 0.85e11],
                                  'Vasiliev2019 <100': [8.5e11, 2.65e11], 'CorreaMagnusVasiliev2022 <200': [1.10e12, 0.25e12]},
                      v_esc={'Piffl2014 (3 R340)': [533, 48], 'Monari2018': [580, 63], 'Deason2019 (2 r200)': [528, 25], 'KoppelmanHelmi2021 (3 R340)': [497, 8],
                             'NecibLin2022 (2 r200)': [445, 17], 'Prudil2022 (2 r200)': [512, 66], 'Roche2024 (2 R200c)': [486, 8]}),
        seconds=time.monotonic() - t0), indent=1) + '\n')


if __name__ == '__main__':
    main()

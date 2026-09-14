"""Frozen baryon models for RPG-1 on the AQUAL grid.

SPARC galaxies use the capture-to-orbit receiver construction (stellar disk from SBdisk with
Upsilon_disk = 0.5 and an exponential tail, spherical bulge from 0.7 Vbul^2, exponential gas disk
of 1.33 M_HI fitted to Vgas) with a declared exponential vertical profile. The Milky Way uses the
archived variant I/II components. No observed rotation speed enters any model.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent if (HERE.parent/'capture-to-orbit').exists() else \
    Path('C:/Users/henry/Documents/Codex/photon-graviton/research_work/results')
sys.path.insert(0, str(RESULTS/'capture-to-orbit'))
sys.path.insert(0, str(HERE))
import inputs as I  # noqa: E402
import aqual as Q  # noqa: E402

A_STAR = 8.563335193921255e-11*Q.SI_TO_CODE          # archived simple-mu a0, (km/s)^2/kpc
DLNR, NTH = .035, 80


def sparc_components(d):
    """Surface densities (Msun/kpc^2) of the stellar and gas disks and the bulge's enclosed mass, per galaxy."""
    rows, cat = d['rotmod'], d['catalog']
    R = rows[:, 0]
    use = (R > 0) & np.isfinite(rows[:, 6]) & (rows[:, 6] > 0)
    Rs, lnS = R[use], np.log(I.UPSILON_DISK*rows[use, 6]*1e6)

    def sigma_star(x):
        x = np.asarray(x, float)
        out = np.interp(x, Rs, lnS)
        beyond = x > Rs[-1]
        out[beyond] = lnS[-1] - (x[beyond] - Rs[-1])/cat['rd']
        return np.exp(out)

    bul = (R > 0) & (rows[:, 5] > 0)
    if bul.any():
        mb = np.maximum.accumulate(I.UPSILON_BULGE*R[bul]*rows[bul, 5]**2/I.G)
        Rb = R[bul]

        def m_bulge(r):
            r = np.asarray(r, float)
            out = np.interp(np.log(r), np.log(Rb), mb)
            inner = r < Rb[0]
            out[inner] = mb[0]*(r[inner]/Rb[0])**3
            return out
    else:
        m_bulge = None
    Mgas = I.HELIUM*cat['MHI9']*1e9
    good = R > 0
    if Mgas > 0 and np.any(rows[good, 3] != 0):
        hg = I.fit_gas_scale(R[good], rows[good, 3], Mgas, cat['rd'])['h_kpc']

        def sigma_gas(x):
            return Mgas/(2*np.pi*hg*hg)*np.exp(-np.asarray(x, float)/hg)
    else:
        hg, sigma_gas = 0., None
    return dict(sigma_star=sigma_star, sigma_gas=sigma_gas, m_bulge=m_bulge, h_gas=hg, rd=cat['rd'])


def _cylinder_mass(sigma, r_out):
    x = np.geomspace(1e-6, r_out, 20001)
    f = 2*np.pi*x*sigma(x)
    return float(np.sum(.5*(f[1:] + f[:-1])*np.diff(x)))


def sparc_grid(d, comp, total_mass, refine=1, r_out_factor=1.):
    r_obs = d['r']
    r_m = np.sqrt(I.G*total_mass/A_STAR)
    r_in = .01*min(comp['rd'], r_obs[0])
    r_out = 1000*max(r_obs[-1], r_m, 10*comp['rd'], 5*comp['h_gas'])*r_out_factor
    nr = int(np.ceil(np.log(r_out/r_in)/DLNR))*refine
    return Q.Grid(r_in, r_out, nr, NTH*refine), r_m


def sparc_masses(grid, comp, hz):
    def sigma(x):
        out = comp['sigma_star'](x)
        if comp['sigma_gas'] is not None:
            out = out + comp['sigma_gas'](x)
        return out
    m, m_in = Q.exponential_disk_cell_masses(grid, sigma, hz)
    if comp['m_bulge'] is not None:
        mb, mb_in = Q.spherical_cell_masses(grid, comp['m_bulge'])
        m, m_in = m + mb, m_in + mb_in
    return m, m_in


def sparc_total_mass(comp):
    total = _cylinder_mass(comp['sigma_star'], 1e4*comp['rd'])
    if comp['sigma_gas'] is not None:
        total += _cylinder_mass(comp['sigma_gas'], 1e4*max(comp['rd'], comp['h_gas']))
    if comp['m_bulge'] is not None:
        total += float(comp['m_bulge'](np.array([1e6]))[0])
    return total


def milky_way_masses(grid, variant, hz_gas=.1):
    components, gas, _ = I._milky_way_definitions()
    m = np.zeros((grid.nr, grid.nth))
    m_in = 0.
    for M, a, b in components[variant]:
        if a == 0:
            mm, mi = Q.spherical_cell_masses(grid, lambda r, M=M, b=b: M*r**3/(r*r + b*b)**1.5)
        else:
            def rho(R, z, M=M, a=a, b=b):
                zb = np.sqrt(z*z + b*b)
                return b*b*M/(4*np.pi)*(a*R*R + (a + 3*zb)*(a + zb)**2)/((R*R + (a + zb)**2)**2.5*zb**3)
            mm, mi = Q.smooth_cell_masses(grid, rho)
        m, m_in = m + mm, m_in + mi

    def sigma_gas(x):
        x = np.maximum(np.asarray(x, float), 1e-9)
        return sum(S*1e6*np.exp(-hole/x - x/h) for S, h, hole in gas)
    mm, mi = Q.exponential_disk_cell_masses(grid, sigma_gas, hz_gas)
    return m + mm, m_in + mi, dict(stellar=float(sum(c[0] for c in components[variant])),
                                   gas=_cylinder_mass(sigma_gas, 1e4))

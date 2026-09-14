"""Frozen inputs for the capture-to-orbit test.

SPARC and Milky Way receivers, the exact-third reference reservoir and its
incident companion field, plus the archived original and MOND-guided
predictions used for side-by-side comparison. No observed rotation speed
enters any receiver model, potential or population.
"""
from pathlib import Path
import ast
import hashlib
import io
import json
import zipfile
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar
from scipy.special import gamma, gammainc, i0e, i1e, k0e, k1e

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
ROOT = HERE.parents[2]
OLD = RESULTS/'isotropic-galaxy-transfer'
EXT = RESULTS/'companion-extensions'
DATA = ROOT/'temporal_candidate_audit/data'
G = 4.30091727003628e-6
UPSILON_DISK, UPSILON_BULGE, HELIUM = .5, .7, 1.33     # archived baryon model; SPARC gas factor
CP = json.loads((OLD/'third-radiation-retention-results.json').read_text())['models']['attenuated']
assert CP['q'] == 1/3
GRID = np.geomspace(1e-4, 1e6, 6000)                  # model-defined radial domain, kpc

INPUT_FILES = [OLD/'third-radiation-retention-results.json', OLD/'model-comparison-predictions.json',
               DATA/'Rotmod_LTG.zip', DATA/'SPARC_Lelli2016c.mrt', EXT/'mond-inventory-results.json',
               EXT/'mond-cross-scale-results.json', OLD/'milky-way-current.py',
               RESULTS/'joint-galaxy-audit/milky-way-predictions.json', RESULTS/'milky-way-capture/inputs.json']


def input_hashes():
    return {str(p.relative_to(ROOT)).replace('\\', '/'): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in INPUT_FILES}


def reference_field(rd, L9, r=GRID, order=192):
    """Exact-third reference (paper equations 5-7) on r.

    Returns the attenuated mean intensity J, the net outward flux fraction
    zeta=int I mu dmu/int I dmu (negative: inward), the density, enclosed mass
    and the analytic total inventory.
    """
    a = rd*CP['scale_to_disk']
    x = r/a
    mu, w = leggauss(order)
    t = x[:, None]*mu
    B2 = 1 + x[:, None]**2*(1 - mu*mu)
    B = np.sqrt(B2)
    tau = CP['k0_per_kpc']*a*(t/(2*B2*(B2 + t*t)) + (np.arctan(t/B) + np.pi/2)/(2*B**3))
    I = np.exp(-np.maximum(tau, 0))
    J = I@w/2
    zeta = (I@(w*mu))/(I@w)
    X = L9/rd**2
    eta = X**(1/3)/(1 + X**(1/3))
    C = 2*CP['C_Msun_kpc3']*eta
    rho = C*J/(1 + x*x)**2
    mass = 4*np.pi*cumulative_trapezoid(r*r*rho, r, initial=0) + 4*np.pi*rho[0]*r[0]**3/3
    T = np.pi*CP['k0_per_kpc']*a/2
    total = C/CP['k0_per_kpc']*np.pi*a*a*(T**(2/3)*gamma(1/3)*gammainc(1/3, T) - 1 + np.exp(-T))
    return dict(J=J, zeta=zeta, rho=rho, mass=mass, total=float(total), eta=float(eta), a=a)


def _catalog():
    cat = {}
    for line in (DATA/'SPARC_Lelli2016c.mrt').read_text().splitlines():
        f = line.split()
        if len(f) == 19:
            try:
                cat[f[0]] = dict(rd=float(f[11]), L9=float(f[7]), MHI9=float(f[13]))
            except ValueError:
                pass
    return cat


def exponential_disk_v2(R, M, h):
    """Freeman thin exponential disk: v^2 = 2 G M/h y^2 [I0K0 - I1K1](y), y=R/2h."""
    y = np.maximum(R/(2*h), 1e-12)
    return 2*G*M/h*y*y*(i0e(y)*k0e(y) - i1e(y)*k1e(y))


def fit_gas_scale(R, vgas, M, rd):
    """Gas scale length reproducing SPARC's Vgas|Vgas| at fixed catalog gas mass.

    Ordinary-matter data reduction only; the observed rotation speed is not used.
    """
    target = vgas*np.abs(vgas)
    lo, hi = np.log(.1*rd), np.log(30*rd)

    def loss(lh):
        return float(np.mean((exponential_disk_v2(R, M, np.exp(lh)) - target)**2))

    grid = np.linspace(lo, hi, 241)
    i = int(np.argmin([loss(v) for v in grid]))
    best = minimize_scalar(loss, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, 240)]), method='bounded',
                           options={'xatol': 1e-10})
    h = float(np.exp(best.x))
    return dict(h_kpc=h, boundary=bool(i in (0, 240)), rms_v2=float(np.sqrt(best.fun)),
                target_rms_v2=float(np.sqrt(np.mean(target**2))))


def _disk_enclosed(sigma, r):
    """Mass of a thin disk inside cylindrical radius r (its monopole enclosed mass)."""
    return cumulative_trapezoid(2*np.pi*r*sigma, r, initial=0) + np.pi*r[0]**2*sigma[0]


def sparc_receivers(rows, cat, r=GRID):
    """Stellar disk (SBdisk), spherical bulge (Vbul) and exponential gas receivers.

    rows: rotmod columns Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul.
    """
    R = rows[:, 0]
    use = (R > 0) & np.isfinite(rows[:, 6]) & (rows[:, 6] > 0)
    Rs, S = R[use], UPSILON_DISK*rows[use, 6]*1e6                 # Msun/kpc^2
    lnS = np.interp(r, Rs, np.log(S))
    beyond = r > Rs[-1]
    lnS[beyond] = np.log(S[-1]) - (r[beyond] - Rs[-1])/cat['rd']
    sigma_star = np.exp(lnS)
    bul = (R > 0) & (rows[:, 5] > 0)
    if bul.any():
        mb = np.maximum.accumulate(UPSILON_BULGE*R[bul]*rows[bul, 5]**2/G)
        Rb = R[bul]
        mbul = np.interp(np.log(r), np.log(Rb), mb)
        mbul[r < Rb[0]] = mb[0]*(r[r < Rb[0]]/Rb[0])**3
    else:
        mbul = np.zeros_like(r)
    Mgas = HELIUM*cat['MHI9']*1e9
    good = R > 0
    if Mgas > 0 and np.any(rows[good, 3] != 0):
        fit = fit_gas_scale(R[good], rows[good, 3], Mgas, cat['rd'])
        sigma_gas = Mgas/(2*np.pi*fit['h_kpc']**2)*np.exp(-r/fit['h_kpc'])
    else:
        fit, sigma_gas = None, np.zeros_like(r)
    m_star, m_gas = _disk_enclosed(sigma_star, r), _disk_enclosed(sigma_gas, r)
    disk_dm = 2*np.pi*r*r*(sigma_star + sigma_gas)                 # dM/dln r
    bulge_dm = np.gradient(mbul, np.log(r))
    return dict(mass=m_star + m_gas + mbul, disk_dm_dlnr=disk_dm, bulge_dm_dlnr=np.maximum(bulge_dm, 0),
                stellar_disk_mass=float(m_star[-1]), gas_mass=float(m_gas[-1]), bulge_mass=float(mbul[-1]),
                gas_fit=fit)


def sparc_galaxies():
    """All 149 archived SPARC galaxies with receivers and frozen comparison rows."""
    cat = _catalog()
    guided = {v['galaxy']: v for v in json.loads((EXT/'mond-inventory-results.json').read_text())['rows']}
    out = []
    with zipfile.ZipFile(DATA/'Rotmod_LTG.zip') as archive:
        for name, g in guided.items():
            rows = np.atleast_2d(np.loadtxt(io.BytesIO(archive.read(name + '_rotmod.dat'))))
            c = cat[name]
            out.append(dict(name=name, split=g['split'], catalog=c, rotmod=rows, r=np.array(g['R_kpc']),
                            y=np.array(g['observed_kms']), reference=np.array(g['reference_kms']),
                            mond_guided=np.array(g['selected_kms']), mond_raw=np.array(g['MOND_raw_kms']),
                            total_inventory=g['total_inventory_Msun'], capped=g['capped_radii'] > 0,
                            reference_enclosed=np.array(g['reference_enclosed_Msun'])))
    assert len(out) == 149 and sum(len(d['r']) for d in out) == 3150
    return out


def _milky_way_definitions():
    """Archived Milky Way components and the reference capture function, parsed
    from milky-way-current.py without executing its driver."""
    source = OLD/'milky-way-current.py'
    tree = ast.parse(source.read_text())
    keep = [n for n in tree.body if isinstance(n, ast.Assign) and any(
        isinstance(t, ast.Name) and t.id in ('components', 'gas') for t in n.targets)]
    keep += [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'capture']
    ns = dict(np=np, cp=CP, leggauss=leggauss, cumulative_trapezoid=cumulative_trapezoid,
              gamma=gamma, gammainc=gammainc)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(source), 'exec'), ns)
    return ns['components'], ns['gas'], ns['capture']


def milky_way_receivers(variant, r=GRID, order=192):
    """Spherically averaged archived baryons (Miyamoto-Nagai/Plummer plus thin gas disks).

    Stellar enclosed masses follow from Gauss's law with the angle-averaged
    radial force, which stays smooth in angle even for thin disks. Thin gas disks
    enclose their mass inside cylindrical R<r. The a=0 Plummer component is the
    nonrotating bulge.
    """
    components, gas, _ = _milky_way_definitions()
    mu, w = leggauss(order)
    R, z = r[:, None]*np.sqrt(1 - mu*mu), r[:, None]*mu
    disk, bulge = np.zeros(r.shape), np.zeros(r.shape)
    for M, a, b in components[variant]:
        B = np.sqrt(z*z + b*b)
        d = R*R + (a + B)**2
        radial = (M*R*R/d**1.5 + M*(a + B)*z*z/(B*d**1.5))/r[:, None]      # r.g/(G r)
        enclosed = r*r*(radial@w/2)
        if a == 0:
            bulge += enclosed
        else:
            disk += enclosed
    sigma_gas = sum(S*1e6*np.exp(-hole/r - r/h) for S, h, hole in gas)
    m_gas = _disk_enclosed(sigma_gas, r)
    lr = np.log(r)
    return dict(mass=disk + bulge + m_gas,
                disk_dm_dlnr=np.maximum(np.gradient(disk, lr), 0) + 2*np.pi*r*r*sigma_gas,
                bulge_dm_dlnr=np.maximum(np.gradient(bulge, lr), 0), stellar_mass=float(disk[-1] + bulge[-1]),
                gas_mass=float(m_gas[-1]), stellar_mass_nominal=float(sum(c[0] for c in components[variant])))


def milky_way_runs():
    """The 18 archived Milky Way scenarios with frozen reference/MOND-guided rows."""
    components, _, _ = _milky_way_definitions()
    base = json.loads((RESULTS/'joint-galaxy-audit/milky-way-predictions.json').read_text())
    archived = json.loads((EXT/'mond-cross-scale-results.json').read_text())['runs']
    runs = []
    for run in archived:
        rows = run['rows']
        v = [b for b in base if b['baryons'] == run['baryons'] and b['model'] == 'baryons'
             and b['observable'] == 'vc_kms']
        assert len(v) == 38
        mstar = sum(c[0] for c in components[run['baryons']])
        runs.append(dict(baryons=run['baryons'], rd=run['Rd_kpc'], lf=run['luminosity_proxy_factor'],
                         L9=mstar/.5*run['luminosity_proxy_factor']/1e9, total_inventory=run['total_inventory_Msun'],
                         R=np.array([x['R_kpc'] for x in rows]), y=np.array([x['observed_kms'] for x in rows]),
                         vb=np.array([x['predicted'] for x in v]), split=np.array([x['split'] for x in v]),
                         reference=np.array([x['reference_kms'] for x in rows]),
                         mond_guided=np.array([x['redistributed_kms'] for x in rows]),
                         mond_raw=np.array([x['MOND_raw_kms'] for x in rows])))
    assert len(runs) == 18
    return runs

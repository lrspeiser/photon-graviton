"""Conservative 3D response via real spherical harmonics; no stellar fit."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import time
import argparse

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline
from scipy.special import roots_legendre, sph_legendre_p_all

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent / 'rotating-bar-orbits'))
sys.path.insert(0, str(HERE.parent / 'bar-field-foundation'))
from fast_multipole import FastMultipole

OLD = HERE.parent / 'conservative-field-completion/run.py'
spec = importlib.util.spec_from_file_location('axisymmetric_completion', OLD)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
CACHE = ROOT / 'research_work/data-cache/full-bar-completion'
CACHE.mkdir(parents=True, exist_ok=True)


def harmonics(ellmax, mu, phi):
    ell = np.array([l for l in range(0, ellmax+1, 2) for m in range(0, l+1, 2)])
    emm = np.array([m for l in range(0, ellmax+1, 2) for m in range(0, l+1, 2)])
    leg = sph_legendre_p_all(ellmax, ellmax, np.arccos(mu), diff_n=1)
    norm = np.where(emm > 0, np.sqrt(2), 1.)[:, None, None]
    co = np.cos(emm[:, None] * phi)
    si = np.sin(emm[:, None] * phi)
    y = (norm * leg[0, ell, emm, :, None] * co[:, None, :]).reshape(len(ell), -1)
    dt = (norm * leg[1, ell, emm, :, None] * co[:, None, :]).reshape(len(ell), -1)
    dp = (-norm * emm[:, None, None] * leg[0, ell, emm, :, None] * si[:, None, :]).reshape(len(ell), -1)
    return ell, emm, y, dt, dp


def integrate_coefficients(r, ell, qr, u):
    # Integration by parts of div(Q), with Q=-B. u=integral Q_t.grad_Omega(Y).
    inner, outer = old.scaled_integrals(r, -u-ell*qr, -u+(ell+1)*qr, ell)
    coeff = -(inner+outer)/(2*ell+1)
    # The uncut empirical power law has no finite monopole potential at infinity.
    # Fix an additive gauge, preserving exact angular-mean radial flux.
    coeff[:, 0] = cumulative_trapezoid(qr[:, 0], r, initial=0)
    return coeff


def disk_model():
    fingerprint = hashlib.sha256(OLD.read_bytes()).hexdigest()
    path = CACHE / f'axis-disk-{fingerprint[:12]}.npz'
    if path.exists():
        f = np.load(path)
        assert str(f['source_hash']) == fingerprint
        return old.Expansion(f['r'], f['ell'], f['phi'], f['grad']), path
    print('Building shared disk harmonic expansion', flush=True)
    d = old.Disk(False)
    np.savez_compressed(path, r=d.r, ell=d.ells, phi=d.phi(np.log(d.r)),
                        grad=d.grad(np.log(d.r)), source_hash=fingerprint)
    return d, path


def build(label, nr, nmu, nphi, ellmax, disk, rmax=80., split=True):
    started = time.monotonic()
    r = np.geomspace(.001, rmax, nr)
    mu, w = roots_legendre(nmu)
    phi = np.arange(nphi) * 2*np.pi/nphi
    ell, emm, Y, T, F = harmonics(64, mu, phi)
    pick = ell <= ellmax
    # Resolve both signs of every retained Fourier harmonic, including m=lmax.
    assert nphi > 2*ellmax
    le, me = ell[pick], emm[pick]
    weights = np.repeat(w, nphi) * 2*np.pi/nphi
    assert np.max(abs(np.sum(Y[pick]**2*weights,axis=1)-1)) < 1e-10
    by = (Y[pick] * weights).T
    bt = (T[pick] * weights).T
    bf = (F[pick] * (weights/np.repeat(np.sqrt(1-mu*mu), nphi))).T
    f = np.load(ROOT / 'research_work/data-cache/bar-field/bar-L64.npz')
    assert np.array_equal(f['l'], ell) and np.array_equal(f['m'], emm)
    spline = CubicSpline(np.log(f['r']), f['coeff'], axis=0)
    m0 = emm == 0
    nuclei = old.CachedAxisymmetric(ROOT / 'research_work/data-cache/bar-field/nuclei-L16.npz')
    qr = np.empty((nr, len(le)))
    uu = np.empty_like(qr)
    for start in range(0, nr, 8):
        rr = r[start:start+8]
        shapes = (len(rr), nmu)
        rad = np.broadcast_to(rr[:, None], shapes).ravel()
        mus = np.broadcast_to(mu, shapes).ravel()
        ax = np.array(disk.evaluate(rad, mus)) + np.array(nuclei.evaluate(rad, mus))
        ax[1] += old.G * 4.1e6 * rad / (rad*rad+.001**2)**1.5
        gr0 = np.repeat(ax[1].reshape(shapes), nphi, axis=1)
        gt0 = np.repeat(ax[2].reshape(shapes), nphi, axis=1)
        c = spline(np.log(rr))
        dc = spline(np.log(rr), 1) / rr[:, None]
        gr = gr0 + dc @ Y
        gt = gt0 + (c @ T) / rr[:, None]
        gp = (c @ F) / (rr[:, None] * np.repeat(np.sqrt(1-mu*mu), nphi))
        strength = np.sqrt(gr*gr+gt*gt+gp*gp)
        mult = old.A * (strength/old.ASTAR)**(old.P-1)
        # Subtract the full axisymmetric reference BEFORE projection. Its
        # thin-disk response is solved separately with many more angular modes.
        if split:
            ar = gr0 + dc[:,m0] @ Y[m0]
            at = gt0 + (c[:,m0] @ T[m0]) / rr[:,None]
            amult = old.A * (np.hypot(ar,at)/old.ASTAR)**(old.P-1)
            mult_gr, mult_gt = mult*gr-amult*ar, mult*gt-amult*at
        else:
            mult_gr, mult_gt = mult*gr, mult*gt
        qr[start:start+len(rr)] = mult_gr @ by
        uu[start:start+len(rr)] = mult_gt @ bt + (mult*gp) @ bf
        if start % 64 == 0:
            print(label, start, '/', nr, 'seconds', round(time.monotonic()-started, 1), flush=True)
    coeff = integrate_coefficients(r, le, qr, uu)
    assert np.isfinite(coeff).all()
    path = CACHE / f'{label}.npz'
    np.savez_compressed(path, r=r, l=le, m=me, coeff=coeff, qr=qr, u=uu,
                        nr=nr, nmu=nmu, nphi=nphi, ellmax=ellmax,
                        source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    return FastMultipole(r, coeff, le, me), path


class SumField:
    def __init__(self, *parts):
        self.parts = parts

    def evaluate(self, xyz):
        pairs = [part.evaluate(xyz) for part in self.parts]
        return sum(p[0] for p in pairs), sum(p[1] for p in pairs)


def load_field(path):
    f=np.load(path)
    return FastMultipole(f['r'],f['coeff'],f['l'],f['m'])


def axis_field(disk, rmax=80., ellmax=128, angular_nodes=512):
    baryons = old.Baryons.__new__(old.Baryons)
    baryons.disk = disk
    baryons.bar = old.CachedAxisymmetric(ROOT/'research_work/data-cache/bar-field/bar-L64.npz')
    baryons.nuclei = old.CachedAxisymmetric(ROOT/'research_work/data-cache/bar-field/nuclei-L16.npz')
    comp = old.Completion(baryons, True, rmax=rmax, ellmax=ellmax, angular_nodes=angular_nodes)
    coeff = comp.phi(np.log(comp.r)) * np.sqrt(4*np.pi/(2*comp.ells+1))
    field = FastMultipole(comp.r, coeff, comp.ells, np.zeros_like(comp.ells))
    path = CACHE / f'axisymmetric-reference-{rmax:g}-L{ellmax}-N{angular_nodes}.npz'
    np.savez_compressed(path,r=comp.r,coeff=coeff,l=comp.ells,m=np.zeros_like(comp.ells))
    return field, path


def analytical_control():
    # Known manufactured gradient Q=grad[f_l(r)Y_lm], with f_l=r^l exp(-r²/2).
    # The Poisson solution must recover its potential (up to monopole gauge).
    checks = []
    for nr in [600, 1200]:
        r = np.geomspace(1e-5, 30., nr)
        ell = np.array([0, 2, 4])
        f = r[:, None]**ell * np.exp(-r[:, None]**2/2)
        fp = f*(ell/r[:, None]-r[:, None])
        u = ell*(ell+1)*f/r[:, None]
        solved = integrate_coefficients(r, ell, fp, u)
        spline = CubicSpline(np.log(r), solved, axis=0)
        x = np.geomspace(.1, 6, 80)
        exact = x[:, None]**ell*np.exp(-x[:, None]**2/2)*(ell/x[:, None]-x[:, None])
        got = spline(np.log(x), 1)/x[:, None]
        err = float(np.max(abs(got-exact))/np.max(abs(exact)))
        checks.append(dict(radial_nodes=nr, radial_gradient_error_over_peak=err))
    assert checks[1]['radial_gradient_error_over_peak'] < .001
    assert checks[1]['radial_gradient_error_over_peak'] < checks[0]['radial_gradient_error_over_peak']*.4
    return checks


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--unsplit',action='store_true')
    args=parser.parse_args()
    control = analytical_control()
    disk, dp = disk_model()
    ap=None
    if not args.unsplit:
        print('Solving high-angular-order axisymmetric reference',flush=True)
        axis, ap = axis_field(disk)
    configs = [('coarse', 192, 64, 96, 24), ('fine', 384, 128, 128, 48)]
    if not args.unsplit:
        configs.append(('finer',768,256,192,64))
    models = []
    paths = []
    for config in configs:
        label,*shape=config
        model, path = build(('unsplit-' if args.unsplit else '')+label,*shape, disk=disk,split=not args.unsplit)
        models.append(model if args.unsplit else SumField(axis,model))
        paths.append(path)
    probes = np.array([[R*np.cos(phi),R*np.sin(phi),z] for R in [.5,1.,3.,8.,20.]
                       for z in [0.,.1,.8,1.5] for phi in [0.,np.pi/6,np.pi/2]])
    ac = models[-2].evaluate(probes)[1]
    af = models[-1].evaluate(probes)[1]
    relative = np.linalg.norm(af-ac, axis=1)/np.linalg.norm(af, axis=1)
    # Component accelerations are differentiated from the same potential.
    dx = 1e-4
    fd = np.column_stack([-(models[-1].evaluate(probes+np.eye(3)[j]*dx)[0]-
                            models[-1].evaluate(probes-np.eye(3)[j]*dx)[0])/(2*dx) for j in range(3)])
    derivative_error = float(np.max(np.linalg.norm(fd-af,axis=1)/np.linalg.norm(af,axis=1)))
    assert derivative_error < 1e-5
    inputs = [Path(__file__), OLD, dp] + ([ap] if ap else []) + paths
    inputs += [ROOT/f'research_work/data-cache/bar-field/{name}' for name in ['bar-L64.npz','nuclei-L16.npz']]
    results = dict(scope='Full-bar conservative potential prototype; no stellar fit or holdout score.',
                   formula_provenance='Known spherical-harmonic Poisson solution with integration by parts; frozen empirical response.',
                   parameters=old.PARAM, radial_domain_kpc=[.001,80.], configurations=configs,
                   numerical_split=not args.unsplit,
                   axisymmetric_angular_order=None if args.unsplit else 128,
                   axisymmetric_angular_nodes=None if args.unsplit else 512,
                   force_comparison_pair=[configs[-2][0],configs[-1][0]],
                   manufactured_gradient_controls=control, potential_derivative_error=derivative_error,
                   force_refinement_maximum=float(relative.max()), force_refinement_median=float(np.median(relative)),
                   provisional_extra_force_refinement_target=.01,
                   refinement_target_met=bool(relative.max()<.01),
                   rows=[dict(position_kpc=x.tolist(),coarse_extra_acceleration=a.tolist(),fine_extra_acceleration=b.tolist(),
                              refinement_fraction=float(e)) for x,a,b,e in zip(probes,ac,af,relative)],
                   input_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                   outer_boundary_refined=False, disk_field_refined=False,
                   full_orbit_inference_ready=False, holdouts_opened=False)
    name='verification-unsplit-results.json' if args.unsplit else 'results.json'
    (HERE/name).write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('Force refinement max/median',relative.max(),np.median(relative),flush=True)
    print('Potential derivative check',derivative_error,flush=True)


if __name__ == '__main__':
    main()

"""Reconstruct the source of the full empirical extra potential, without fitting data."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
from scipy.special import sph_legendre_p_all

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DRIVER = HERE.parent/'full-bar-orbits/run.py'
spec = importlib.util.spec_from_file_location('source_driver', DRIVER)
orbit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(orbit)
full = orbit.full
G = full.old.G

def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def laplacian(field, xyz):
    r = np.linalg.norm(xyz, axis=1)
    theta = np.arctan2(np.hypot(xyz[:, 0], xyz[:, 1]), xyz[:, 2])
    phi = np.arctan2(xyz[:, 1], xyz[:, 0])
    leg = sph_legendre_p_all(int(field.l.max()), int(field.m.max()), theta, diff_n=0)
    Y = leg[0, field.l, field.m]*np.cos(field.m[:, None]*phi)*np.where(field.m[:, None] > 0, np.sqrt(2), 1.)
    log = np.log(r)
    terms = field.spline(log, 2)+field.spline(log, 1)-field.l*(field.l+1)*field.spline(log)
    return (terms*Y.T).sum(axis=1)/r**2

def divergence(fun, xyz, h):
    out = np.zeros(len(xyz))
    for j in range(3):
        delta = np.zeros(3)
        delta[j] = h
        out += (fun(xyz+delta)[:, j]-fun(xyz-delta)[:, j])/(2*h)
    return out

def chunks(fun, xyz):
    return np.concatenate([fun(part) for part in np.array_split(xyz, max(1, (len(xyz)+15)//16))])

def main():
    ordinary = orbit.Field('ordinary')
    axispath = full.CACHE/'axisymmetric-reference-80-L128-N512.npz'
    finepath, finerpath = full.CACHE/'fine.npz', full.CACHE/'finer.npz'
    fields = [full.load_field(p) for p in [axispath, finepath, finerpath]]
    axis, fine, finer = fields
    files = [Path(__file__), DRIVER, HERE.parent/'full-bar-completion/run.py', axispath, finepath, finerpath, *ordinary.paths]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in files}
    locations = [(R,z,phi) for R in [.5,1,2,3,5,8,12,20] for z in [0,.1,.5,1,2,4] for phi in np.linspace(0,np.pi/2,5)]
    xyz = np.array([[R*np.cos(phi), R*np.sin(phi), z] for R,z,phi in locations])
    def density(field):
        return chunks(lambda x: (laplacian(axis,x)+laplacian(field,x))/(4*np.pi*G),xyz)
    rho_fine, rho_finer = density(fine), density(finer)
    def extra(x):
        return axis.evaluate(x)[1]+finer.evaluate(x)[1]
    def Q(x):
        grad = -ordinary.evaluate(x)[1]
        mult = full.old.A*(np.linalg.norm(grad,axis=1)/full.old.ASTAR)**(full.old.P-1)
        return mult[:,None]*grad
    steps = [.004,.002,.001]
    lap_fd = [chunks(lambda x: -divergence(extra,x,h)/(4*np.pi*G),xyz) for h in steps]
    source_fd = [chunks(lambda x: divergence(Q,x,h)/(4*np.pi*G),xyz) for h in steps]
    assert np.isfinite(np.array([rho_fine,rho_finer,*lap_fd,*source_fd])).all()
    # Control the divergence sign and normalization with a known quadratic potential.
    np.testing.assert_allclose(divergence(lambda x: -2*x,xyz,.001), -6, rtol=1e-10)
    rows=[]
    for i,(R,z,phi) in enumerate(locations):
        rows.append(dict(R_kpc=R,z_kpc=z,phi_rad=float(phi),xyz_kpc=xyz[i].tolist(),
            rho_fine_Msun_kpc3=float(rho_fine[i]),rho_finer_Msun_kpc3=float(rho_finer[i]),
            rho_from_acceleration_steps=[float(v[i]) for v in lap_fd],
            rho_from_requested_divQ_steps=[float(v[i]) for v in source_fd],
            harmonic_resolution_difference=float(abs(rho_finer[i]-rho_fine[i])),
            requested_minus_reconstructed=float(source_fd[-1][i]-rho_finer[i])))
    out=dict(scope='Conditional source audit of the full-bar empirical potential; not an energy budget or observational fit.',
        samples=len(rows),steps_kpc=steps,rows=rows,
        negative_reconstructed_finer=int((rho_finer<0).sum()),negative_requested_divQ=int((source_fd[-1]<0).sum()),
        reconstructed_range_Msun_kpc3=[float(rho_finer.min()),float(rho_finer.max())],
        requested_range_Msun_kpc3=[float(source_fd[-1].min()),float(source_fd[-1].max())],
        source_hashes=hashes,holdouts_opened=False,photon_deposit_identification_proved=False)
    for p in files:
        assert digest(p)==hashes[str(p.relative_to(ROOT))]
    (HERE/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['rows','source_hashes']},indent=2))

if __name__=='__main__':
    main()

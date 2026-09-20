"""RW-1 caches: the whirlpool columns of every SPARC galaxy (all splits) on the declared kernel grid, with the K7 doublings;
the cluster pressure columns on stage 2's 6000-point grid with the K8 doublings; the Milky Way columns. Writes
research_work/generated/routes/rw1_galaxies.npz, rw1_clusters.npz, rw1_milky_way.npz. Usage: python rw1_build.py {galaxies|clusters|milky_way}"""
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cl2_sources as CS      # noqa: E402
import cl2s2_lib as L2        # noqa: E402
import cl2s2 as D             # noqa: E402
import rw1_lib as RW          # noqa: E402

GEN = HERE.parents[2]/'research_work/generated/routes'
GEN.mkdir(parents=True, exist_ok=True)


def galaxies():
    t0 = time.time()
    store = {'kernels': np.array(RW.KERNELS_GAL)}
    gals = sorted(CS.I.sparc_galaxies(), key=lambda g: (g['split'], g['name']))
    for i, gal in enumerate(gals):
        rec = L2.reconstructed_newton(gal)
        R = rec['R']
        cols, src = RW.galaxy_whirl_columns(gal, R)
        pcols, _ = RW.galaxy_whirl_columns(gal, R, point_source=True)
        key = f"{gal['split']}/{gal['name']}/"
        store[key + 'cols'] = cols
        store[key + 'point'] = pcols
        store[key + 'R'] = R
        store[key + 'geometry'] = np.array([src['R_half'], src['M_grid'], src['rd'], rec['mass']])
        if gal['name'] in ('NGC2403', 'NGC3198', 'DDO154'):
            c256, _ = RW.galaxy_whirl_columns(gal, R, nodes=RW.N256)
            c6000, _ = RW.galaxy_whirl_columns(gal, R, n_grid=6000)
            store[key + 'nodes256'] = c256
            store[key + 'grid6000'] = c6000
        if i % 10 == 0:
            print(f'{i + 1}/{len(gals)} {gal["name"]} done at {time.time() - t0:.0f} s', flush=True)
    np.savez_compressed(GEN/'rw1_galaxies.npz', **store)
    print('galaxies saved at %.0f s' % (time.time() - t0), flush=True)


def clusters():
    t0 = time.time()
    ext = CS.load_xcop()
    frac = CS.stellar_fraction_profile(ext)
    cls, _ = D.load_clusters(ext, frac)
    store = {'kernels': np.array(RW.KERNELS_CL)}
    for name, cl in cls.items():
        ops, gs = RW.cluster_whirl_columns(cl)
        store[name + '/ops'] = ops
        store[name + '/g'] = gs
        store[name + '/rp'] = cl['rp']
        store[name + '/r'] = cl['r']
        print(f'{name} done at {time.time() - t0:.0f} s', flush=True)
    # K8: the full-grid evaluation and the 12000-point source on A1795 and A2255 at w in {1, 16, 256} and every p
    sub = [(p, w) for p in RW.P_GRID for w in (1., 16., 256.)]
    for name in ('A1795', 'A2255'):
        cl = cls[name]
        ops_i, _ = RW.cluster_whirl_columns(cl, sub)
        ops_f, _ = RW.cluster_whirl_columns(cl, sub, full=True)
        cl12 = CS.build_cluster(name, ext['clusters'][name], np.array([100.]), n_grid=12000, frac_profile=frac)
        ops_12, _ = RW.cluster_whirl_columns(cl12, sub)
        store[name + '/K8_sub'] = np.array(sub)
        store[name + '/K8_interp'] = ops_i
        store[name + '/K8_full'] = ops_f
        store[name + '/K8_grid12000'] = ops_12
        print(f'{name} K8 done at {time.time() - t0:.0f} s', flush=True)
    np.savez_compressed(GEN/'rw1_clusters.npz', **store)
    print('clusters saved at %.0f s' % (time.time() - t0), flush=True)


def milky_way():
    t0 = time.time()
    store = {'kernels': np.array(RW.KERNELS_GAL)}
    for b in CS.milky_way_blocks(np.array([1.])):
        variant = b.meta['variant']
        R = b.meta['R']
        sigma, m_bulge = CS.milky_way_source(variant)
        Rg = np.geomspace(1e-3, 300., 3000)
        sig = sigma(Rg)
        Rm = np.sqrt(Rg[1:]*Rg[:-1])
        mass = np.pi*(Rg[1:]**2 - Rg[:-1]**2)*.5*(sig[1:] + sig[:-1])
        cols = RW.whirl_columns_disk(R, Rm, mass, RW.KERNELS_GAL)
        cum = np.concatenate([[0.], np.cumsum(mass)])
        if m_bulge is not None:
            Mb = m_bulge(Rg)
            cols = cols + np.array([RW.shell_apply(R, Rm, np.diff(Mb), p, w) for (p, w) in RW.KERNELS_GAL]).T
            cum = cum + Mb
        store[variant + '/cols'] = cols
        store[variant + '/R'] = R
        store[variant + '/y'] = b.meta['y']
        store[variant + '/gN'] = b.meta['gN']
        store[variant + '/geometry'] = np.array([float(np.interp(.5*cum[-1], cum, Rg)), float(cum[-1])])
        print(f'Milky Way {variant} done at {time.time() - t0:.0f} s', flush=True)
    np.savez_compressed(GEN/'rw1_milky_way.npz', **store)


if __name__ == '__main__':
    {'galaxies': galaxies, 'clusters': clusters, 'milky_way': milky_way}[sys.argv[1]]()

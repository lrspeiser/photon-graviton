"""Acquire the X-COP public data release and write the compact extract CL-2 uses (protocol-cl2.md, CL-2A).

    python cl2_xcop_acquire.py [--from-file PATH]

The release is one archive shared from the X-COP data page (https://dominiqueeckert.wixsite.com/xcop/data):
per cluster, the electron-density shells (the L1-penalty deprojection), the X-ray and X/SZ temperature
profiles, the X-ray and SZ pressure profiles, the hydrostatic mass reconstructions, the gas-mass and
gas-fraction profiles, and for seven clusters the cumulative stellar-mass profiles of Ghizzardi et al.
The archive is 315 MB, mostly image mosaics, and is cached under research_work/generated (ignored); only
the profile tables are read, and only their numbers are written to cl2-inputs-xcop-profiles.json, which is
committed. The archive's SHA-256 is pinned below; a different archive is refused. Needs astropy.

Every quantity here is a data-release product with its own modelling (deprojection, spectral fitting,
hydrostatic reconstruction, the release's H0 = 70 distances); the extract records that, and CL-2 never
treats the hydrostatic masses as raw observations.
"""
import argparse
import hashlib
import io
import json
import os
import shutil
import sys
import tarfile
import time
import urllib.request
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT/'research_work/generated/xcop-release'
OUT = HERE/'cl2-inputs-xcop-profiles.json'
URL = 'https://drive.switch.ch/index.php/s/j3WUOYXWgv9Jbnz/download'
SHA256 = '0edf5038b419b70d070b73b22f4801e27f318b0854db61eec52142c27c140d94'
SIZE = 315080566
PRODUCTS = ('density_L1', 'temperature', 'pressure', 'hydro_mass', 'fgas_profile', 'mstar')
SOURCES = dict(
    release='X-COP data release (D. Eckert and collaborators), the archive linked from the data page',
    density='Ghirardini et al. 2019, A&A 621, A41 (arXiv:1805.00042): electron density by the L1-penalty deprojection',
    thermodynamics='Ghirardini et al. 2019: X-ray temperature and pressure, SZ pressure (Planck), X/SZ temperature; self-similar scalings T500, P500 in the headers',
    masses='Ettori et al. 2019, A&A 621, A39 (arXiv:1805.00035): hydrostatic mass profiles for several mass models (forward, NFW, Einasto)',
    stars='Ghizzardi et al. 2021 (as cited in the release headers): cumulative stellar mass profiles, seven clusters',
    cosmology='the release\'s own H0 = 70, Omega_m = 0.3 distances, adopted as scenario inputs (universe contract)')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def fetch(from_file=None):
    CACHE.mkdir(parents=True, exist_ok=True)
    target = CACHE/'allfiles.tar.gz'
    if target.exists() and sha(target) == SHA256:
        return target
    if from_file:
        shutil.copyfile(from_file, target)
    else:
        print('downloading', URL, flush=True)
        urllib.request.urlretrieve(URL, target)
    digest = sha(target)
    if digest != SHA256:
        raise SystemExit(f'archive hash {digest} does not match the pinned {SHA256}; refusing to use it')
    return target


def col(t, name):
    return [float(v) for v in t.data[name]]


def read_cluster(tar, name):
    import astropy.io.fits as F
    tables = {}
    for prod in PRODUCTS:
        member = f'{name}/{name}_{prod}.fits'
        try:
            f = tar.extractfile(member)
        except KeyError:
            f = None
        if f is None:
            continue
        tables[prod] = F.open(io.BytesIO(f.read()))
    h = tables['density_L1']
    hd = h['DENSITY'].header
    c = dict(header=dict(z=float(hd['REDSHIFT']), R500_kpc=float(hd['R500']), R500_err=float(hd['ERR_R500']),
                         M500_1e14=float(hd['M500']), M500_err=float(hd['ERR_M500'])))
    t = h['DENSITY']
    c['density'] = dict(r_in_kpc=col(t, 'R_IN'), r_out_kpc=col(t, 'R_OUT'), ne_cm3=col(t, 'NE'), ne_lo=col(t, 'NE_LOW'), ne_hi=col(t, 'NE_HIGH'),
                        note='electron density in shells, L1-penalty deprojection; lo/hi are the release\'s bounds')
    h = tables['temperature']
    c['temperature'] = dict(T500_keV=float(h['XRAY'].header['T500']),
                            xray=dict(r_over_R500=col(h['XRAY'], 'RW_X'), T_over_T500=col(h['XRAY'], 'T_X'), err=col(h['XRAY'], 'eT_X')),
                            xsz=dict(r_over_R500=col(h['XSZ'], 'RW_SZ'), T_over_T500=col(h['XSZ'], 'T_SZ'), err=col(h['XSZ'], 'eT_SZ')))
    h = tables['pressure']
    c['pressure'] = dict(P500_keV_cm3=float(h['XRAY'].header['P500']),
                         xray=dict(r_over_R500=col(h['XRAY'], 'RW_X'), P_over_P500=col(h['XRAY'], 'P_X'), err=col(h['XRAY'], 'eP_X')),
                         sz=dict(r_over_R500=col(h['SZ'], 'RW_SZ'), P_over_P500=col(h['SZ'], 'P_SZ'), err=col(h['SZ'], 'eP_SZ')),
                         note='electron pressure scaled by P500; X-ray points from n_e kT, SZ points from Planck')
    h = tables['hydro_mass']
    t = h['HYDRO_MASS']
    c['hydro_mass'] = {k: col(t, k) for k in ('RADIUS', 'M_FORW', 'EM_FORW', 'M_NFW', 'EM_NFW', 'M_EIN', 'EM_EIN')}
    c['hydro_params'] = [dict(model=str(r['MODEL']).strip(), rs_kpc=float(r['RS']), rs_err=float(r['eRS']), c200=float(r['C200']), c200_err=float(r['eC200']))
                         for r in h['PARAMS'].data]
    t = tables['fgas_profile']['FGAS']
    c['gas_mass'] = {k: col(t, k) for k in ('RADIUS', 'MGAS', 'MGAS_LO', 'MGAS_HI', 'M_NFW', 'FGAS', 'FGAS_LO', 'FGAS_HI')}
    c['gas_mass']['note'] = 'RADIUS in units of R500; MGAS the release\'s gas mass, M_NFW its hydrostatic NFW mass'
    if 'mstar' in tables:
        t = tables['mstar']['MSTAR_SMOOTHED']
        idx = np.unique(np.concatenate([[0], np.linspace(0, len(t.data) - 1, 80).astype(int)]))
        c['stellar_mass'] = dict(radius_kpc=[float(t.data['RADIUS'][i]) for i in idx], Mstar=[float(t.data['MSTAR'][i]) for i in idx],
                                 Mstar_lo=[float(t.data['MSTAR_LO'][i]) for i in idx], Mstar_hi=[float(t.data['MSTAR_HI'][i]) for i in idx],
                                 rows_in_release=int(len(t.data)), table='MSTAR_SMOOTHED', note='cumulative stellar mass, smoothed, total uncertainties')
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--from-file', help='use this copy of the archive instead of downloading (it must match the pinned hash)')
    args = ap.parse_args()
    target = fetch(args.from_file)
    out = dict(scope='Compact extract of the X-COP data release for CL-2A; every number is a release product, not a raw observation',
               provenance=dict(url=URL, sha256=SHA256, size_bytes=SIZE, products=list(PRODUCTS), sources=SOURCES,
                               extracted=time.strftime('%Y-%m-%d'), mosaics_and_images='not read'),
               units=dict(radius='kpc unless a column is named r_over_R500', ne='cm^-3', T='T500 units, T500 in keV', P='P500 units, P500 in keV cm^-3',
                          masses='Msun', molecular_weights_used_by_the_release='mu = 0.6 per particle, mu_e = 1.14 per electron'),
               clusters={})
    with tarfile.open(target, 'r:gz') as tar:
        names = sorted({m.name.split('/')[0] for m in tar.getmembers() if '/' in m.name and m.name.endswith('_density_L1.fits')})
        for name in names:
            out['clusters'][name] = read_cluster(tar, name)
            c = out['clusters'][name]
            print(name, 'z', c['header']['z'], 'R500', c['header']['R500_kpc'], 'shells', len(c['density']['ne_cm3']),
                  'pressure points', len(c['pressure']['xray']['r_over_R500']) + len(c['pressure']['sz']['r_over_R500']),
                  'stars', 'measured' if 'stellar_mass' in c else 'none', flush=True)
    text = json.dumps(out, indent=0) + '\n'
    OUT.write_text(text, encoding='utf-8', newline='\n')
    print('wrote', OUT, len(text), 'bytes; clusters', len(out['clusters']))


if __name__ == '__main__':
    main()

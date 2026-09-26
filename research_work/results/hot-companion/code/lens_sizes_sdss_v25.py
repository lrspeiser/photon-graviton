"""Round 25: how large are the galaxies that lens in KiDS? Their sizes, measured on the same SDSS galaxies as their heat.

The trapped-companion candidate (round25-trapped-suite.md) needs each lens's size: the self-energy repair keeps a
galaxy's hot glow inside the region where its own pull is strong enough, so the share that escapes depends on how
concentrated the stars are. This script repeats the selection of code/lens_heat_sdss_v12.py (SDSS DR17 spectra,
0.1 < z < 0.25, MPA-JHU stellar masses 10^9.8-10^11.3, the same red / blue and bulge / disc splits) and adds the
galaxies' sizes:
    early types (red, bulge)   the de Vaucouleurs fit's half-light radius, circularised: deVRad_r sqrt(deVAB_r)
    late types (blue, disc)    the exponential fit's scale length on the major axis: expRad_r / 1.678 (SDSS gives
                               the half-light radius)
converted to kpc with the stellar masses' own cosmology (flat, H0 = 70, Omega_m = 0.3). The medians per sample and
0.1-dex bin of log M* go to data/lens_sizes_sdss_v25.json. (TOP without ORDER BY, as in round 12: the server returns
its first 120,000 matches, so the set can differ slightly from round 12's; medians of thousands per bin.)

    python code/lens_sizes_sdss_v25.py --output data/lens_sizes_sdss_v25.json
"""
from __future__ import annotations
import argparse, csv, io, json, urllib.parse, urllib.request
from pathlib import Path
import numpy as np

SQL = ("SELECT TOP 120000 s.z, s.velDisp, s.velDispErr, g.lgm_tot_p50, p.modelMag_u, p.modelMag_r, p.fracDeV_r, p.deVRad_r, "
       "p.deVAB_r, p.expRad_r, p.expAB_r "
       "FROM SpecObj s JOIN galSpecExtra g ON g.specObjID = s.specObjID JOIN PhotoObj p ON p.objID = s.bestObjID "
       "WHERE s.class = 'GALAXY' AND s.zWarning = 0 AND s.z BETWEEN 0.1 AND 0.25 AND g.lgm_tot_p50 BETWEEN 9.8 AND 11.3 "
       "AND s.velDisp > 0 AND s.velDispErr > 0 AND s.velDispErr < 40")
URL = 'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch'
BINS = np.round(np.arange(9.8, 11.3001, 0.1), 2)
C_KMS = 299792.458


def kpc_per_arcsec(z, H0=70.0, Om=0.3):
    zz = np.linspace(0.0, 0.3, 3001)
    Ez = np.sqrt(Om * (1 + zz) ** 3 + 1 - Om)
    Dc = C_KMS / H0 * np.concatenate([[0.0], np.cumsum(0.5 * (1 / Ez[1:] + 1 / Ez[:-1]) * np.diff(zz))])   # Mpc
    DA = np.interp(z, zz, Dc) / (1 + z)
    return DA * 1e3 / 206264.80624709636


def fetch(cache: Path | None):
    if cache is not None and cache.exists():
        text = cache.read_text()
    else:
        q = urllib.parse.urlencode(dict(cmd=SQL, format='csv'))
        text = urllib.request.urlopen(f'{URL}?{q}', timeout=900).read().decode()
        if cache is not None:
            cache.write_text(text)
    rows = list(csv.DictReader(io.StringIO('\n'.join(text.splitlines()[1:]))))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def binned(d):
    ur = d['modelMag_u'] - d['modelMag_r']
    ok = np.isfinite(ur) & (ur > 0) & (ur < 5) & (d['deVRad_r'] > 0.3) & (d['velDisp'] > 0)     # round 12's cuts
    f = np.clip(d['fracDeV_r'], 0, 1)
    kpc = kpc_per_arcsec(d['z'])
    Re = d['deVRad_r'] * np.sqrt(np.clip(d['deVAB_r'], 0.05, 1.0)) * kpc
    Rd = d['expRad_r'] / 1.678 * kpc
    samples = dict(red=ok & (ur > 2.5), blue=ok & (ur <= 2.5), bulge=ok & (f > 0.5), disc=ok & (f <= 0.5))
    out = {}
    for name, sel in samples.items():
        rows = []
        for lo in BINS[:-1]:
            m = sel & (d['lgm_tot_p50'] >= lo) & (d['lgm_tot_p50'] < lo + 0.1) & (d['expRad_r'] > 0)
            if m.sum() < 20:
                continue
            rows.append(dict(logM_lo=float(lo), n=int(m.sum()), median_Re_deV_circ_kpc=float(np.median(Re[m])),
                             median_Rd_exp_kpc=float(np.median(Rd[m])), median_fracDeV=float(np.median(f[m]))))
        out[name] = rows
    return out


def size(table, sample, logM):
    """The size the trapped-companion candidate uses for a lens of this sample and log M*: early types (red, bulge)
    a Hernquist profile with the circularised de Vaucouleurs half-light radius, late types (blue, disc) an exponential
    disk with the exponential scale length (nearest 0.1-dex bin)."""
    rows = table['samples'][sample]
    r = min(rows, key=lambda r: abs(r['logM_lo'] + 0.05 - logM))
    if sample in ('red', 'bulge'):
        return 'hernquist', r['median_Re_deV_circ_kpc']
    return 'disk', r['median_Rd_exp_kpc']


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--table', type=Path, default=None, help='a cached copy of the SDSS query result (CSV)')
    args = ap.parse_args()
    d = fetch(args.table)
    res = dict(experiment='round 25: the sizes of KiDS-like lenses, from SDSS DR17 (the galaxies whose heat round 12 measured)', sql=SQL,
               n_galaxies=int(len(d['z'])), cosmology='flat, H0 = 70, Omega_m = 0.3 (the MPA-JHU masses\' units)',
               early='deVRad_r sqrt(deVAB_r): circularised half-light radius of the de Vaucouleurs fit',
               late='expRad_r / 1.678: scale length of the exponential fit (major axis)', samples=binned(d))
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    for s, rows in res['samples'].items():
        print(s, ' '.join(f"{r['logM_lo']:.1f}:{(r['median_Re_deV_circ_kpc'] if s in ('red', 'bulge') else r['median_Rd_exp_kpc']):.2f}" for r in rows))
    print(f"wrote {args.output} ({res['n_galaxies']} galaxies)")


if __name__ == '__main__':
    main()

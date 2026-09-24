"""Round 12: how hot are the galaxies that lens in KiDS? Their stars' speeds, measured.

Our law gives each lens a heat weight k = 3 sigma^2 / u^2 for the part of its stars that moves at
random, and k about 0.1 for a cold disk (sigma about 30 km/s). Until round 11 the KiDS and Mistele
comparisons assumed one number per class: 160 km/s for every star of a red (or bulge-dominated) lens
and k = 0.1 for every blue (or disk-dominated) one. This script measures instead, from SDSS DR17 spectra
of galaxies with the same stellar masses and the same splits:
    red / blue           observed u - r above / below 2.5 (Brouwer et al. 2021's KiDS split)
    bulge / disc         fracDeV_r above / below 0.5 (a proxy for their Sersic n above / below 2)
For each galaxy the hot part is fracDeV_r of its light, with the fibre dispersion corrected to one
effective radius (sigma_R / sigma_e = (R / R_e)^-0.066, Cappellari et al. 2006); the rest is a disk
at 30 km/s. Stellar masses are the MPA-JHU values (Kroupa, close to Chabrier), 0.1 < z < 0.25.

Output: for each sample and 0.1-dex bin of log M*, the mean of fracDeV x sigma_e^2 and of
(1 - fracDeV); k_eff(u) = 3 [<f sigma_e^2> + <1 - f> 30^2] / u^2. Written to
data/lens_heat_sdss_v12.json (the SDSS table is cached next to it only if --keep-table is given).

    python code/lens_heat_sdss_v12.py --output ../data/lens_heat_sdss_v12.json
"""
from __future__ import annotations
import argparse, csv, io, json, urllib.parse, urllib.request
from pathlib import Path
import numpy as np

SQL = ("SELECT TOP 120000 s.z, s.velDisp, s.velDispErr, g.lgm_tot_p50, p.modelMag_u, p.modelMag_r, p.fracDeV_r, p.deVRad_r, p.expRad_r "
       "FROM SpecObj s JOIN galSpecExtra g ON g.specObjID = s.specObjID JOIN PhotoObj p ON p.objID = s.bestObjID "
       "WHERE s.class = 'GALAXY' AND s.zWarning = 0 AND s.z BETWEEN 0.1 AND 0.25 AND g.lgm_tot_p50 BETWEEN 9.8 AND 11.3 "
       "AND s.velDisp > 0 AND s.velDispErr > 0 AND s.velDispErr < 40")
URL = 'https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch'
SIGMA_DISK = 30.0
BINS = np.round(np.arange(9.8, 11.3001, 0.1), 2)


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
    ok = np.isfinite(ur) & (ur > 0) & (ur < 5) & (d['deVRad_r'] > 0.3) & (d['velDisp'] > 0)
    sig_e = d['velDisp'] * (1.5 / np.where(ok, d['deVRad_r'], 1.0)) ** 0.066
    f = np.clip(d['fracDeV_r'], 0, 1)
    samples = dict(red=ok & (ur > 2.5), blue=ok & (ur <= 2.5), bulge=ok & (f > 0.5), disc=ok & (f <= 0.5))
    out = {}
    for name, sel in samples.items():
        rows = []
        for lo in BINS[:-1]:
            m = sel & (d['lgm_tot_p50'] >= lo) & (d['lgm_tot_p50'] < lo + 0.1)
            if m.sum() < 20:
                continue
            rows.append(dict(logM_lo=float(lo), n=int(m.sum()), mean_f_sigma2=float(np.mean(f[m] * sig_e[m] ** 2)),
                             mean_1_minus_f=float(np.mean(1 - f[m])), median_sigma_e=float(np.median(sig_e[m]))))
        out[name] = rows
    return out


def k_eff(table, sample, u_kms, logM=None, lo=10.3, hi=10.9):
    """Mean heat weight of a sample: at one log M* (nearest bin), or averaged over bins in [lo, hi)."""
    rows = table['samples'][sample]
    import law as L
    def k(r):
        if L.HEAT_P == 2.0:
            return 3 * (r['mean_f_sigma2'] + r['mean_1_minus_f'] * table['sigma_disk_kms'] ** 2) / u_kms ** 2
        # round 14, k = 3 (sigma / u)^p: bulge and disk weighted separately, the bulge at its f-weighted rms speed
        f = max(1.0 - r['mean_1_minus_f'], 1e-9)
        return f * L.heat_weight(np.sqrt(r['mean_f_sigma2'] / f), u_kms) + r['mean_1_minus_f'] * L.heat_weight(table['sigma_disk_kms'], u_kms)
    if logM is not None:
        r = min(rows, key=lambda r: abs(r['logM_lo'] + 0.05 - logM))
        return k(r)
    sel = [r for r in rows if lo - 1e-9 <= r['logM_lo'] < hi - 1e-9]
    return float(np.mean([k(r) for r in sel]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--table', type=Path, default=None, help='a cached copy of the SDSS query result (CSV)')
    args = ap.parse_args()
    d = fetch(args.table)
    res = dict(experiment='round 12: the stars\' speeds of KiDS-like lenses, from SDSS DR17 spectra', sql=SQL, n_galaxies=int(len(d['z'])),
               sigma_disk_kms=SIGMA_DISK, aperture_correction='sigma_e = sigma_fibre (1.5 arcsec / R_e)^0.066 (Cappellari et al. 2006)',
               samples=binned(d))
    for u in (162.59, 197.41):
        res.setdefault('k_eff_10p3_10p9', {})[f'u={u}'] = {s: k_eff(res, s, u) for s in res['samples']}
    args.output.write_text(json.dumps(res, indent=1) + '\n')
    for u, v in res['k_eff_10p3_10p9'].items():
        print(u, {s: round(x, 2) for s, x in v.items()})
    print('wrote', args.output)


if __name__ == '__main__':
    main()

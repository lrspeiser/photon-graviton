"""Read-only catalog audit; recover BRAVA into ignored cache if absent."""
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd
from astropy.io import fits
from astropy.table import Table

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT / 'research_work/data-cache/stellar-catalogs'

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

records = []
for expected in json.loads((HERE / 'download-manifest.json').read_text()):
    path = CACHE / expected['file']
    actual = digest(path)
    assert actual == expected['sha256'], path
    assert path.stat().st_size == expected['bytes'], path
    with fits.open(path, memmap=True) as hdus:
        rows = hdus[1].header['NAXIS2']
    records.append(dict(path=str(path), rows=rows, sha256=actual,
                        original_download_hash_matches=True))

inventory = json.loads((HERE / 'local-inventory.json').read_text())
brava = next(row for row in inventory if row['path'].endswith('brava_catalog.tbl'))
source = Path(brava['path'])
assert digest(source) == brava['sha256']
destination = CACHE / 'brava_catalog.tbl'
if not destination.exists():
    shutil.copy2(source, destination)
assert digest(destination) == brava['sha256']
rows = len(Table.read(destination, format='ascii.ipac'))
assert rows == brava['rows']
records.append(dict(path=str(destination), rows=rows,
                    sha256=brava['sha256'], recovered_from=str(source)))

path = CACHE / 'gaia-dr3-quality-covariance.fits'
with fits.open(path, memmap=True) as hdus:
    gaia_rows = hdus[1].header['NAXIS2']
records.append(dict(path=str(path), rows=gaia_rows, sha256=digest(path)))
path = CACHE / 'matched-with-gaia-covariance.parquet'
data = pd.read_parquet(path, columns=['source_id', 'gaia_quality_candidate', 'R_kpc', 'z_kpc'])
assert data.source_id.is_unique and data.source_id.dtype.kind in 'iu'
assert len(data) == gaia_rows == 140407
quality = data.loc[data.gaia_quality_candidate]
regions = {}
for label, low, high, zlow, zhigh in [
    ('plane_beneath_bulge', .5, 3.5, 0, .2),
    ('above_below_bulge', .5, 3.5, .5, 1.5),
    ('disk_plane_control', 5, 9, 0, .2),
    ('disk_off_plane_control', 5, 9, .5, 1.5),
]:
    mask = ((quality.R_kpc >= low) & (quality.R_kpc < high)
            & (quality.z_kpc.abs() >= zlow) & (quality.z_kpc.abs() < zhigh))
    regions[label] = int(mask.sum())
result = dict(audit_date='2026-09-10', catalogs=records,
              prepared_rows=len(data), quality_candidates=len(quality),
              unique_integer_source_ids=True, region_coverage=regions,
              prepared_sha256=digest(path),
              scope='File integrity and coverage only; no new dynamical fit or holdout score.')
with (HERE / 'local-recheck.json').open('w', newline='\n') as stream:
    json.dump(result, stream, indent=2)
    stream.write('\n')
print(json.dumps(result, indent=2))

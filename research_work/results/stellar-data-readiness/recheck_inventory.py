"""Recheck local file integrity without opening held-out stellar outcomes."""
from pathlib import Path
import hashlib
import json
from astropy.io import fits
from astropy.table import Table
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
CACHE = HERE.parents[1] / 'data-cache' / 'stellar-catalogs'
previous = json.loads((HERE / 'local-recheck.json').read_text())
records = []
for old in previous['catalogs']:
    path = CACHE / Path(old['path']).name
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    if path.suffix == '.fits':
        with fits.open(path, memmap=True) as hdus:
            rows = int(hdus[1].header['NAXIS2'])
    else:
        rows = len(Table.read(path, format='ascii.ipac'))
    records.append(dict(file=path.name, bytes=path.stat().st_size, rows=rows,
                        sha256=digest, matches_previous=digest == old['sha256']))
path = CACHE / 'matched-with-gaia-covariance.parquet'
with path.open('rb') as stream:
    digest = hashlib.file_digest(stream, 'sha256').hexdigest()
result = dict(catalogs=records, prepared_rows=pq.read_metadata(path).num_rows,
              prepared_sha256=digest,
              prepared_matches_previous=digest == previous['prepared_sha256'],
              scope='Integrity and metadata only; no new holdout outcomes or gravity fit.')
assert all(r['matches_previous'] for r in records)
assert result['prepared_matches_previous']
(HERE / 'inventory-recheck.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
print(json.dumps(result, indent=2))

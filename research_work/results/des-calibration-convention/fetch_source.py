"""Download pinned source as read-only evidence, never execute it."""
from pathlib import Path
import hashlib,json,urllib.request
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
m=json.loads((HERE/'source-manifest.json').read_text(encoding='utf-8'))
out=ROOT/'research_work/generated/des-photometric-calibration/snana-source';out.mkdir(parents=True,exist_ok=True)
for row in m['files']:
    data=urllib.request.urlopen(row['url'],timeout=60).read()
    assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256']
    (out/row['path']).write_bytes(data)
print('Verified four pinned source files; no source code executed.')

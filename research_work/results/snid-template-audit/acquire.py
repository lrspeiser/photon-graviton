"""Acquire official SNID archives for inspection; never run downloaded code."""
from pathlib import Path
import hashlib
import io
import json
import tarfile
import urllib.request

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = ROOT/'research_work/generated/snid-template-audit'
CACHE.mkdir(parents=True, exist_ok=True)
base = 'https://people.lam.fr/blondin.stephane/software/snid/'
archives = []
for filename in ['templates-2.0.tgz', 'snid-5.0.tar.gz']:
    destination = CACHE/filename
    if not destination.exists():
        with urllib.request.urlopen(base+filename, timeout=60) as response:
            blob = response.read(30_000_001)
        assert len(blob) <= 30_000_000
        destination.write_bytes(blob)
    blob = destination.read_bytes()
    members = []
    with tarfile.open(fileobj=io.BytesIO(blob), mode='r:gz') as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            raw = archive.extractfile(member).read()
            row = dict(name=member.name, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            if member.name.endswith('.lnw'):
                text = raw.decode('ascii', errors='replace')
                row['header'] = text.splitlines()[0]
            members.append(row)
    archives.append(dict(url=base+filename, cache=str(destination.relative_to(ROOT)),
                         bytes=len(blob), sha256=hashlib.sha256(blob).hexdigest(), members=members))
out = dict(source_page=base+'index.html', archives=archives,
           scope='Archive inventory and template headers only; no downloaded program executed',
           script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'manifest.json').write_text(json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
print(json.dumps([dict(url=a['url'],bytes=a['bytes'],files=len(a['members']),templates=sum('header' in m for m in a['members'])) for a in archives]))

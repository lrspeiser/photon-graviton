"""Restore the losslessly compressed large covariance input before legacy analyses."""
from pathlib import Path
import gzip,hashlib
p=Path(__file__).resolve().parent/'shared_interaction_test/data/Pantheon+SH0ES_STAT+SYS.cov'
expected=next(x['sha256'] for x in __import__('json').loads((Path(__file__).parent/'SNAPSHOT_MANIFEST.json').read_text())['files'] if x['path']==str(p.relative_to(Path(__file__).parent)))
if not p.exists():
 data=gzip.decompress(p.with_suffix(p.suffix+'.gz').read_bytes())
 if hashlib.sha256(data).hexdigest()!=expected:raise RuntimeError('Covariance checksum mismatch')
 p.write_bytes(data)
if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise RuntimeError('Covariance checksum mismatch')
print('Covariance input restored and verified.')

"""Read template arrays and identity overlap, without invoking SNID."""
from pathlib import Path
import csv
import hashlib
import io
import json
import tarfile
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
manifest=json.loads((HERE/'manifest.json').read_text(encoding='utf-8'))
catalog=HERE.parent/'electromagnetic-audit/spectral-aging-predictions.csv'
with catalog.open(newline='',encoding='utf-8') as f:
    objects=list(csv.DictReader(f))
sets=[]
for entry, prefix in [(manifest['archives'][0],'templates-2.0/'),
                      (manifest['archives'][1],'snid-5.0/templates/')]:
    path=ROOT/entry['cache']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256']
    templates=[]
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            if not member.name.startswith(prefix) or not member.name.endswith('.lnw'):
                continue
            lines=archive.extractfile(member).read().decode('ascii').splitlines()
            header=lines[0].split()
            nspec,npix,nknots=int(header[0]),int(header[1]),int(header[4])
            age=np.array([float(v) for v in lines[nknots+2].split()])
            spectra=np.loadtxt(io.StringIO('\n'.join(lines[nknots+3:])))
            assert len(age)==nspec+1
            assert spectra.shape==(npix,nspec+1),member.name
            assert np.isfinite(spectra).all() and np.isfinite(age).all()
            assert np.all(np.diff(spectra[:,0])>0)
            templates.append(dict(file=member.name,object=header[5],subtype=header[7],
                                  spectra=nspec,pixels=npix,age_flag=float(age[0]),
                                  phase_min=float(age[1:].min()),phase_max=float(age[1:].max()),
                                  phases=age[1:].tolist()))
    overlaps=[]
    for obj in objects:
        name=obj['object'].lower()
        candidates={'sn'+name}
        if name[:2] in ['19','20']:
            candidates.add('sn'+name[2:])
        matches=[t['file'] for t in templates if t['object'].lower() in candidates]
        overlaps.append(dict(object=obj['object'],z=float(obj['z']),template_files=matches))
    sets.append(dict(prefix=prefix,templates=templates,template_count=len(templates),
                     spectral_columns=sum(t['spectra'] for t in templates),
                     overlap_with_published_aging_sample=overlaps,
                     overlapping_objects=sum(bool(o['template_files']) for o in overlaps)))
out=dict(scope='Structural template audit; phases not independently recalibrated',sets=sets,
         hashes={'manifest.json':hashlib.sha256((HERE/'manifest.json').read_bytes()).hexdigest(),
                 'aging_catalog':hashlib.sha256(catalog.read_bytes()).hexdigest(),
                 'inspect_templates.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
(HERE/'template-audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([{k:v for k,v in s.items() if k not in ['templates','overlap_with_published_aging_sample']} for s in sets]))

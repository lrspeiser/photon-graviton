from pathlib import Path
import csv
import hashlib
import io
import json
import tarfile
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
proto=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
manifest_path=HERE.parent/'snid-template-audit/manifest.json'
manifest=json.loads(manifest_path.read_text(encoding='utf-8'))
archive_path=ROOT/manifest['archives'][0]['cache']
assert hashlib.sha256(archive_path.read_bytes()).hexdigest()==manifest['archives'][0]['sha256']
catalog=HERE.parent/'electromagnetic-audit/spectral-aging-predictions.csv'
with catalog.open(newline='',encoding='utf-8') as f:
    target_names={'sn'+r['object'].lower() for r in csv.DictReader(f) if float(r['z'])<.1}
library=[];excluded=[];wave=None
with tarfile.open(archive_path) as archive:
    for member in archive.getmembers():
        if not member.name.endswith('.lnw'):continue
        lines=archive.extractfile(member).read().decode('ascii').splitlines()
        h=lines[0].split()
        if not h[7].startswith('Ia'):continue
        nk=int(h[4]);ages=np.array([float(v) for v in lines[nk+2].split()])
        if ages[0]!=0:continue
        data=np.loadtxt(io.StringIO('\n'.join(lines[nk+3:])))
        if wave is None:wave=data[:,0]
        assert np.array_equal(wave,data[:,0])
        for j,phase in enumerate(ages[1:]):
            if not proto['training_phase_days'][0]<=phase<=proto['training_phase_days'][1]:continue
            row=dict(object=h[5].lower(),column=j,phase=float(phase),file=member.name)
            values=data[:,j+1];support=np.flatnonzero(values!=0)
            if not len(support) or wave[support[0]]>4000 or wave[support[-1]]<7000:
                excluded.append(dict(row,reason='Insufficient full-window support'));continue
            library.append(dict(row,values=values))
assert len(target_names)==22
names=np.array([r['object'] for r in library]);phase=np.array([r['phase'] for r in library])
targets=[i for i,r in enumerate(library) if r['object'] in target_names and -10<=r['phase']<=30]
rows=[]
for lower,upper in proto['windows_angstrom']:
    mask=(wave>=lower)&(wave<=upper)
    flux=np.array([r['values'][mask] for r in library]);flux-=flux.mean(axis=1,keepdims=True)
    norms=np.linalg.norm(flux,axis=1)
    assert np.all(norms>0)
    normalized=flux/norms[:,None]
    scores=normalized[targets]@normalized.T
    for k,index in enumerate(targets):
        obj=names[index];candidates=[]
        for other in np.unique(names):
            if other==obj:continue
            ix=np.flatnonzero(names==other)
            best=ix[np.argmax(scores[k,ix])]
            candidates.append((float(scores[k,best]),int(best)))
        top=sorted(candidates,reverse=True)[:5]
        assert len(top)==5 and all(names[j]!=obj for _,j in top)
        rows.append(dict(object=obj,column=library[index]['column'],phase=float(phase[index]),
                         window=[lower,upper],estimated_phase=float(np.median([phase[j] for _,j in top])),
                         matches=[dict(object=str(names[j]),phase=float(phase[j]),score=score) for score,j in top]))
summaries=[]
for window in proto['windows_angstrom']:
    selected=[r for r in rows if r['window']==window]
    residual=np.array([r['estimated_phase']-r['phase'] for r in selected])
    objects=[]
    for obj in sorted(target_names):
        sample=[r for r in selected if r['object']==obj]
        if not sample:
            objects.append(dict(object=obj,n=0,status='No eligible spectra'));continue
        x=np.array([r['phase'] for r in sample]);y=np.array([r['estimated_phase'] for r in sample])
        item=dict(object=obj,n=len(sample),mean_error_days=float(np.mean(y-x)),
                  rms_error_days=float(np.sqrt(np.mean((y-x)**2))),span_days=float(np.ptp(x)))
        if len(sample)>=3 and np.ptp(x)>=15:
            item['slope']=float(np.linalg.lstsq(np.column_stack([np.ones(len(x)),x]),y,rcond=None)[0][1])
        objects.append(item)
    summaries.append(dict(window=window,n_spectra=len(selected),mean_error_days=float(residual.mean()),
                          median_absolute_error_days=float(np.median(abs(residual))),
                          rms_error_days=float(np.sqrt(np.mean(residual**2))),objects=objects))
out=dict(scope=proto['scope'],library_spectra=len(library),library_objects=len(set(names)),
         excluded=excluded,rows=rows,summaries=summaries,
         hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [catalog,manifest_path,HERE/'protocol.json',Path(__file__)]})
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(dict(library_spectra=len(library),targets=len(targets),summaries=[{k:v for k,v in s.items() if k!='objects'} for s in summaries])))

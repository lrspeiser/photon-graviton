from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
proto=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
loader=HERE.parent/'spectral-object-exclusion/run.py'
source=loader.read_text(encoding='utf-8')
assert source.count('rows=[]')==1
# Reuse only the verified data-loading portion, without rerunning its matcher.
ns={'__file__':str(loader),'__name__':'library_loader'}
exec(compile(source.split('rows=[]')[0],str(loader),'exec'),ns)
library,names,phase,targets,wave=[ns[k] for k in ['library','names','phase','targets','wave']]
flux=np.array([r['values'] for r in library]);target_flux=flux[targets]
target_names=names[targets];truth=phase[targets]
objects=np.unique(names);groups=[np.flatnonzero(names==obj) for obj in objects]
broad=(wave>=4000)&(wave<=7000)
noise_scale=np.std(target_flux[:,broad],axis=1)
assert np.all(noise_scale>0)
cache={}
for window in proto['windows_angstrom']:
    mask=(wave>=window[0])&(wave<=window[1]);basis=flux[:,mask].copy()
    basis-=basis.mean(axis=1,keepdims=True);basis/=np.linalg.norm(basis,axis=1)[:,None]
    cache[tuple(window)]=(mask,basis)

def estimate(values,window):
    mask,basis=cache[tuple(window)]
    x=values[:,mask].copy();x-=x.mean(axis=1,keepdims=True)
    x/=np.linalg.norm(x,axis=1)[:,None]
    score=x@basis.T
    object_score=np.empty((len(targets),len(objects)))
    object_index=np.empty(object_score.shape,dtype=int)
    for j,ix in enumerate(groups):
        best=ix[np.argmax(score[:,ix],axis=1)]
        object_index[:,j]=best;object_score[:,j]=score[np.arange(len(targets)),best]
    object_score[target_names[:,None]==objects[None,:]]=-np.inf
    top=np.argsort(object_score,axis=1)[:,-5:]
    selected=np.take_along_axis(object_index,top,axis=1)
    assert np.all(names[selected]!=target_names[:,None])
    return np.median(phase[selected],axis=1)

old_path=HERE.parent/'spectral-object-exclusion/results.json'
old=json.loads(old_path.read_text(encoding='utf-8'))
checks=[]
for window in proto['windows_angstrom']:
    actual=estimate(target_flux,window)
    expected={(r['object'],r['column']):r['estimated_phase'] for r in old['rows'] if r['window']==window}
    reference=np.array([expected[(library[i]['object'],library[i]['column'])] for i in targets])
    error=float(np.max(abs(actual-reference)));assert error<=proto['zero_noise_gate_days']
    checks.append(dict(window=window,max_zero_noise_difference_days=error))
rows=[]
for seed in proto['seeds']:
    draw=np.random.default_rng(seed).normal(size=target_flux.shape)
    for snr in proto['contrast_snr']:
        noisy=target_flux+draw*noise_scale[:,None]/snr
        for window in proto['windows_angstrom']:
            recovered=estimate(noisy,window);error=recovered-truth
            per_object=[]
            for obj in sorted(set(target_names)):
                ix=np.flatnonzero(target_names==obj);x=truth[ix];y=recovered[ix]
                item=dict(object=str(obj),n=len(ix),mean_error_days=float(np.mean(y-x)))
                if len(ix)>=3 and np.ptp(x)>=15:
                    item['slope']=float(np.linalg.lstsq(np.column_stack([np.ones(len(ix)),x]),y,rcond=None)[0][1])
                per_object.append(item)
            slopes=[r['slope'] for r in per_object if 'slope' in r]
            rows.append(dict(seed=seed,contrast_snr=snr,window=window,rms_error_days=float(np.sqrt(np.mean(error**2))),
                             mean_error_days=float(np.mean(error)),median_object_slope=float(np.median(slopes)),
                             object_summaries=per_object,recovered_phase_days=recovered.tolist()))
    print(json.dumps(dict(seed_completed=seed)),flush=True)
summaries=[]
for snr in proto['contrast_snr']:
    for window in proto['windows_angstrom']:
        selected=[r for r in rows if r['contrast_snr']==snr and r['window']==window]
        v=np.array([r['median_object_slope'] for r in selected])
        summaries.append(dict(contrast_snr=snr,window=window,mean_rms_error_days=float(np.mean([r['rms_error_days'] for r in selected])),
                              mean_seed_median_slope=float(v.mean()),seed_median_slope_range=[float(v.min()),float(v.max())]))
out=dict(scope=proto['scope'],checks=checks,summaries=summaries,rows=rows,
         targets=[dict(object=library[i]['object'],column=library[i]['column'],phase=library[i]['phase']) for i in targets],
         hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [loader,old_path,HERE/'protocol.json',Path(__file__)]})
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summaries))

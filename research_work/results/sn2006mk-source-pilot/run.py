from pathlib import Path
import json,hashlib
import numpy as np
from scipy.interpolate import RegularGridInterpolator
H=Path(__file__).resolve().parent;R=H.parents[2]
sp=H.parent/'spectral-real-pair/results.json';ph=H.parent/'sn2006mk-joint-photometry/results.json';source=R/'research_work/generated/des-photometric-calibration/selected/snsed/Hsiao07.dat'
spec=json.loads(sp.read_text(encoding='utf-8'));phot=json.loads(ph.read_text(encoding='utf-8'));a=np.loadtxt(source);times=np.unique(a[:,0]);wave=np.unique(a[:,1]);flux=a[:,2].reshape(len(times),len(wave));interp=RegularGridInterpolator((times,wave),flux,bounds_error=True)
assert hashlib.sha256(source.read_bytes()).hexdigest()==json.loads((H.parent/'hsiao-source-audit/results.json').read_text())['cached_sha256']
z=.4754;S=1+z;dt=spec['observer_span_days'];primary=sorted([x for x in spec['real_epoch_fits'] if x['degree']==2 and not x['host_basis']],key=lambda x:x['epoch']);early=primary[0]['estimated_phase']
cases={'spectral_pair':(early,primary[1]['estimated_phase']),'event_stretch_1_plus_z':(early,early+dt/S),'unchanged_arrival_intervals':(early,early+dt)}
rows=[];hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [sp,ph,source,Path(__file__)]}
for band in ['R','I']:
 p=R/'research_work/generated/essence-passbands'/f'CTIO4m_{band}.dat';hashes[str(p.relative_to(R))]=hashlib.sha256(p.read_bytes()).hexdigest();data=np.loadtxt(p);w,energy=data.T
 # Legacy energy response already includes the photon-counting wavelength factor.
 # For a physical shifted response: move photon transmission, then multiply lambda once.
 photon=energy/w
 observed=next(x for x in phot['spectral_links'] if x['band']==band+'4m')
 for shift in [-50.,0.,50.]:
  wo=w+shift;response=photon*wo
  for clip in [False,True]:
   response_used=np.maximum(response,0) if clip else response
   for label,(t1,t2) in cases.items():
    values=[]
    for time in [t1,t2]:
     sed=interp(np.column_stack([np.full(len(wo),time),wo/S]));values.append(float(np.trapezoid(sed*response_used,wo)))
    ratio=values[1]/values[0];rows.append(dict(band=band,case=label,source_phases=[t1,t2],filter_shift_angstrom=shift,clip_negative_response=clip,predicted_flux_ratio=ratio,observed_flux_ratio=observed['later_to_earlier_flux_ratio'],observed_diagonal_sigma=observed['diagonal_linearized_ratio_sigma']))
# Independent interpolation at one exact tabulated phase reproduces its row.
check=interp(np.column_stack([np.zeros(len(wave)),wave]));assert np.array_equal(check,flux[np.flatnonzero(times==0)[0]])
out=dict(scope='Conditional source/legacy-filter pilot with no photometric fit; not release-calibrated validation',source_model='Empirical Hsiao author template',cases=cases,rows=rows,hashes=hashes,limitations=['Legacy filters are not release matched','50-Angstrom shifts are illustrative, not a calibrated error bound','Spectral ages and empirical source model have unresolved shared calibration assumptions','No dust or intrinsic source diversity fitted','Close photometric/spectral epoch offsets are not corrected','Flux ratios cancel constant propagation amplitude; absolute brightness not tested'])
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([x for x in rows if x['filter_shift_angstrom']==0 and not x['clip_negative_response']]))

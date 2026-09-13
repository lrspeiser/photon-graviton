from pathlib import Path
import numpy as np,json,hashlib,tarfile
H=Path(__file__).resolve().parent;R=H.parents[2];C=R/'research_work/generated/essence-passbands'
C.mkdir(parents=True,exist_ok=True)
needed=['README','CTIO4m_R.dat','CTIO4m_I.dat']
if any(not (C/name).exists() for name in needed):
 with tarfile.open(R/'research_work/generated/des-photometric-calibration/SNDATA_ROOT_2024-07-03.tar.gz') as archive:
  for name in needed:
   (C/name).write_bytes(archive.extractfile('filters/ESSENCE/'+name).read())
if (H/'results.json').exists():
 previous=json.loads((H/'results.json').read_text(encoding='utf-8'))
 for row in previous['rows']:assert hashlib.sha256((C/row['filename']).read_bytes()).hexdigest()==row['sha256']
 assert hashlib.sha256((C/'README').read_bytes()).hexdigest()==previous['readme_sha256']
# Printed positive-throughput samples of Narayan et al. 2016 Table B2.
anchors={'R':[(6240,.5338),(7005,.3378),(7775,.1092),(8540,.0185)],'I':[(7635,.3580),(8330,.4729),(9030,.0023),(9725,.0008)]}
rows=[]
for band,points in anchors.items():
 p=C/f'CTIO4m_{band}.dat';d=np.loadtxt(p);w,e=d.T;assert np.all(np.diff(w)>0)
 photon=e/w # README states energy convention has one additional lambda.
 x,y=np.array(points).T;old=np.interp(x,w,photon,left=0,right=0);a=float(old@y/(old@old));pred=a*old
 # Retain negative interpolation lobes; do not silently clip stored inputs.
 rows.append(dict(band=band,filename=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),range_angstrom=[float(w[0]),float(w[-1])],negative_samples=int(sum(e<0)),samples=len(w),best_amplitude_to_printed_samples=a,anchors=[dict(wavelength=float(xx),published_total=float(yy),rescaled_old=float(pp),relative_difference=float(pp/yy-1)) for xx,yy,pp in zip(x,y,pred)],rms_absolute_anchor_difference=float(np.sqrt(np.mean((pred-y)**2)))))
out=dict(scope='Audit of locally available legacy response; not adoption as the 2016 photon throughput',archive='research_work/generated/des-photometric-calibration/SNDATA_ROOT_2024-07-03.tar.gz',archive_members=['filters/ESSENCE/README','filters/ESSENCE/CTIO4m_R.dat','filters/ESSENCE/CTIO4m_I.dat'],readme=(C/'README').read_text(),readme_sha256=hashlib.sha256((C/'README').read_bytes()).hexdigest(),table_source='https://lss.fnal.gov/archive/2016/pub/fermilab-pub-16-402-ae.pdf',table='B2, printed selected samples only',rows=rows)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(rows))

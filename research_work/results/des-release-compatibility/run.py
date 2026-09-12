"""Compare paired DES release measurement arrays without fitting observations."""
from pathlib import Path
import hashlib,json,tarfile
import numpy as np
from astropy.io import fits
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/generated/des-release-compatibility';CACHE.mkdir(exist_ok=True)
a=json.loads((HERE.parent/'des-photometric-calibration/acquisition.json').read_text(encoding='utf-8'))
archive=ROOT/'research_work/generated/des-photometric-calibration/SNDATA_ROOT_2024-07-03.tar.gz'
hash=hashlib.sha256()
with archive.open('rb') as f:
    while block:=f.read(4*1024*1024):hash.update(block)
assert hash.hexdigest()==a['archive_sha256']
prefix='lcmerge/DES-SN5YR/DES-SN5YR_DES/'
names=[prefix+'DES-SN5YR_DES'+s for s in ['.README','_HEAD.FITS.gz','_PHOT.FITS.gz']]
manifest=[]
with tarfile.open(archive,'r|gz') as tar:
    for member in tar:
        if not member.isfile() or member.name not in names:continue
        data=tar.extractfile(member).read();target=CACHE/Path(member.name).name;target.write_bytes(data)
        manifest.append({'member':member.name,'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)})
assert len(manifest)==3
current=json.loads((HERE.parent/'timing-foundation/input-audit.json').read_text(encoding='utf-8'))
paths={Path(v['path']).name:ROOT/v['cache_relative_path'] for v in current['input_manifest']}
for v in current['input_manifest']:
    if v['path'].endswith('FITS.gz'):assert hashlib.sha256((ROOT/v['cache_relative_path']).read_bytes()).hexdigest()==v['sha256']
oldh=fits.getdata(CACHE/'DES-SN5YR_DES_HEAD.FITS.gz',1);oldp=fits.getdata(CACHE/'DES-SN5YR_DES_PHOT.FITS.gz',1)
newh=fits.getdata(paths['DES-SN5YR_DES_HEAD.FITS.gz'],1);newp=fits.getdata(paths['DES-SN5YR_DES_PHOT.FITS.gz'],1)
old={str(v['SNID']).strip():v for v in oldh};new={str(v['SNID']).strip():v for v in newh}
assert len(old)==len(oldh) and len(new)==len(newh)
common=sorted(set(old)&set(new));bad_order=[];different_count=[];total=0
fields=[n for n in oldp.names if n in newp.names and n not in ['MJD','BAND']]
stats={n:{'different_entries':0,'max_absolute_difference':0.} for n in fields}
headstats={n:0 for n in oldh.names if n in newh.names and n not in ['SNID','PTROBS_MIN','PTROBS_MAX']}
def band(x):return np.char.replace(np.char.strip(x.astype(str)),'DES-','')
for snid in common:
    h0,h1=old[snid],new[snid];p0=oldp[int(h0['PTROBS_MIN'])-1:int(h0['PTROBS_MAX'])];p1=newp[int(h1['PTROBS_MIN'])-1:int(h1['PTROBS_MAX'])]
    for name in headstats:
        v0=np.asarray(h0[name]);v1=np.asarray(h1[name])
        eq=np.array_equal(v0,v1,equal_nan=True) if v0.dtype.kind in 'fc' and v1.dtype.kind in 'fc' else np.array_equal(v0,v1)
        if not eq:headstats[name]+=1
    if len(p0)!=len(p1):different_count.append(snid);continue
    if not np.array_equal(p0['MJD'],p1['MJD']) or not np.array_equal(band(p0['BAND']),band(p1['BAND'])):bad_order.append(snid);continue
    total+=len(p0)
    for name in fields:
        v0=np.asarray(p0[name]);v1=np.asarray(p1[name]);same=v0==v1
        if v0.dtype.kind in 'fc' and v1.dtype.kind in 'fc':same|=np.isnan(v0)&np.isnan(v1)
        stats[name]['different_entries']+=int(np.count_nonzero(~same))
        if v0.dtype.kind in 'fciu' and v1.dtype.kind in 'fciu':
            delta=abs(v0.astype(float)-v1.astype(float));finite=delta[np.isfinite(delta)]
            stats[name]['max_absolute_difference']=max(stats[name]['max_absolute_difference'],float(finite.max(initial=0)))
result={'scope':'Version comparison of already exposed DES source data; no light-curve fit or physical model inference.','archive_sha256':a['archive_sha256'],'old_files':manifest,'new_source_commit':current['input_manifest'][0]['commit'],'old_events':len(old),'new_events':len(new),'matched_events':len(common),'old_only_count':len(set(old)-set(new)),'new_only_count':len(set(new)-set(old)),'different_point_count_ids':different_count,'unaligned_observation_ids':bad_order,'aligned_measurement_rows':total,'measurement_field_comparisons':stats,'changed_header_fields':{k:v for k,v in headstats.items() if v},'old_only_header_fields':sorted(set(oldh.names)-set(newh.names)),'new_only_header_fields':sorted(set(newh.names)-set(oldh.names)),'old_only_photometry_fields':sorted(set(oldp.names)-set(newp.names)),'new_only_photometry_fields':sorted(set(newp.names)-set(oldp.names))}
(HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['old_files','scope','archive_sha256']}))

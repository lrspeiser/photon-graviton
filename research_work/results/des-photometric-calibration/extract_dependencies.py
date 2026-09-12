"""Extract the exact dependencies named by the verified calibration input."""
from pathlib import Path
import hashlib,json,tarfile
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CACHE=ROOT/'research_work/generated/des-photometric-calibration'
a=json.loads((HERE/'acquisition.json').read_text(encoding='utf-8'))
archive=CACHE/'SNDATA_ROOT_2024-07-03.tar.gz'
sha=hashlib.sha256()
with archive.open('rb') as stream:
    while block:=stream.read(4*1024*1024):sha.update(block)
assert sha.hexdigest()==a['archive_sha256']
names=['filters/DES/DES-SN3YR_DECam/DECam_'+b+'.dat' for b in 'griz']+['snsed/Hsiao07.dat']
selected=[];out=(CACHE/'selected').resolve()
with tarfile.open(archive,'r|gz') as tar:
    for member in tar:
        if not member.isfile() or member.name not in names:continue
        path=(out/member.name).resolve();assert path.is_relative_to(out)
        data=tar.extractfile(member).read();path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
        selected.append({'archive_path':member.name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'cache_path':path.relative_to(ROOT).as_posix()})
assert {v['archive_path'] for v in selected}==set(names)
(HERE/'dependencies.json').write_text(json.dumps({'archive_sha256':a['archive_sha256'],'selected':selected,'template_use':'Hsiao template acquired as a declared input dependency, not adopted as an independent timing or luminosity prior.'},indent=2)+'\n',encoding='utf-8',newline='\n')
print('Extracted '+str(len(selected))+' dependencies from verified archive.')

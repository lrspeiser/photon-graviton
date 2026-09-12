"""Acquire checksum-pinned original DES calibration resources; execute no external code."""
from pathlib import Path
import hashlib,json,tarfile,urllib.request
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/generated/des-photometric-calibration'
URL='https://zenodo.org/api/records/12655677/files/SNDATA_ROOT_2024-07-03.tar.gz/content'
NAME='SNDATA_ROOT_2024-07-03.tar.gz';SIZE=1653076974;MD5='b75e8fec76e53b9b9e51901a09226b46'
def main():
    CACHE.mkdir(parents=True,exist_ok=True);archive=CACHE/NAME
    if not archive.exists():
        temporary=CACHE/(NAME+'.partial');n=0;milestone=100_000_000
        with urllib.request.urlopen(URL,timeout=60) as response,temporary.open('wb') as target:
            while block:=response.read(4*1024*1024):
                target.write(block);n+=len(block)
                if n>=milestone:print(f'Downloaded {n}/{SIZE} bytes',flush=True);milestone+=100_000_000
        if n!=SIZE:raise RuntimeError('Archive size mismatch')
        temporary.replace(archive)
    md5=hashlib.md5();sha=hashlib.sha256()
    with archive.open('rb') as stream:
        while block:=stream.read(4*1024*1024):md5.update(block);sha.update(block)
    assert archive.stat().st_size==SIZE and md5.hexdigest()==MD5
    print('Full archive checksum verified.',flush=True)
    selected=[];names=[];out=(CACHE/'selected').resolve()
    with tarfile.open(archive,'r|gz') as tar:
        for member in tar:
            name=member.name;names.append(name)
            key='/'+name.lstrip('./')
            if not member.isfile() or not ('/kcor/DES/' in key or '/standards/' in key):continue
            path=(out/name).resolve()
            if not path.is_relative_to(out):raise RuntimeError('Unsafe archive path')
            data=tar.extractfile(member).read();path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
            selected.append({'archive_path':name,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),'cache_path':str(path.relative_to(ROOT)).replace('\\','/')})
    result={'source_doi':'10.5281/zenodo.12655677','archive_url':URL,'archive_size':SIZE,'archive_md5':MD5,'archive_sha256':sha.hexdigest(),'archive_member_count':len(names),'selected':selected,'policy':'Only regular files under kcor/DES or standards extracted. No external code executed or transient-flux products inspected. This older archive requires compatibility audit with the pinned later DES metadata.'}
    (HERE/'acquisition.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    (CACHE/'archive-members.json').write_text(json.dumps(names,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'selected_count':len(selected),'target_calibration':[v['cache_path'] for v in selected if v['archive_path'].endswith('calib_DES-SN5YR_DES.input')]}),flush=True)
if __name__=='__main__':main()

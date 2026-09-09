"""Pinned public-input acquisition only; no timing fit and no downloaded code execution."""
from pathlib import Path
import csv,gzip,hashlib,json,urllib.request
from astropy.io import fits
import astropy
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
CACHE=ROOT/'research_work/generated/timing-inputs'

def main():
    protocol=json.loads((HERE/'acquisition-protocol.json').read_text(encoding='utf-8'))
    CACHE.mkdir(parents=True,exist_ok=True);records=[]
    for repo,commit,paths,group in [(protocol['release_repository'],protocol['release_commit'],protocol['inputs'],'release'),
                                    (protocol['method_repository'],protocol['method_commit'],protocol['method_inputs'],'methods')]:
        tree=json.load(urllib.request.urlopen(f'https://api.github.com/repos/{repo}/git/trees/{commit}?recursive=1'))
        assert not tree.get('truncated',False)
        lookup={x['path']:x for x in tree['tree'] if x['type']=='blob'}
        for path in paths:
            meta=lookup[path];dest=CACHE/group/path;dest.parent.mkdir(parents=True,exist_ok=True)
            url=f'https://raw.githubusercontent.com/{repo}/{commit}/{path}'
            if not dest.exists():
                with urllib.request.urlopen(url,timeout=120) as response:
                    payload=response.read()
                dest.write_bytes(payload)
            payload=dest.read_bytes()
            blob=hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()
            assert blob==meta['sha'] and len(payload)==meta['size'],path
            records.append({'repository':repo,'commit':commit,'path':path,'url':url,'bytes':len(payload),
                            'git_blob_sha1':blob,'sha256':hashlib.sha256(payload).hexdigest(),
                            'cache_relative_path':dest.relative_to(ROOT).as_posix()})
            print('Verified',path,len(payload),flush=True)
    base=CACHE/'release/0_DATA/DES-SN5YR_DES'
    with fits.open(base/'DES-SN5YR_DES_HEAD.FITS.gz',memmap=False) as h, fits.open(base/'DES-SN5YR_DES_PHOT.FITS.gz',memmap=False) as p:
        hd=h[1].data;pd=p[1].data
        start=np.asarray(hd['PTROBS_MIN'],dtype=np.int64);end=np.asarray(hd['PTROBS_MAX'],dtype=np.int64)
        assert np.all(start>=1) and np.all(end>=start) and np.all(end<=len(pd))
        lengths=end-start+1
        assert np.all(start[1:]>end[:-1])
        assert np.all(lengths==hd['NOBS'])
        assert len(set(hd['SNID']))==len(hd)
        physical=np.zeros(len(pd),dtype=bool)
        for lo,hi in zip(start,end):physical[lo-1:hi]=True
        times=np.asarray(pd['MJD'])[physical];flux=np.asarray(pd['FLUXCAL'])[physical];err=np.asarray(pd['FLUXCALERR'])[physical]
        band_name='BAND' if 'BAND' in pd.names else 'FLT'
        bands,counts=np.unique(pd[band_name][physical],return_counts=True)
        fit_audit={'head_rows':len(hd),'photometry_table_rows_including_separators':len(pd),
                   'pointer_selected_photometry_rows':int(np.sum(physical)),
                   'head_columns':list(hd.names),'photometry_columns':list(pd.names),
                   'min_observations_per_head_row':int(np.min(lengths)),'max_observations_per_head_row':int(np.max(lengths)),
                   'observer_MJD_range':[float(np.min(times)),float(np.max(times))],
                   'nonfinite_MJD_rows':int(np.sum(~np.isfinite(times))),
                   'nonfinite_flux_or_error_rows':int(np.sum(~np.isfinite(flux)|~np.isfinite(err))),
                   'nonpositive_error_rows':int(np.sum(err<=0)),
                   'bands':{str(b).strip():int(n) for b,n in zip(bands,counts)},
                   'pointer_validation':'1-based inclusive nonoverlapping ranges inside photometry table; separator rows excluded by pointers',
                   'timing_transform_applied':False,'sample_quality_cuts_applied':False}
        types,type_counts=np.unique(hd['SNTYPE'],return_counts=True)
        flags,flag_counts=np.unique(hd['FAKE'],return_counts=True)
        iz=np.asarray(hd['REDSHIFT_HELIO'])[hd['SNTYPE']==1]
        fit_audit['survey_type_counts']={str(int(t)):int(n) for t,n in zip(types,type_counts)}
        fit_audit['FAKE_flag_counts']={str(int(t)):int(n) for t,n in zip(flags,flag_counts)}
        fit_audit['spectroscopic_Ia_inventory']={'count':len(iz),'positive_finite_heliocentric_redshifts':int(np.sum(np.isfinite(iz)&(iz>0))),
                                                'heliocentric_redshift_range':[float(np.min(iz)),float(np.max(iz))],
                                                'scope':'Inventory only; no analysis sample or quality cuts selected.'}
    with (CACHE/'release/3_CLASSIFICATION/DES_classification.csv').open(encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f);names=reader.fieldnames;count=sum(1 for _ in reader)
    audit={'scope':protocol['scope'],'exposure':protocol['exposure'],'input_manifest':records,
           'protocol_sha256':hashlib.sha256((HERE/'acquisition-protocol.json').read_bytes()).hexdigest(),
           'fits_audit':fit_audit,'classification_audit':{'columns':names,'rows':count},
           'inference_performed':False,'runtime':{'numpy':np.__version__,'astropy':astropy.__version__},
           'missing_before_inference':['Frozen inference/selection protocol','Observer-time peak/normalization fit or justified audit of SALT-derived alternatives','Filter throughput/wavelength matching','Injection-recovery tests covering b=0 and b=1','Event/shared-reference covariance treatment','Explicit source-duration evolution and selection sensitivities']}
    (HERE/'input-audit.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(fit_audit,indent=2),flush=True)

if __name__=='__main__':main()

"""Acquire Gaia DR3 covariance and quality fields for the prepared sample only."""
from pathlib import Path
import json
import pandas as pd
from astropy.table import Table
from astroquery.gaia import Gaia
from prepare import CACHE, HERE, digest

QUERY = '''SELECT g.source_id, g.ra, g.dec, g.ra_error, g.dec_error,
 g.parallax, g.parallax_error, g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
 g.ra_dec_corr, g.ra_parallax_corr, g.ra_pmra_corr, g.ra_pmdec_corr,
 g.dec_parallax_corr, g.dec_pmra_corr, g.dec_pmdec_corr,
 g.parallax_pmra_corr, g.parallax_pmdec_corr, g.pmra_pmdec_corr,
 g.ruwe, g.visibility_periods_used, g.astrometric_params_solved,
 g.duplicated_source, g.phot_g_mean_mag, g.bp_rp,
 g.radial_velocity, g.radial_velocity_error, g.non_single_star
 FROM gaiadr3.gaia_source AS g
 JOIN TAP_UPLOAD.targets AS t ON g.source_id=t.source_id'''

if __name__=='__main__':
    ids=pd.read_parquet(CACHE/'matched-exploratory.parquet',columns=['source_id'])
    assert ids.source_id.dtype.kind in 'iu' and ids.source_id.is_unique
    (HERE/'gaia-query.adql').write_text(QUERY+'\n',encoding='utf-8',newline='\n')
    upload=CACHE/'gaia-targets.xml'
    Table.from_pandas(ids).write(upload,format='votable',overwrite=True)
    job=Gaia.launch_job_async(QUERY,upload_resource=str(upload),upload_table_name='targets',
                              background=True)
    print('Job:',job.jobid,flush=True)
    (HERE/'gaia-job.json').write_text(json.dumps({'job_id':str(job.jobid),'requested_sources':len(ids)},indent=2)+'\n',encoding='utf-8',newline='\n')
    result=job.get_results()
    assert len(result)==len(ids), 'Archive returned fewer rows than requested'
    assert len(set(result['source_id']))==len(ids)
    p=CACHE/'gaia-dr3-quality-covariance.fits'
    result.write(p,overwrite=True)
    manifest={'rows':len(result),'sha256':digest(p),'query':'gaia-query.adql',
              'job_id':str(job.jobid),'archive':'https://gea.esac.esa.int/tap-server/tap'}
    (HERE/'gaia-download.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(manifest,flush=True)

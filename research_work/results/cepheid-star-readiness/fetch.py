"""Fetch the full Gaia DR3 DCEP parent join, preserving rejected candidates."""
from pathlib import Path
import csv
import io
import json
import hashlib
import requests
from astropy.table import Table
from astroquery.gaia import Gaia

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
CACHE=ROOT/'research_work/data-cache/cepheid-stars';CACHE.mkdir(parents=True,exist_ok=True)
QUERY='''SELECT g.source_id, g.ref_epoch, g.ra, g.dec, g.l, g.b,
g.ra_error, g.dec_error, g.parallax, g.parallax_error,
g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
g.ra_dec_corr, g.ra_parallax_corr, g.ra_pmra_corr, g.ra_pmdec_corr,
g.dec_parallax_corr, g.dec_pmra_corr, g.dec_pmdec_corr,
g.parallax_pmra_corr, g.parallax_pmdec_corr, g.pmra_pmdec_corr,
g.ruwe, g.visibility_periods_used, g.astrometric_params_solved,
g.duplicated_source, g.non_single_star,
g.radial_velocity, g.radial_velocity_error, g.rv_nb_transits,
g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
c.type_best_classification, c.mode_best_classification, c.multi_mode_best_classification,
c.pf, c.pf_error, c.p1_o, c.p1_o_error, c.p2_o, c.p2_o_error,
c.int_average_g, c.int_average_g_error, c.int_average_bp, c.int_average_bp_error,
c.int_average_rp, c.int_average_rp_error, c.average_rv, c.average_rv_error,
c.num_clean_epochs_rv, c.num_clean_epochs_g, c.peak_to_peak_rv,
c.metallicity, c.metallicity_error
FROM gaiadr3.vari_cepheid AS c
JOIN gaiadr3.gaia_source AS g ON c.source_id=g.source_id
WHERE c.type_best_classification='DCEP' '''
COUNT="SELECT COUNT(*) AS n FROM gaiadr3.vari_cepheid AS c JOIN gaiadr3.gaia_source AS g ON c.source_id=g.source_id WHERE c.type_best_classification='DCEP'"
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
(HERE/'query.adql').write_text(QUERY.strip()+'\n',newline='\n')
query_hash=hashlib.sha256(QUERY.encode()).hexdigest()
path=CACHE/'gaia-dr3-dcep-parent.votable'
state=HERE/'job.json'
if not path.exists():
    if state.exists():
        info=json.loads(state.read_text());assert info['query_sha256']==query_hash
        job=Gaia.load_async_job(jobid=info['job_id'])
    else:
        job=Gaia.launch_job_async(QUERY,background=True)
        state.write_text(json.dumps(dict(job_id=str(job.jobid),query_sha256=query_hash),indent=2)+'\n',newline='\n')
    print('Archive job',job.jobid,flush=True)
    table=job.get_results();table.write(path,format='votable',overwrite=False)
else:
    table=Table.read(path,format='votable')
response=requests.get('https://gea.esac.esa.int/tap-server/tap/sync',params={
    'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':COUNT},timeout=45)
response.raise_for_status()
expected=int(next(csv.DictReader(io.StringIO(response.text)))['n'])
assert len(table)==expected
id_column=next(name for name in table.colnames if name.lower()=='source_id')
assert table[id_column].dtype.kind in 'iu' and len(set(table[id_column]))==len(table)
manifest=dict(rows=len(table),independent_archive_count=expected,columns=len(table.colnames),
    path=str(path),sha256=digest(path),query_sha256=query_hash,
    job=json.loads(state.read_text()),archive='https://gea.esac.esa.int/tap-server/tap',
    no_top_limit_or_quality_filter=True,unique_integer_ids=True,archive_id_column=id_column,
    role='Parent catalog for a future common-frame analysis; not the reconstructed 903-star sample')
(HERE/'download.json').write_text(json.dumps(manifest,indent=2)+'\n',newline='\n')
print(json.dumps(manifest,indent=2),flush=True)

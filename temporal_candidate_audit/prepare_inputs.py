from pathlib import Path
import shutil,json,urllib.request,hashlib
from concurrent.futures import ThreadPoolExecutor
P=Path(__file__).resolve().parent; root=P.parent; D=P/'data'
files={'redshift_features.json':'temporal_hypothesis_checks/source_inputs/features.json','train_labels.json':'temporal_hypothesis_checks/source_inputs/train_labels.json','validation_labels.json':'temporal_hypothesis_checks/source_inputs/validation_labels.json','test_labels.json':'temporal_hypothesis_checks/source_inputs/test_labels.json','cf4_table2.dat':'redshift_paper_sources/expanded/time-redshift-expanded/table2.dat','cf4_ReadMe.txt':'redshift_paper_sources/expanded/time-redshift-expanded/ReadMe','spectral_aging.json':'time_first_principles/data/spectral_aging.json','firas.txt':'time_first_principles/data/firas.txt','DES_event_averages.csv':'time_revision/data/DES_event_averages.csv','SPARC_Lelli2016c.mrt':'option3_test/raw/SPARC_Lelli2016c.mrt','Rotmod_LTG.zip':'option3_test/raw/Rotmod_LTG.zip','sparc_frozen.json':'option3_test/frozen.json','prior_cliff_results.json':'option3_cliff/results.json','prior_cliff_models.json':'option3_cliff/fold_models.json','prior_redshift_scores.json':'redshift_paper_sources/expanded/time-redshift-expanded/scores.json'}
for dst,src in files.items():shutil.copy2(root/src,D/dst)
urls={f'planck_{kind}.txt':f'https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_PowerSpect_CMB-{kind}-binned_R3.02.txt' if kind=='EE' else f'https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_PowerSpect_CMB-{kind}-binned_R3.01.txt' for kind in ['TT','EE']}
def fetch(item):
 name,url=item
 try:
  b=urllib.request.urlopen(url,timeout=40).read()
  if not b.startswith(b'#'):raise ValueError('Unexpected non-table content')
  (D/name).write_bytes(b);return {'file':name,'url':url,'bytes':len(b),'status':'downloaded'}
 except Exception as e:return {'file':name,'url':url,'status':'failed','error':str(e)}
with ThreadPoolExecutor(max_workers=2) as ex:download=list(ex.map(fetch,urls.items()))
external={'GW170817':{'gamma_minus_GW_seconds':1.74,'sigma_seconds':.05,'nominal_distance_Mpc':40.,'distance_sensitivity_Mpc':[30.,40.,50.],'source':'https://arxiv.org/abs/1710.05834','note':'Published timing summary; distances are illustrative sensitivity values, not a refit of a standard siren in modified gravity.'},'Galileo':{'epsilonG':.19e-5,'sigma':2.48e-5,'source':'https://doi.org/10.1103/PhysRevLett.121.231101'},'alpha_clock':{'alpha_dot_per_year':1e-18,'sigma_per_year':1.1e-18,'source':'https://doi.org/10.1103/PhysRevLett.126.011102'},'HFLS3':{'z':6.34,'T_interval_1sigma_K':[16.4,30.2],'source':'https://arxiv.org/abs/2202.00693'},'T0':{'K':2.72548,'sigma_K':.00057,'source':'https://arxiv.org/abs/0911.1955'}}
(D/'published_measurements.json').write_text(json.dumps(external,indent=2))
manifest={'downloads':download,'copied_inputs':files,'sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in D.iterdir() if f.is_file()}}
(P/'input_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(download,indent=2))

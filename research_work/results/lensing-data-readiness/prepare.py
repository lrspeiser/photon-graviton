"""SLACS observation/inference separation and frozen candidate roles."""
from pathlib import Path
import json,hashlib,shutil
import numpy as np
from astropy.table import Table
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
C=ROOT/'research_work/data-cache/lensing-observations'
H.mkdir(exist_ok=True,parents=True)
obs=Table.read(C/'table4.dat',format='ascii.cds',readme=str(C/'ReadMe')).to_pandas()
lens=Table.read(C/'table5.dat',format='ascii.cds',readme=str(C/'ReadMe')).to_pandas()
assert len(obs)==131 and len(lens)==63 and obs.Name.is_unique and lens.Name.is_unique
d=lens.merge(obs,on='Name',validate='one_to_one')
assert len(d)==63 and d.Lens.eq('A').all()
assert (d.zBG>d.zFG).all() and (d.zFG>0).all()
salt='photon-graviton-slacs-object-v1'
def role(name):
    bucket=int.from_bytes(hashlib.sha256((salt+'|'+name).encode()).digest()[:8],'big')%10
    return 'training' if bucket<6 else ('validation' if bucket<8 else 'test')
d['role']=d.Name.map(role)
cols=['Name','SDSS','Plate','MJD','Fiber','zFG','zBG','Imag','AI','Re','b/a','PA','sigma','e_sigma','Mph','Mul','Lens','n_Lens','bSIE','qSIE','PASIE','bLTM','gammaLTM','PAgLTM','Nsrc','Ring?','Good?','role']
clean=d[cols].copy()
clean['spectroscopic_dispersion_available']=d.sigma.notna()&d.e_sigma.notna()&(d.sigma>0)&(d.e_sigma>0)
# This is conditional model geometry, not an independent distance catalog.
fit=ROOT/'research_work/results/joint-galaxy-audit/results.json'
alpha=json.loads(fit.read_text())['redshift']['parameters']['alpha_per_Mpc']
geometry=[]
for r in d.itertuples():
    Dl=np.log1p(r.zFG)/alpha;Ds=np.log1p(r.zBG)/alpha
    geometry.append(dict(Name=r.Name,role=r.role,conditional_Dl_Mpc=Dl,conditional_Ds_Mpc=Ds,
        conditional_Dls_over_Ds=1-Dl/Ds,geometry='Static Euclidean distances inferred using frozen exponential redshift rule; peculiar and endpoint redshifts omitted.'))
# JSON null preserves missing values. No imputation of unavailable dispersions.
(H/'lens-observations-and-image-models.json').write_text(clean.to_json(orient='records',indent=2)+'\n',encoding='utf-8')
(H/'conditional-geometry.json').write_text(json.dumps(geometry,indent=2)+'\n',encoding='utf-8')
summary=dict(candidates=len(obs),modeled_lenses=len(d),spectroscopic_dispersion_available=int(clean.spectroscopic_dispersion_available.sum()),
    published_good_dispersion_flag={str(k):int(v) for k,v in d['Good?'].value_counts().items()},
    morphology={str(k):int(v) for k,v in d.Mph.value_counts().items()},roles={str(k):int(v) for k,v in d.role.value_counts().items()},
    roles_with_available_dispersion={str(k):int(v) for k,v in clean[clean.spectroscopic_dispersion_available].role.value_counts().items()},
    redshift_ranges={'lens':[float(d.zFG.min()),float(d.zFG.max())],'source':[float(d.zBG.min()),float(d.zBG.max())]},
    frozen_alpha_per_Mpc=alpha,split_salt=salt,split_method='First eight SHA256 bytes, big endian, modulo ten; 0–5 train, 6–7 validation, 8–9 test. Split by system ID.',
    provenance='Publicly available 2008 sample; basic catalog fields and aggregate summaries inspected. Assigned roles do not imply pristine blind validation or an independent galaxy population.',
    formulas_fitted=False,observed_lensing_scores_evaluated=False,
    raw_sha256={n:hashlib.sha256((C/n).read_bytes()).hexdigest() for n in ['ReadMe','table4.dat','table5.dat','notes.dat']},
    origin='https://cdsarc.cds.unistra.fr/ftp/J/ApJ/682/964/',retrieved='2026-09-10')
(H/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
raw=H/'source-tables';raw.mkdir(exist_ok=True)
for n in ['ReadMe','table4.dat','table5.dat','notes.dat']:shutil.copyfile(C/n,raw/n)
assert np.all(np.array([r['conditional_Dls_over_Ds'] for r in geometry])>0)
print(json.dumps(summary,indent=2))

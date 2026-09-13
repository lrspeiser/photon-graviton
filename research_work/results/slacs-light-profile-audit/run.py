from pathlib import Path
import json,csv,hashlib
import numpy as np
from scipy.special import gammainc,gamma
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
source=ROOT/'research_work/generated/slacs-resolved-data/ExternalLenses/SLACS/slacs_all_params.csv'
audit=json.loads((HERE.parent/'slacs-resolved-input-audit/results.json').read_text(encoding='utf-8'))
expected=next(r['sha256'] for r in audit['additional_sources'] if r['path'].endswith('slacs_all_params.csv'));assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
fits=json.loads((HERE.parent/'slacs-resolved-fit/results.json').read_text(encoding='utf-8'));names={'SDSS'+r['Name'] for r in fits['rows']}
old={r['Name']:r for r in json.loads((HERE.parent/'lensing-data-readiness/lens-observations-and-image-models.json').read_text(encoding='utf-8'))}
rows=[]
for row in csv.DictReader(source.open(encoding='utf-8')):
 if row['name'] not in names:continue
 name=row['name'].removeprefix('SDSS');kind=row['light_model']
 if kind not in ['single_sersic','double_sersic']:
  rows.append(dict(Name=name,status='missing_published_components',previous_Re_arcsec=old[name]['Re'],fallback_radius_not_adopted=row['r_eff_auger']));continue
 components=[]
 for j in range(1,2 if kind=='single_sersic' else 3):
  R=float(row[f'r_sersic_{j}']);n=float(row[f'n_sersic_{j}']);amp=float(row[f'amp_sersic_{j}']);e1=float(row[f'e1_sersic{j}']);e2=float(row[f'e2_sersic{j}']);e=np.hypot(e1,e2);q=(1-e)/(1+e)
  assert R>0 and n>0 and amp>0 and 0<q<=1
  bn=1.9992*n-.3271
  luminosity=2*np.pi*amp*R*R*n*np.exp(bn)*gamma(2*n)/bn**(2*n)
  components.append(dict(R_arcsec=R,n=n,amp_at_R=amp,e1=e1,e2=e2,q=float(q),bn=bn,relative_total_light=float(luminosity)))
 total=sum(r['relative_total_light'] for r in components)
 def enclosed(R):return sum(r['relative_total_light']*gammainc(2*r['n'],r['bn']*(R/r['R_arcsec'])**(1/r['n'])) for r in components)/total
 half=brentq(lambda R:enclosed(R)-.5,1e-6,1e4)
 for r in components:r['total_light_fraction']=r['relative_total_light']/total
 rows.append(dict(Name=name,status='conditional_components_available',components=components,computed_equal_area_half_light_arcsec=half,tabulated_r_eff_arcsec=float(row['r_eff']),computed_over_tabulated=half/float(row['r_eff']),previous_Re_arcsec=old[name]['Re'],computed_over_previous=half/old[name]['Re']))
(HERE/'results.json').write_text(json.dumps(dict(scope='Image-profile interpretation audit only, no mass or kinematic refit',source_sha256=expected,conventions='Conditional product-average/equal-area radius, amplitude at R; lenstronomy documented bn approximation; historical image-fitting version unverified',rows=rows),indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([{k:v for k,v in r.items() if k!='components'} for r in rows]))

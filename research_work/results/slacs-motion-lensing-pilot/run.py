"""Training-only motion-to-lensing transfer under explicit approximations."""
from pathlib import Path
import json,hashlib
import numpy as np
from model import *
HERE=Path(__file__).resolve().parent;DATA=HERE.parent/'lensing-data-readiness'
files=[DATA/'lens-observations-and-image-models.json',DATA/'conditional-geometry.json',HERE.parent/'joint-galaxy-audit/results.json',HERE/'protocol.json',HERE/'model.py',HERE/'run.py']
protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
obs=json.loads(files[0].read_text(encoding='utf-8'));geo={r['Name']:r for r in json.loads(files[1].read_text(encoding='utf-8'))}
par=json.loads(files[2].read_text(encoding='utf-8'))['sparc']['parameters']
A,p,astar=par['A'],par['p'],par['a_star_m_s2']*3.085677581491367e19/1e6
selected=[r for r in obs if r['role']=='training' and r['Mph']=='E' and r['sigma'] is not None and r['sigma']>0]
rows=[]
for i,r in enumerate(selected):
 g=geo[r['Name']];Dl=1000*g['conditional_Dl_Mpc'];ratio=g['conditional_Dls_over_Ds'];a=Dl*r['Re']/ARCSEC/1.8153
 rap=Dl*protocol['aperture_radius_arcsec']/ARCSEC;psf=Dl*protocol['gaussian_seeing_fwhm_arcsec']/ARCSEC/np.sqrt(8*np.log(2))
 for cut in protocol['outer_cutoffs_in_a']:
  cb,cc=aperture_coefficients(a,rap,psf,A,p,astar,cut)
  for model,coefficient,amp in [('baryons',0.,0.),('empirical_extra',cc,A)]:
   masses=[];predictions=[]
   for sigma in [r['sigma']-r['e_sigma'],r['sigma'],r['sigma']+r['e_sigma']]:
    mass=mass_from_sigma(sigma,cb,coefficient,p);masses.append(mass);predictions.append(angle(mass,a,Dl,ratio,amp,p,astar,cut))
   rows.append(dict(Name=r['Name'],role=r['role'],model=model,cutoff_in_a=cut,light_axis_ratio=r['b/a'],catalog_good_dispersion=r['Good?'],measured_sigma_kms=r['sigma'],sigma_error_kms=r['e_sigma'],scale_a_kpc=a,inferred_stellar_mass_Msun=masses[1],mass_at_sigma_minus_plus=[masses[0],masses[2]],predicted_einstein_arcsec=predictions[1],angles_at_sigma_minus_plus=[predictions[0],predictions[2]],catalog_SIE_arcsec=r['bSIE'],prediction_over_catalog=None if predictions[1] is None else predictions[1]/r['bSIE']))
 print(json.dumps(dict(completed=i+1,total=len(selected),Name=r['Name'])),flush=True)
summary=[]
for model in ['baryons','empirical_extra']:
 for cut in protocol['outer_cutoffs_in_a']:
  group=[r for r in rows if r['model']==model and r['cutoff_in_a']==cut];values=np.array([r['prediction_over_catalog'] for r in group if r['prediction_over_catalog'] is not None])
  summary.append(dict(model=model,cutoff_in_a=cut,n=len(group),missing_roots=sum(r['prediction_over_catalog'] is None for r in group),median_prediction_over_catalog=float(np.median(values)),median_absolute_fractional_discrepancy=float(np.median(abs(values-1))),rms_fractional_discrepancy=float(np.sqrt(np.mean((values-1)**2)))))
(HERE/'results.json').write_text(json.dumps(dict(scope=protocol['scope'],selection=protocol['selection'],parameters=par,summary=summary,rows=rows,hashes={str(f.relative_to(HERE.parents[2])):hashlib.sha256(f.read_bytes()).hexdigest() for f in files}),indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summary),flush=True)

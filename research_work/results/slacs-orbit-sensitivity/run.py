from pathlib import Path
import json,hashlib,importlib.util
import numpy as np
from model import coefficients,G
HERE=Path(__file__).resolve().parent;PILOT=HERE.parent/'slacs-motion-lensing-pilot'
spec=importlib.util.spec_from_file_location('pilot',PILOT/'model.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
protocol=json.loads((HERE/'protocol.json').read_text(encoding='utf-8'));betas=protocol['beta_values']
base=json.loads((PILOT/'results.json').read_text(encoding='utf-8'));par=base['parameters'];A,p,astar=par['A'],par['p'],par['a_star_m_s2']*3.085677581491367e19/1e6
geo={r['Name']:r for r in json.loads((HERE.parent/'lensing-data-readiness/conditional-geometry.json').read_text(encoding='utf-8'))}
selected=[r for r in base['rows'] if r['model']=='empirical_extra' and r['cutoff_in_a']==20]
rows=[];zero_errors=[]
for i,r in enumerate(selected):
 Dl=1000*geo[r['Name']]['conditional_Dl_Mpc'];ratio=geo[r['Name']]['conditional_Dls_over_Ds'];a=r['scale_a_kpc'];rap=Dl*1.5/old.ARCSEC;psf=Dl*1.4/old.ARCSEC/np.sqrt(8*np.log(2))
 coefs=coefficients(a,rap,psf,A,p,astar,20,betas)
 for beta,(cb,cc) in zip(betas,coefs):
  for label,c,amp in [('baryons',0.,0.),('empirical_extra',cc,A)]:
   mass=old.mass_from_sigma(r['measured_sigma_kms'],cb,c,p);theta=old.angle(mass,a,Dl,ratio,amp,p,astar,20)
   rows.append(dict(Name=r['Name'],beta=beta,model=label,mass_Msun=mass,predicted_angle_arcsec=theta,catalog_angle_arcsec=r['catalog_SIE_arcsec'],prediction_over_catalog=None if theta is None else theta/r['catalog_SIE_arcsec']))
   if beta==0 and label=='empirical_extra':zero_errors.append(abs(theta/r['predicted_einstein_arcsec']-1))
 if (i+1)%10==0:print(json.dumps(dict(completed=i+1,total=len(selected))),flush=True)
assert max(zero_errors)<1e-10
summaries=[]
for label in ['baryons','empirical_extra']:
 for beta in betas:
  group=[r for r in rows if r['model']==label and r['beta']==beta];vals=np.array([r['prediction_over_catalog'] for r in group if r['prediction_over_catalog'] is not None])
  summaries.append(dict(model=label,beta=beta,n=len(group),missing_roots=len(group)-len(vals),median_ratio=float(np.median(vals)),median_absolute_fractional_discrepancy=float(np.median(abs(vals-1))),rms_fractional_discrepancy=float(np.sqrt(np.mean((vals-1)**2)))))
checks=[]
for beta,coef in zip(betas,coefficients(3,np.inf,0,A,p,astar,20,betas)):
 error=abs(coef[0]/(G*1e11/(18*3))-1);assert error<1.5e-4
 checks.append(dict(beta=beta,global_virial_relative_error=error))
out=dict(scope=protocol['scope'],summaries=summaries,rows=rows,maximum_isotropic_reproduction_error=max(zero_errors),virial_checks=checks,hashes={str(f.relative_to(HERE.parents[2])):hashlib.sha256(f.read_bytes()).hexdigest() for f in [HERE/'model.py',HERE/'run.py',HERE/'protocol.json',PILOT/'model.py',PILOT/'results.json']})
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(summaries),flush=True)

from pathlib import Path
import os
import json,csv,re,hashlib,numpy as np
root=Path(__file__).resolve().parents[3];out=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',root/'research_work/generated'))/'energy';data=root/'companion_causal_test/data/themis.dat';rows=list(csv.DictReader((root/'companion_causal_test/energy.csv').open()))
out.mkdir(parents=True,exist_ok=True)
canon=lambda s:re.sub(r'(?<=\D)0+(?=\d)','',re.sub(r'[^A-Z0-9]','',s.upper()))
dust={}
for l in data.read_text().splitlines():
 try:
  name=l[:23].strip();ra=float(l[24:33]);dec=float(l[34:43]);dist=float(l[49:61]);lum=float(l[100:109]);a,b=np.deg2rad([ra,dec])
  if lum>0 and dist>0:dust[canon(name)]=(dist*np.array([np.cos(b)*np.cos(a),np.cos(b)*np.sin(a),np.sin(b)]),lum)
 except ValueError:pass
KAPPA=7.7315e-5*3.261563777167433;KPC=3.085677581491367e19;LSUN=3.828e26;C=299792458.;calc=[]
for row in rows:
 target=canon(row['galaxy']);pos=dust[target][0];fluxes=np.zeros(3)
 for name,(xyz,lum) in dust.items():
  if name==target:continue
  sep=max(float(np.linalg.norm(xyz-pos)),.03);x=KAPPA*sep;norm=lum*LSUN/(4*np.pi*(sep*1000*KPC)**2)
  fluxes+=norm*np.array([x*np.exp(-2*x),(-np.expm1(-x))*np.exp(-x),(-np.expm1(-x))])
 reproduced=fluxes[0]/C;original=float(row['external_companion_u_J_m3']);assert np.isclose(reproduced,original,rtol=1e-12)
 ratios=fluxes/fluxes[0];local=float(row['local_conversion_shortfall']);ext=float(row['catalog_external_shortfall']);full=float(row['full_conversion_shortfall'])
 calc.append({'galaxy':row['galaxy'],'shared_flux_W_m2':float(fluxes[0]),'no_loss_arrival_S_flux_gain':float(ratios[1]),'no_loss_arrival_1_flux_gain':float(ratios[2]),'shared_combined_shortfall':1/(1/local+1/ext),'no_loss_arrival_S_combined_shortfall':1/(1/local+ratios[1]/ext),'no_loss_arrival_1_combined_shortfall':1/(1/local+ratios[2]/ext),'full_conversion_shortfall':full,'constant_present_luminosity_100pct_supply_time_yr':full*1e10})
summary={'status':'conditional current-catalog sensitivity, not a derived no-loss propagation model','sources':len(dust),'matched_targets':len(rows),'source_sha256':hashlib.sha256(data.read_bytes()).hexdigest(),'assumptions':['Same 10 Gyr, catalog geometry and luminosities as original; no new cosmic history.','Companion per-emitted-photon energy fraction 1-exp(-x), assuming complete transfer of photon loss and no companion loss.','Compare two stipulated companion arrival factors, S=exp(x) and 1. Neither is yet derived from the candidate interaction.','Local contribution kept at archived already-optimistic first-transfer energy.','No additional sources, capture gains, or halo fits.'],'medians':{k:float(np.median([r[k] for r in calc])) for k in calc[0] if k!='galaxy'},'maximum_flux_gains':{k:max(r[k] for r in calc) for k in ['no_loss_arrival_S_flux_gain','no_loss_arrival_1_flux_gain']}}
(out/'nearby-no-loss-comparison.json').write_text(json.dumps(summary,indent=2));
with (out/'nearby-no-loss-comparison.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(calc[0]));w.writeheader();w.writerows(calc)
print(json.dumps(summary,indent=2))

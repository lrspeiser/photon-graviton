from pathlib import Path
import numpy as np,json,csv
P=Path(__file__).resolve().parent
# Hand transcription: Ettori 1805.00035v3 table1 radii, Eckert 1805.00034 table2 hydrostatic masses/gas fractions.
# Exclude A2029 because mass values differ across the tables; HydraA has no entry in gas table2.
a=[('A85',1.235,5.65,.18,.150,.005),('A644',1.230,5.66,.48,.132,.012),('A1644',1.054,3.48,.20,.128,.008),('A1795',1.153,4.63,.14,.139,.005),('A2142',1.424,8.95,.26,.158,.005),('A2255',1.196,5.26,.34,.153,.011),('A2319',1.346,7.31,.28,.189,.008),('A3158',1.123,4.26,.18,.145,.007),('A3266',1.430,8.80,.57,.132,.009),('RXC1825',1.105,4.08,.13,.133,.005),('ZW1215',1.358,7.66,.52,.106,.008)]
s=json.loads((P/'results.json').read_text());A=10**s['cross_galaxy_outer_mass_scalings']['area']['log10_normalization'];out=[]
for n,R,M,eM,fg,efg in a:
 row=dict(cluster=n,R500_Mpc=R,M500_1e14Msun=M,eM500_1e14Msun=eM,fgas500=fg,efgas500=efg,galaxy_area_prediction_Msun=A*(1000*R)**2)
 for fs in [.01,.02,.03]:
  mx=M*1e14*(1-fg-fs);row[f'extra_mass_fstar{fs}_Msun']=mx;row[f'required_area_gain_fstar{fs}']=mx/row['galaxy_area_prediction_Msun']
 row['energy_J_eta1']=row['extra_mass_fstar0.02_Msun']*1.98847e30*299792458.**2
 out.append(row)
with (P/'cluster_comparison.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
result={'n':len(out),'fstar_assumed':[.01,.02,.03],'median_required_area_gain':{str(fs):float(np.median([r[f'required_area_gain_fstar{fs}'] for r in out])) for fs in [.01,.02,.03]},'range_required_area_gain_fstar002':[min(r['required_area_gain_fstar0.02'] for r in out),max(r['required_area_gain_fstar0.02'] for r in out)],'median_extra_mass_Msun':float(np.median([r['extra_mass_fstar0.02_Msun'] for r in out])),'median_energy_J_eta1':float(np.median([r['energy_J_eta1'] for r in out])),'logmass_RMSE_dex':float(np.sqrt(np.mean([np.log10(r['required_area_gain_fstar0.02'])**2 for r in out]))),'excluded':{'HydraA':'No gas fraction row in Eckert Table2','A2029':'Different mass values in the two source tables; not silently combined'},'status':'External sample scaling check, not a joint likelihood. Published HSE/standard-distance assumptions retained. fstar is a sensitivity parameter, not new measured stellar masses. Galaxy rmax is an observing limit, not the same definition as cluster R500.'}
(P/'cluster_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))

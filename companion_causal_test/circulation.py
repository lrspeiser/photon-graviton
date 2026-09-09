from pathlib import Path
import json,csv,numpy as np
P=Path(__file__).resolve().parent;gamma=7.7315e-11;k=7.7315e-5;p=.26390611655880836;c=299792458.
e=list(csv.DictReader((P/'energy.csv').open()));u_required=np.array([float(x['external_companion_u_J_m3'])*float(x['catalog_external_shortfall']) for x in e]);u_ebl=4*np.pi*50e-9/c
rows=[]
for C in [1000,10000,30000,100000]:
 x=k*p*C;loops=(1/(k*p))/C
 if x<1:
  n=(1-x)**(-1/p);ph=1/n;co=np.log(n)/n;ti=1-(1+np.log(n))/n
 else:n=ph=co=ti=None
 rows.append(dict(circumference_Mly=C,future_circuit_budget=loops,first_circuit_n=n,photon_fraction=ph,free_companion_fraction=co,timing_field_fraction=ti))
results={'original_exponential_at_constant_path_speed_gamma_lifetime_yr':1/gamma,'revised_history_future_path_Mly':1/(k*p),'circulation_examples':rows,'EBL_benchmark_I_nW_m2_sr':50,'EBL_energy_density_J_m3':u_ebl,'required_external_energy_density_median_J_m3':float(np.median(u_required)),'required_over_EBL_median':float(np.median(u_required)/u_ebl),'capture_duration_yr_if_companion_u_equals_EBL':float(1e10*np.median(u_required)/u_ebl),'nonredshifting_companion_build_time_yr_steady_original_gamma_ignoring_capture':float(np.median(u_required)/u_ebl/gamma),'notes':['EBL source Driver2016 arXiv1605.01523 integrated COB24+CIB26 nW/m2/sr; benchmark not companion measurement.','Future path derives from unbounded extrapolation dot n=gamma n^p, v=c/n in static geometry; not an observational horizon measurement.','For shared secondary redshift, no-capture companion fraction is log(n)/n, not 1-1/n.','Constant-gamma uniform steady bath obeys uc=gamma*ugamma/(gamma+lambda); circulating topology alone does not evade this local balance.']}
(P/'circulation_results.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))

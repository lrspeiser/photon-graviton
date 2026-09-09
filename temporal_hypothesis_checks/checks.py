import json, math, hashlib
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import lambertw
P=Path(__file__).resolve().parent/'source_inputs'
O=Path(__file__).resolve().parent
rows=json.loads((P/'features.json').read_text())
labels={}
for split in ['train','validation','test']: labels.update(json.loads((P/(split+'_labels.json')).read_text()))
c=299792.458
# Each new variant fits one slope on the original 104 training groups only.
# Original held-out labels have already been seen; these are exploratory evaluations.
def pred(model,k,rr):
 d=np.array([r['distance_mpc'] for r in rr]); sd=np.array([r['sd_mpc'] for r in rr]); x=k*d/c
 if model=='positive_exponential_history':
  z=2*x/(1+np.sqrt(1+4*x)); der=k/np.sqrt(1+4*x)
 else:
  a=.5 if model=='energy_only' else 1.
  w=lambertw(a*x).real; z=np.expm1(w/a); der=k*np.exp((1-a)*w/a)/(1+w)
 return c*z,np.sqrt(300**2+(der*sd)**2)
train=[r for r in rows if r['split']=='train']
y=np.array([labels[str(r['pgc'])] for r in train])
out={'status':'Exploratory reuse of already examined splits; external checks use published summary measurements, not raw-data reanalyses.','models':{}}
for model in ['full_signal','energy_only','positive_exponential_history']:
 def obj(k):
  v,s=pred(model,k,train)
  return np.sum(np.log(s)+.5*((v-y)/s)**2)
 opt=minimize_scalar(obj,bounds=(10,150),method='bounded'); k=float(opt.x)
 result={'k_kms_mpc':k,'scores':{}}
 for split in ['train','validation','test']:
  rr=[r for r in rows if r['split']==split]; v,s=pred(model,k,rr); yy=np.array([labels[str(r['pgc'])] for r in rr]); e=v-yy
  result['scores'][split]={'n':len(rr),'rmse_kms':float(np.sqrt(np.mean(e**2))),'mae_kms':float(np.mean(abs(e)))}
 out['models'][model]=result
b=1.003; sig=math.hypot(.005,.010)
gamma=7.731496595524618e-11
T0=2.72548; arad=7.5657e-16; u=arad*T0**4
out['external_checks']={'DES_duration':{'b':b,'sigma_quadrature':sig,'b1_standardized_difference':(b-1)/sig,'stretch_z1':2**b,'stretch_z1_sigma':2**b*math.log(2)*sig,'note':'Published method-two estimate is partly a consistency check; do not interpret Gaussian extrapolation to b=0 as rigorous exclusion significance.'},'conditional_alpha_coupling':{'gamma_per_year':gamma,'clock_alpha_mean_per_year':1e-18,'clock_alpha_sigma_per_year':1.1e-18,'conservative_approx_95_abs_q_bound':(1e-18+1.96*1.1e-18)/gamma,'note':'Defines alpha_dot/alpha=q*gamma, only for a model with this coupling. Not a bound on common clock rescaling; epsilon=n epsilon0 and c=c0/n can cancel in alpha.'},'cmb_energy_budget':{'T0_K':T0,'u_J_m3':u,'Q_W_m3_if_all_CMB_couples':u*gamma/(365.25*86400),'energy_transferred_in_1Gyr_fraction_initial':1-math.exp(-gamma*1e9),'note':'Constant present fractional transfer rate and static volume; no replenishment. Not a self-consistent cosmic field evolution.'},'rate_zero_point_degeneracy_mag':5*math.log10(75.59797232699/69.51915447584)}
out['input_sha256']={n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['features.json','train_labels.json','validation_labels.json','test_labels.json']}
(O/'results.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))

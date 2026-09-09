from pathlib import Path
import json,csv,hashlib,statistics
from decimal import Decimal, getcontext
getcontext().prec=90
class mp:
 mpf=Decimal
 exp=staticmethod(lambda x:x.exp())
 log=staticmethod(lambda x:x.ln())
 expm1=staticmethod(lambda x:x.exp()-1)
 log1p=staticmethod(lambda x:(1+x).ln())
 diff=staticmethod(lambda f,x:(f(x+Decimal(1))-f(x-Decimal(1)))/2)
P=Path(__file__).resolve().parent
c=mp.mpf('299792458');au=mp.mpf('149597870700');year=mp.mpf('31557600');g=mp.mpf('7.731496595524618e-11')/year
# Static material-coordinate endpoints, negligible delay, fixed coherent ratio.
# clock frequency per coordinate second f_atom=f0*n**(-p).
def observable(R,p,t2=mp.mpf(0),history='rolling'):
 if history=='rolling':
  n2=mp.exp(g*t2);S=1+2*g*R*n2/c
  tau=mp.log(S)/g if p==0 else n2**(-p)*mp.expm1(p*mp.log(S))/(p*g)
 else:
  n2=1+g*t2;S=mp.exp(2*g*R/c)
  tau=mp.log(S)/g if p==1 else n2**(1-p)*(-mp.expm1(-(1-p)*mp.log(S)))/(g*(1-p))
 y=mp.expm1((p-1)*mp.log(S))
 return S,y,c*tau/2
rows=[]
for h in ['rolling','affine']:
 for p in [0,1,2]:
  for a in [1,10,30,150]:
   R=a*au;S,y,r=observable(R,p,history=h)
   _,_,rplus=observable(R,p,year/2,h);_,_,rminus=observable(R,p,-year/2,h)
   rows.append({'history':h,'clock_exponent_p':p,'distance_AU':a,'two_way_fractional_frequency_residual':float(y),'inferred_range_rate_mm_s':float(-c*y/2*1000),'fixed_geometry_range_change_per_year_m':float(rplus-rminus),'range_offset_at_reference_epoch_m':float(r-R)})
# Independent propagation + phase check for constant clocks using arbitrary precision.
R=10*au;te=-mp.log1p(2*g*R/c)/g;te2=te+1
tr2=-mp.log(mp.exp(-g*te2)-2*g*R/c)/g
S,y,r=observable(R,0)
checks={'one_second_pulse_stretch':str(tr2),'differential_frequency_stretch':str(S),'finite_pulse_relative_difference':float(tr2/S-1),'p1_cancellation_exact_to_precision':all(abs(observable(a*au,1,history=h)[1])<mp.mpf('1e-50') and abs(observable(a*au,1,history=h)[2]-a*au)<mp.mpf('1e-40') for h in ['rolling','affine'] for a in [1,10,30]),'range_rate_identity_residual_m_per_s':float(mp.diff(lambda t:observable(R,0,t)[2],mp.mpf(0))+c*y/2)}
(P/'static_link_results.json').write_text(json.dumps({'gamma_per_year':float(g*year),'rows':rows,'checks':checks},indent=2))
# Real reconstructed trajectory geometry only. No tracking observations or residual fit.
s=json.loads((P/'data/horizons_response.json').read_text())['result'];geometry=[]
for line in s.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines():
 f=[x.strip() for x in line.split(',')];R=mp.mpf(f[9])*1000;S,y,r=observable(R,0)
 geometry.append({'JD_TDB':float(f[0]),'date':f[1],'geometric_range_AU':float(R/au),'geometric_range_rate_kms':float(f[10]),'stationary_endpoint_two_way_y_forecast':float(y),'stationary_endpoint_equivalent_mm_s':float(-c*y/2*1000),'stationary_endpoint_range_drift_m_per_year':float(g*year*R)})
with (P/'cassini_geometry_forecast.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(geometry[0]));w.writeheader();w.writerows(geometry)
summary={'N_epochs':len(geometry),'range_AU':[min(r['geometric_range_AU'] for r in geometry),max(r['geometric_range_AU'] for r in geometry)],'two_way_y':[min(r['stationary_endpoint_two_way_y_forecast'] for r in geometry),max(r['stationary_endpoint_two_way_y_forecast'] for r in geometry)],'equivalent_mm_s':[min(r['stationary_endpoint_equivalent_mm_s'] for r in geometry),max(r['stationary_endpoint_equivalent_mm_s'] for r in geometry)],'mean_removed_y_rms':statistics.pstdev(r['stationary_endpoint_two_way_y_forecast'] for r in geometry),'status':'Forecast evaluated at genuine JPL reconstructed ranges; not real Doppler residuals, not a moving-endpoint light-time model, and not an observational fit.'}
(P/'geometry_summary.json').write_text(json.dumps(summary,indent=2))
(P/'checksums.json').write_text(json.dumps({str(f.relative_to(P)):hashlib.sha256(f.read_bytes()).hexdigest() for f in P.rglob('*') if f.is_file() and f.name!='checksums.json'},indent=2))
print(json.dumps({'static_p0':[r for r in rows if r['history']=='rolling' and r['clock_exponent_p']==0],'checks':checks,'geometry':summary},indent=2))

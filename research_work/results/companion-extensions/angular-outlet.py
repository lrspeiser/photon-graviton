"""Angular empty-channel fraction versus repeated deflection diagnostic."""
from pathlib import Path
import json,math
from scipy.integrate import quad
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
old=json.loads((P/'quantum-bridge-results.json').read_text())
N=2*old['mean_energy_loss_fraction']/1e-8
RSUN=6.957e8;PC=3.085677581491367e16;ARC=180*3600/math.pi
# Angular translation d in units of the illuminated disk radius.
def outside(d):
 if d>=2:return 1.
 return 1-(2*math.acos(d/2)-d*math.sqrt(1-d*d/4))/math.pi
out=dict(scope='Uniform single-source angular disk, fixed-size isotropic angular kicks, equal inside-mode occupations; local geometry and illustrative repeated-kick diagnostic',minimum_net_steps=N,rows=[],checks=[])
for nh in [1e-12,1.,1000.]:
 for target in [.1,.9,.99]:
  f=(1-target)*(1+nh)/(1+(1-target)*nh)
  d=brentq(lambda d:outside(d)-f,0,2)
  # Independent intersection area by integrating intersecting vertical chords.
  area=quad(lambda x:2*min(math.sqrt(max(0,1-x*x)),math.sqrt(max(0,1-(x-d)**2))),d-1,1,points=[d/2],epsabs=1e-10)[0]
  assert abs(1-area/math.pi-f)<1e-8
  r=(1-f)*(1+nh)/(1+(1-f)*nh)
  assert abs(r-target)<1e-12
  for distance in [1/206264.80624709636,1.,1e3,1e6]:
   radius=math.asin(RSUN/(distance*PC));theta=d*radius
   # Fixed independent isotropic kicks: <cos total angle> = (cos theta)^N.
   logcorr=N*math.log1p(-2*math.sin(theta/2)**2)
   proxy=math.sqrt(N)*theta
   out['rows'].append(dict(input_occupation=nh,target_reverse_forward=target,empty_fraction=f,step_over_source_radius=d,source_distance_pc=distance,step_arcsec=theta*ARC,small_angle_direction_rms_arcsec=proxy*ARC if proxy<.1 else None,direction_correlation=math.exp(logcorr),one_minus_direction_correlation=-math.expm1(logcorr)))
  # Source distance needed for a one-arcsec RMS proxy, under fixed local geometry.
  min_distance=RSUN/(math.sin((1/ARC)/(d*math.sqrt(N))))/PC
  out['checks'].append(dict(nh=nh,target=target,quadrature_error=abs(1-area/math.pi-f),minimum_fixed_geometry_distance_pc_for_one_arcsec=min_distance))
(P/'angular-outlet-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print('steps',N)
for r in out['rows']:
 if r['input_occupation']==1e-12 and r['target_reverse_forward']==.9:print(r)
print(out['checks'][1])

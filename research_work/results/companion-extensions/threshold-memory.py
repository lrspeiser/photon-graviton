"""Rate-independent filling versus post-source retention tradeoff."""
from pathlib import Path
import math,json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
P=Path(__file__).resolve().parent
history=json.loads((P/'threshold-history-results.json').read_text())
base=[r for r in history['rows'] if r['duration'] is None]
L=6*math.log(10)
def w(y):return math.sin(math.pi/3)/(2*math.pi*(math.cosh(y/3)+.5))
norm=quad(w,-L,L)[0]
out=dict(scope='Stipulated common-rate threshold kinetics, source removed after 90 percent equilibrium formation; no physical ages or observed fading histories',rows=[])
for row in base:
 X=row['X'];u=row['duration_to_90_percent']
 def retained(v):
  return quad(lambda y: X/(X+math.exp(y))*(-math.expm1(-(X+math.exp(y))*u))*math.exp(-math.exp(y)*v)*w(y),-L,L,epsabs=1e-11)[0]/norm
 initial=retained(0);assert abs(initial/.9/row['equilibrium_fraction']-1)<1e-8
 vals={}
 for fraction in [.9,.5]:
  upper=1.
  while retained(upper)>fraction*initial:upper*=10
  v=brentq(lambda v:retained(v)-fraction*initial,0,upper,xtol=1e-10)
  assert abs(retained(v)/initial-fraction)<1e-8
  vals[str(fraction)]=v
 out['rows'].append(dict(galaxy=row['galaxy'],split=row['split'],X=X,formation_duration=u,initial_occupied_fraction=initial,post_source_duration_to_90_percent=vals['0.9'],post_source_half_life=vals['0.5'],max_dark_to_build_ratio_for_90_percent=vals['0.9']/u,half_life_to_build_ratio=vals['0.5']/u))
slow=max(out['rows'],key=lambda r:r['formation_duration']);fastloss=min(out['rows'],key=lambda r:r['post_source_duration_to_90_percent'])
out['common_history_bound']=dict(slowest_formation_galaxy=slow['galaxy'],largest_formation_duration=slow['formation_duration'],fastest_loss_galaxy=fastloss['galaxy'],smallest_loss_duration=fastloss['post_source_duration_to_90_percent'],maximum_common_dark_to_build_ratio=fastloss['post_source_duration_to_90_percent']/slow['formation_duration'])
out['individual_ratio_range']=[min(r['max_dark_to_build_ratio_for_90_percent'] for r in out['rows']),max(r['max_dark_to_build_ratio_for_90_percent'] for r in out['rows'])]
(P/'threshold-memory-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(out['common_history_bound']);print(out['individual_ratio_range']);print(min(out['rows'],key=lambda r:r['X']));print(max(out['rows'],key=lambda r:r['X']))

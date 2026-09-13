"""Analytic statistics of a proposed Poisson inelastic loss process."""
from pathlib import Path
import json,math
P=Path(__file__).resolve().parent;alpha=.0002488993286382367
rows=[]
for D in [30.660139,100.,1000.]:
 A=alpha*D
 for width in [1e-3,1e-5,1e-6]:
  delta=math.log1p(width*width)/A;events=A/delta
  rows.append(dict(distance_Mpc=D,mean_energy_loss_fraction=-math.expm1(-A),mean_energy_redshift=math.expm1(A),illustrative_relative_line_width=width,maximum_fraction_lost_per_event=delta,minimum_mean_events=events,mean_log_energy_loss=-events*math.log1p(-delta),variance_log_energy=events*math.log1p(-delta)**2))
  assert abs(math.sqrt(math.expm1(A*delta))/width-1)<1e-12
(P/'conversion-results.json').write_text(json.dumps(dict(alpha_per_Mpc=alpha,scope='Illustrative linewidth tolerances, not measured exclusions; ensemble mean energy redshift is not identical to mean individual redshift',rows=rows),indent=2),encoding="utf-8",newline="\n")
print(json.dumps(rows[:3],indent=2))

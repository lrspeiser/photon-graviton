"""Rate cost of matching source-size angular moments with Gaussian overlap."""
from pathlib import Path
import math,json,importlib.util
from scipy.constants import hbar,c,electron_volt
from scipy.special import erfcx
from scipy.integrate import quad
P=Path(__file__).resolve().parent;PC=3.085677581491367e16
old=json.loads((P/'angular-outlet-results.json').read_text())
d=next(r['step_over_source_radius'] for r in old['rows'] if r['input_occupation']==1e-12 and r['target_reverse_forward']==.9)
R=6.957e8/PC;E=2.;D=30660139.;alpha=.0002488993286382367/1e6
spec=importlib.util.spec_from_file_location('spatial',P/'spatial-response.py');spatial=importlib.util.module_from_spec(spec);spec.loader.exec_module(spatial)
out=dict(scope='Match narrow-angle RMS of prescribed source-tracking kernel using earlier Gaussian overlap; fixed optical energy and unchanged coupling/target density',rows=[],path_checks=[])
for delta in [1e-8,2e-10]:
 spatial.gap=delta
 sc=d*R*math.sqrt(E*(E-delta))/delta
 def logS(s):
  theta=d*R/s;a=1/(theta*math.sqrt(E*(E-delta)))
  return -(a*delta)**2+math.log(3/(8*a*a*E*(E-delta)))
 for s in [.1,1.,10.,1000.,1e6]:
  theta=d*R/s;a=1/(theta*math.sqrt(E*(E-delta)))
  LS=logS(s)
  if s<=1:
   direct,u=spatial.angular(E,a)
   assert abs(math.log(direct)-LS)<1e-10
   assert abs(2*u/theta**2-1)<1e-10
  out['rows'].append(dict(transfer_eV=delta,source_distance_pc=s,coherence_sigma_m=a*hbar*c/electron_volt,longitudinal_mismatch=a*delta,log10_rate_suppression=LS/math.log(10),log10_rate_relative_to_one_pc=(LS-logS(1))/math.log(10),cutoff_distance_pc=sc))
 # Rate ratio normalized at 1 pc, integrated to the old path length.
 effective=quad(lambda s:math.exp(-((s/sc)**2-(1/sc)**2))/s**2,1,max(1+20*sc,20*sc),epsabs=1e-12)[0]
 analytic=1-math.sqrt(math.pi)/sc*erfcx(1/sc)
 assert abs(effective/analytic-1)<1e-9
 z=math.expm1(alpha*effective)
 out['path_checks'].append(dict(transfer_eV=delta,cutoff_distance_pc=sc,effective_path_pc=effective,infinite_path_analytic_pc=analytic,redshift_if_alpha_at_one_pc_is_reference=z,uniform_reference_redshift=math.expm1(alpha*(D-1)),normalization_increase_to_match_reference_path=(D-1)/effective))
(P/'kernel-rate-closure-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
print(out['path_checks'])
for r in out['rows']:
 if r['source_distance_pc'] in [1,1000]:print(r)

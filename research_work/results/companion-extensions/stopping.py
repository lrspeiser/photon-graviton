"""Energy landscape along frozen circular endpoints; no dynamical claim."""
import importlib.util
from pathlib import Path
spec=importlib.util.spec_from_file_location('coupled',Path(__file__).with_name('coupled-torque.py'))
ct=importlib.util.module_from_spec(spec);spec.loader.exec_module(ct)
np=ct.np;json=ct.json;P=ct.P

def energy(q):
 return q['final_energy']['total']*ct.MSUN*1e6

def derivatives(rows,h,s):
 E={round(q['s'],9):energy(q) for q in rows}
 a,b,c=[E[round(x,9)] for x in [s-h,s,s+h]]
 return (c-a)/(2*h),(c-2*b+a)/h**2

if __name__=='__main__':
 s=float(ct.selected['s']);out={'scope':'Constrained endpoint-family energy, not evolution or a derived restoring interaction','s_reference':s,'cases':[]}
 for baseline in ['I','II']:
  inp=ct.prepare(baseline,8192);fineinp=ct.prepare(baseline,16384)
  for lo,hi in [(15,30),(30,60),(60,120)]:
   rows=[]
   for x in sorted(set([.4,.5,.6,.7,.8,.9,1.,s-.02,s-.01,s,s+.01,s+.02])):
    q,r=ct.solve(inp,lo,hi,x);q['s']=x;rows.append(q)
   fine=[]
   for x in [s-.01,s,s+.01]:
    q,r=ct.solve(fineinp,lo,hi,x,mesh=32768);q['s']=x;fine.append(q)
   good=all(q['converged'] for q in rows+fine)
   case=dict(baryons=baseline,receiver_band=[lo,hi],all_converged=good,rows=rows,refined_rows=fine)
   if good:
    d1,d2=derivatives(rows,.01,s);coarse1,coarse2=derivatives(rows,.02,s);f1,f2=derivatives(fine,.01,s)
    case['slope_J']=d1;case['curvature_J']=d2
    case['half_step_slope_relative_change']=d1/coarse1-1
    case['refined_slope_relative_change']=f1/d1-1
    case['refined_curvature_J']=f2
    case['energy_increases_with_s_on_sample']=all(energy(b)>energy(a) for a,b in zip(rows,rows[1:]))
    case['inverse_support']=[]
    reference=next(q for q in rows if q['s']==s)
    for p in [1,2,3]:
     A=d1*s**(p+1)/p;Aref=f1*s**(p+1)/p
     storage=A*(s**(-p)-1)
     case['inverse_support'].append(dict(p=p,A_J=A,energy_at_reference_J=A/s**p,increase_from_uncontracted_J=storage,conditional_remaining_release_J=reference['released_energy_J']-storage,local_curvature_J=d2+(p+1)*d1/s,refined_local_curvature_J=f2+(p+1)*f1/s,positive_coefficient=bool(A>0),refined_A_J=Aref))
   out['cases'].append(case)
   print(json.dumps({k:v for k,v in case.items() if k not in ['rows','refined_rows']}),flush=True)
 # Analytic power-law finite differences independently check sign/convention.
 h=1e-5;f=lambda x:-2/x
 assert abs((f(s+h)-f(s-h))/(2*h)-2/s**2)<1e-7
 (P/'stopping-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

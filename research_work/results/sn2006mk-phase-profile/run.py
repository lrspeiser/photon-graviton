from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import minimize_scalar
H=Path(__file__).resolve().parent;R=H.parents[2]
loader=H.parent/'sn2006mk-lightcurve-pilot/run.py';s=loader.read_text(encoding='utf-8');cut='\nresults=[]\n';assert s.count(cut)==1
ns={'__file__':str(loader),'__name__':'input_loader'};exec(compile(s.split(cut)[0],str(loader),'exec'),ns)
phot,eligible,interp,epoch,stretches,start,t=[ns[k] for k in ['phot','eligible','interp','epoch','stretches','start','t']]
oldpath=H.parent/'sn2006mk-lightcurve-pilot/results.json';old=json.loads(oldpath.read_text(encoding='utf-8'))
cache=[]
for band in ['R','I']:
 p=R/'research_work/generated/essence-passbands'/f'CTIO4m_{band}.dat';d=np.loadtxt(p);wo,energy=d.T
 curve=np.array([float(np.trapezoid(interp(np.column_stack([np.full(len(wo),phase),wo/1.4754]))*energy,wo)) for phase in t])
 # Integration and linear phase interpolation commute; verify at a non-grid phase.
 probe=2.123;direct=float(np.trapezoid(interp(np.column_stack([np.full(len(wo),probe),wo/1.4754]))*energy,wo));assert np.isclose(np.interp(probe,t,curve),direct,rtol=1e-12)
 obs=sorted([x for x in eligible if x['band']==band+'4m'],key=lambda x:x['mjd']);anchor=min(obs,key=lambda x:abs(x['mjd']-epoch));others=[x for x in obs if x!=anchor]
 cache.append(dict(band=band,curve=curve,anchor=anchor,dates=np.array([x['mjd'] for x in others]),y=np.array([x['flux'] for x in others]),variance=np.array([x['error_upper']**2 for x in others])))

def evaluate(q,A):
 parts=[]
 for c in cache:
  phases=q+(c['dates']-epoch)/A;pa=q+(c['anchor']['mjd']-epoch)/A
  assert np.min(phases)>=-19 and np.max(phases)<=85 and -19<=pa<=85
  ratio=np.interp(phases,t,c['curve'])/np.interp(pa,t,c['curve']);res=c['y']-ratio*c['anchor']['flux'];v=c['variance'];u=ratio*c['anchor']['error_upper']
  factor=1+np.sum(u*u/v);chi=float(np.sum(res*res/v)-np.sum(res*u/v)**2/factor);logdet=float(np.sum(np.log(v))+np.log(factor))
  parts.append(dict(band=c['band'],chi_square=chi,log_determinant=logdet,score=chi+logdet))
 return dict(start_phase=float(q),chi_square=sum(x['chi_square'] for x in parts),score=sum(x['score'] for x in parts),parts=parts)
results=[]
for name,A in stretches.items():
 fixed=evaluate(start,A);expected=sum(x['chi_square_conditional'] for x in old['results'] if x['case']==name);assert abs(fixed['chi_square']-expected)<1e-8
 scans=[]
 for n in [301,601]:
  grid=np.linspace(-10,5,n);vals=[evaluate(q,A) for q in grid];j=int(np.argmin([x['score'] for x in vals]));candidates=[vals[0],vals[-1],vals[j]]
  if 0<j<n-1:
   opt=minimize_scalar(lambda q:evaluate(q,A)['score'],bounds=(grid[j-1],grid[j+1]),method='bounded',options={'xatol':1e-9});assert opt.success;candidates.append(evaluate(opt.x,A))
  best=min(candidates,key=lambda x:x['score']);scans.append(dict(nodes=n,best=best))
 assert abs(scans[0]['best']['score']-scans[1]['best']['score'])<1e-6
 best=scans[1]['best'];results.append(dict(case=name,event_stretch=A,fixed=fixed,profiled=best,at_boundary=abs(best['start_phase']+10)<1e-5 or abs(best['start_phase']-5)<1e-5,refinement_score_difference=abs(scans[0]['best']['score']-best['score']),grid_columns=['start_phase','chi_square','score'],grid=[[v['start_phase'],v['chi_square'],v['score']] for v in vals]))
out=dict(scope='Exploratory common source-start-phase fit to exposed photometry; not calibrated propagation inference',phase_bounds=[-10,5],fixed_phase=start,scored_points=sum(len(x['y']) for x in cache),score_definition='chi-square plus log determinant; common Gaussian constant omitted',results=results,hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [loader,oldpath,Path(__file__)]})
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps([{k:v for k,v in x.items() if k!='grid'} for x in results]))

"""Real SPARC rotation data; MOND-guided positive fixed-inventory diagnostic."""
from pathlib import Path
import json,hashlib
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
P=Path(__file__).resolve().parent;ROOT=P.parents[2];OLD=P.parent/'isotropic-galaxy-transfer'
paths=[OLD/'model-comparison-predictions.json',OLD/'model-comparison-results.json',OLD/'third-radiation-retention-results.json',ROOT/'temporal_candidate_audit/data/SPARC_Lelli2016c.mrt',P/'mond-inventory-protocol.md',Path(__file__)]
saved=json.loads(paths[0].read_text());controls=json.loads(paths[1].read_text());cp=json.loads(paths[2].read_text())['models']['attenuated']
G=4.30091727003628e-6;a0=controls['MOND_a0_SI']*3.085677581491367e19/1e6
cat={}
for line in paths[3].read_text().splitlines():
 t=line.split()
 if len(t)==19:
  try:cat[t[0]]=(float(t[11]),float(t[7]))
  except ValueError:pass
data=[]
for row in saved:
 if row['model']!='companion_third':continue
 b=next(v for v in saved if v['model']=='baryons' and v['galaxy']==row['galaxy'])
 r=np.array(row['R_kpc']);y=np.array(row['observed_kms']);vb=np.array(b['predicted_kms'])**2;ref=np.array(row['predicted_kms'])
 rd,L=cat[row['galaxy']];eta=(L/rd**2)**(1/3)/(1+(L/rd**2)**(1/3));a=cp['scale_to_disk']*rd
 gb=vb/r;gc=2*a0*gb/(np.sqrt(gb*gb+4*a0*gb)+gb) # stable g_MOND - g_b
 target=r*r*gc/G;mass=r*(ref*ref-vb)/G
 assert np.min(mass)>0 and np.all(np.diff(mass)>=-1e-6*max(mass))
 data.append(dict(name=row['galaxy'],split=row['split'],r=r,y=y,vb=vb,ref=ref,rd=rd,a=a,eta=eta,target=target,mass=mass,mond=np.sqrt(vb+r*gc)))
assert len(data)==149 and sum(len(d['r']) for d in data)==3150

def inventory(d,order):
 mu,w=leggauss(order);T=cp['k0_per_kpc']*d['a']
 def integrand(theta):
  x=np.tan(theta);b2=1+x*x*(1-mu*mu);b=np.sqrt(b2);z=x*mu
  tau=T*(z/(2*b2*(b2+z*z))+(np.arctan(z/b)+np.pi/2)/(2*b**3))
  return np.sin(theta)**2*(np.exp(-np.maximum(tau,0))@w/2)
 value,err=quad(integrand,0,np.pi/2,epsabs=1e-10,epsrel=1e-9,limit=200)
 return 4*np.pi*2*cp['C_Msun_kpc3']*d['eta']*d['a']**3*value
for d in data:
 d['total']=inventory(d,192);d['coarse_total']=inventory(d,96)
 assert d['total']>=max(d['mass'])*(1-1e-5)
 d['endpoint']=np.minimum(np.maximum.accumulate(d['target']),d['total'])
 d['coarse_endpoint']=np.minimum(np.maximum.accumulate(d['target']),d['coarse_total'])
 assert np.all(np.diff(d['endpoint'])>=0) and max(d['endpoint'])<=d['total']
def velocities(f,key='endpoint'):
 return [np.sqrt(d['vb']+G*((1-f)*d['mass']+f*d[key])/d['r']) for d in data]
def score(v):
 return {s:dict(n=sum(d['split']==s for d in data),RMSE_kms=float(np.sqrt(np.mean([np.mean((p-d['y'])**2) for d,p in zip(data,v) if d['split']==s]))),log_RMS=float(np.sqrt(np.mean([np.mean(np.log10(p/d['y'])**2) for d,p in zip(data,v) if d['split']==s])))) for s in ['train','validation','test']}
def loss(f):return score(velocities(f))['train']['log_RMS']**2
grid=np.linspace(0,1,101);scan=[loss(f) for f in grid];candidates=[(0.,scan[0]),(1.,scan[-1])]
for i in range(1,100):
 if scan[i]<=scan[i-1] and scan[i]<=scan[i+1]:
  result=minimize_scalar(loss,bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-10});assert result.success;candidates.append((float(result.x),float(result.fun)))
f,_=min(candidates,key=lambda t:t[1]);v=velocities(f);vcoarse=velocities(f,'coarse_endpoint')
maxdiff=max(float(max(abs(a-b))) for a,b in zip(v,vcoarse));assert maxdiff<.05,maxdiff
models=dict(reference=[d['ref'] for d in data],MOND_raw=[d['mond'] for d in data],positive_capped_endpoint=velocities(1),selected=v)
rows=[]
for d,p in zip(data,v):
 decline=np.diff(d['target'])< -1e-8*max(d['target'])
 rows.append(dict(galaxy=d['name'],split=d['split'],R_kpc=d['r'].tolist(),observed_kms=d['y'].tolist(),reference_kms=d['ref'].tolist(),selected_kms=p.tolist(),MOND_raw_kms=d['mond'].tolist(),total_inventory_Msun=d['total'],raw_target_Msun=d['target'].tolist(),positive_capped_target_Msun=d['endpoint'].tolist(),reference_enclosed_Msun=d['mass'].tolist(),negative_shell_intervals=int(sum(decline)),largest_target_drop_fraction=float(max(0,max(-np.diff(d['target'])))/max(d['target'])) if len(d['r'])>1 else 0.,capped_radii=int(sum(d['target']>d['total'])),endpoint_inventory_outside_last_fraction=float(1-d['endpoint'][-1]/d['total'])))
radial=[]
for name,values in models.items():
 for i,label in enumerate(['inner','middle','outer']):
  errors=[p[np.digitize(d['r']/d['rd'],[1,3])==i]-d['y'][np.digitize(d['r']/d['rd'],[1,3])==i] for d,p in zip(data,values)]
  errors=[e for e in errors if len(e)]
  radial.append(dict(model=name,bin=label,mean_equal_galaxy_error_kms=float(np.mean([e.mean() for e in errors]))))
out=dict(scope='Positive spherical equivalent redistribution of existing exact-third inventory guided by known simple MOND; exposed SPARC data; not a capture derivation or full 3D MOND solution.',f=f,a0_SI=controls['MOND_a0_SI'],scores={name:score(values) for name,values in models.items()},negative_shell_galaxies=sum(r['negative_shell_intervals']>0 for r in rows),negative_shell_intervals=sum(r['negative_shell_intervals'] for r in rows),capped_galaxies=sum(r['capped_radii']>0 for r in rows),capped_radii=sum(r['capped_radii'] for r in rows),max_inventory_refinement_fraction=max(abs(d['total']/d['coarse_total']-1) for d in data),max_prediction_refinement_kms=maxdiff,radial=radial,rows=rows,scan=[dict(f=float(a),loss=b) for a,b in zip(grid,scan)],candidates=candidates,source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
# Zero mixture must recover saved predictions; raw MOND must recover existing control.
assert max(max(abs(p-d['ref'])) for p,d in zip(velocities(0),data))<1e-10
for d in data:
 old=next(r for r in saved if r['model']=='MOND_simple_fitted' and r['galaxy']==d['name'])
 assert max(abs(d['mond']-old['predicted_kms']))<1e-9
(P/'mond-inventory-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
print('fraction',f,'negative galaxies',out['negative_shell_galaxies'],'caps',out['capped_galaxies'],'refinement',maxdiff)
print(out['scores']);print(radial)

"""Shared phase parameters and an explicit Newtonian cooling-energy ledger."""
from phase import *

def external(g,nmu):
 r=np.geomspace(3e-7,1e5,12000);mu,w=leggauss(nmu);mu=(mu+1)/2;w=w/2
 phi,fr,fz,_,_=g.field(r[:,None]*np.sqrt(1-mu**2),r[:,None]*mu)
 pb=phi@w+g.bindD(np.log(r))
 gb=(fr*np.sqrt(1-mu**2)+fz*mu)@w-G*g.massD(np.log(r))/r**2
 return r,pb,gb

def mechanical(r,dm,ext):
 idx=np.argsort(r);r=r[idx];dm=dm[idx]
 us=-G*np.sum(dm*(np.cumsum(dm)-.5*dm)/r)
 ub=np.sum(dm*np.interp(r,ext[0],ext[1]))
 kb=.5*np.sum(dm*r*np.interp(r,ext[0],ext[2]));k=kb-.5*us
 return dict(self_energy=float(us),external_energy=float(ub),support_kinetic=float(k),total=float(us+ub+k))

def ledger(data,ext,mass,s):
 edges,r,dm,sigma2,rho,vb=data
 m=mass*electron_volt/c**2;n=rho*MSUN/KPC**3/m
 rat=m*m*sigma2*1e6/(2*np.pi*hbar*hbar)*(zeta(1.5)/n)**(2/3)
 f=np.maximum(0,1-np.minimum(rat,1)**1.5)
 initial=mechanical(r,dm,ext)
 final=mechanical(np.r_[r,s*r],np.r_[dm*(1-f),dm*f],ext)
 Q=initial['total']-final['total'];scale=dm.sum()
 return dict(initial=initial,final=final,cooling_energy_Msun_kms2=float(Q),cooling_energy_J=float(Q*MSUN*1e6),cooling_per_deposit_kms2=float(Q/scale),cooling_over_rest_energy=float(Q/(scale*(c/1000)**2)),ledger_relative_residual=float((final['total']+Q-initial['total'])/max(abs(initial['total']),1)),passive=bool(Q>=-1e-10*abs(initial['total'])))

if __name__=='__main__':
 datasets={};exts={};out={'scope':'Shared two-baseline phase fit plus conditional virial cooling ledger; not a solved capture interaction or stable phase equilibrium','runs':[],'controls':[]}
 for b in ['I','II']:
  g=Galaxy(b);datasets[b]=setup(g);exts[b]=external(g,128)
  out['controls'].append(dict(baryons=b,**evaluate(datasets[b],100,1)))
 for mass in np.geomspace(1e-3,100,41):
  for s in np.geomspace(.03,1,33):
   rows=[]
   for b in ['I','II']:
    rows.append(dict(baryons=b,**evaluate(datasets[b],float(mass),float(s)),energy=ledger(datasets[b],exts[b],mass,s)))
   score=sum(q['scores']['inner']['RMSE_kms']**2 for q in rows)/2
   out['runs'].append(dict(mass_eV=float(mass),s=float(s),shared_inner_MSE=score,passive=all(q['energy']['passive'] for q in rows),rows=rows))
 best=min((q for q in out['runs'] if q['passive']),key=lambda q:q['shared_inner_MSE']);out['selected']=best
 out['refinement']=[]
 for row in best['rows']:
  b=row['baryons'];g=Galaxy(b);data=setup(g,4096,256);ext=external(g,256)
  pred=evaluate(data,best['mass_eV'],best['s']);en=ledger(data,ext,best['mass_eV'],best['s'])
  out['refinement'].append(dict(baryons=b,max_speed_change_kms=float(np.max(np.abs(np.array(pred['predicted_kms'])-row['predicted_kms']))),cooling_relative_change=float(en['cooling_energy_J']/row['energy']['cooling_energy_J']-1),energy=en))
 # Independent thin-shell self-energy and zero-change controls.
 test=mechanical(np.array([2.]),np.array([3.]),(np.array([1.,3.]),np.zeros(2),np.zeros(2)))
 assert abs(test['self_energy']-(-G*9/4))<1e-15
 for b in ['I','II']:
  q=ledger(datasets[b],exts[b],10,1)
  assert abs(q['cooling_over_rest_energy'])<1e-15
 assert max(abs(q['energy']['ledger_relative_residual']) for run in out['runs'] for q in run['rows'])<1e-12
 (P/'settling-results.json').write_text(json.dumps(out,allow_nan=False),encoding='utf-8',newline='\n')
 print(json.dumps({'mass_eV':best['mass_eV'],'s':best['s'],'rows':[{k:v for k,v in q.items() if k!='predicted_kms'} for q in best['rows']],'refinement':out['refinement'],'rejected':sum(not q['passive'] for q in out['runs'])},indent=2))

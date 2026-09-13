"""Ideal Bose fraction as a proposed selector for conservative redistribution."""
from common import *
from scipy.constants import hbar, electron_volt, c
from scipy.special import zeta
MSUN=1.98847e30;KPC=3.085677581491367e19

def setup(g,n=2048,nmu=192):
 edges=np.geomspace(g.r[1],g.r[-1],n+1);r=np.sqrt(edges[1:]*edges[:-1])
 mass=g.massD(np.log(edges));dm=np.diff(mass)
 mu,w=leggauss(nmu);mu=(mu+1)/2;w=w/2
 _,fr,fz,_,_=g.field(r[:,None]*np.sqrt(1-mu**2),r[:,None]*mu)
 sigma2=r*np.sum((fr*np.sqrt(1-mu**2)+fz*mu)*w,axis=1)/3
 rho=np.exp(g.logrho(np.log(r)))
 vb2=Robs*g.field(Robs,0)[1]-G*g.massD(np.log(Robs))/Robs
 return edges,r,dm,sigma2,rho,vb2

def evaluate(data,mev,s):
 edges,r,dm,sigma2,rho,vb2=data
 m=mev*electron_volt/c**2;nnumber=rho*MSUN/KPC**3/m
 ratio=m*m*sigma2*1e6/(2*np.pi*hbar*hbar)*(zeta(1.5)/nnumber)**(2/3)
 f=np.maximum(0,1-np.minimum(ratio,1)**1.5)
 compact=np.r_[0,np.cumsum(dm*f)];extended=np.r_[0,np.cumsum(dm*(1-f))]
 total=compact[-1]+extended[-1]
 assert abs(total/dm.sum()-1)<1e-12 and np.all((f>=0)&(f<=1))
 enclosed=np.interp(Robs/s,edges,compact)+np.interp(Robs,edges,extended)
 v=np.sqrt(vb2+G*enclosed/Robs)
 return dict(mass_eV=mev,s=s,compact_fraction=float(compact[-1]/total),scores=scores(v),predicted_kms=v.tolist(),mass_relative_error=float(total/dm.sum()-1))

if __name__=='__main__':
 out={'scope':'Ideal-gas phase selector, not equilibrium or a superfluid field solution; gas monopole matched control','runs':[],'selected':[],'controls':[],'refinement':[],'ultralight':[]}
 for variant in ['I','II']:
  g=Galaxy(variant);data=setup(g)
  control=evaluate(data,100.,1.);control['baryons']=variant;out['controls'].append(control)
  cases=[]
  for m in np.geomspace(1e-3,1e2,41):
   for s in np.geomspace(.03,1.,33):
    q=evaluate(data,float(m),float(s));q['baryons']=variant;cases.append(q)
  out['runs'].extend(cases);best=min(cases,key=lambda q:q['scores']['inner']['RMSE_kms']);out['selected'].append(best)
  refined=evaluate(setup(g,4096,256),best['mass_eV'],best['s'])
  out['refinement'].append({'baryons':variant,'max_speed_change_kms':float(np.max(np.abs(np.array(best['predicted_kms'])-refined['predicted_kms']))),'refined':refined})
  out['ultralight'].append(dict(baryons=variant,**evaluate(data,1e-22,best['s'])))
  print(variant,'control',control['scores'],'selected', {k:v for k,v in best.items() if k!='predicted_kms'},flush=True)
 out['cross_baseline_transfer']=[]
 for best in out['selected']:
  target='II' if best['baryons']=='I' else 'I'
  result=evaluate(setup(Galaxy(target)),best['mass_eV'],best['s'])
  out['cross_baseline_transfer'].append(dict(source_baseline=best['baryons'],target_baseline=target,**result))
 (P/'phase-results.json').write_text(json.dumps(out,allow_nan=False),encoding="utf-8",newline="\n")

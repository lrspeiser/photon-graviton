"""Inverse shared pressure fits followed by mass-normalized equilibria."""
from settling import *
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
REF=1e7;MSCALE=1e11
inputs=json.loads((P/'local-support-results.json').read_text())['rows']
targets={q['baryons']:q for q in inputs if q['settled']}

def integrate(rhoc,gamma,v2,ext=None,tol=1e-7):
 lo=1e-5;hi=1e5;hc=gamma/(gamma-1)*(rhoc/REF)**(gamma-1)
 def rhs(t,y):
  r=np.exp(t);rho=REF*(max(y[0],0)*(gamma-1)/gamma)**(1/(gamma-1))
  gb=0. if ext is None else np.interp(r,ext[0],ext[2])
  return [-(r*gb+G*MSCALE*y[1]/r)/v2,4*np.pi*r**3*rho/MSCALE]
 def surface(t,y):return y[0]
 surface.terminal=True;surface.direction=-1
 sol=solve_ivp(rhs,[np.log(lo),np.log(hi)],[hc,4*np.pi/3*lo**3*rhoc/MSCALE],rtol=tol,atol=tol*1e-3,max_step=.1,events=surface,dense_output=True)
 assert sol.success
 return sol,float(sol.y[1,-1]*MSCALE),bool(len(sol.t_events[0]))

def pressure_fit(gamma,residual=False):
 offsets=[];data={}
 for b,q in targets.items():
  p=q['profile'];r=np.array(p['r_kpc']);rho=np.array(p['rho_Msun_kpc3']);pr=np.array(p['pressure_Msun_kpc3_kms2'])
  if residual:
   pr=pr-np.divide(pr,np.array(p['pressure_ratio']),out=np.zeros_like(pr),where=np.array(p['pressure_ratio'])>0)*(1-np.array(p['assigned_fraction']))**(5/3)
   pr=np.maximum(pr,1e-100)
  at=lambda values,x:np.interp(np.log(x),np.log(r),np.log(np.maximum(values,1e-100)))
  x=np.geomspace(.1,10,100);offsets.extend(at(pr,x)-np.log(REF)-gamma*(at(rho,x)-np.log(REF)))
  data[b]=(r,rho,pr,at)
 v2=float(np.exp(np.mean(offsets)));rows=[]
 for b,(r,rho,pr,at) in data.items():
  errors={}
  for label,lo,hi in [('inner',.1,10),('outer',10,30)]:
   x=np.geomspace(lo,hi,100);err=(np.log(REF*v2)+gamma*(at(rho,x)-np.log(REF))-at(pr,x))/np.log(10)
   errors[label+'_log10_RMS']=float(np.sqrt(np.mean(err*err)))
  slope=np.gradient(np.log(rho),np.log(r));mask=(r>=.1)&(r<=30)
  rows.append(dict(baryons=b,positive_density_slope_fraction=float(np.mean(slope[mask]>0)),**errors))
 return v2,rows

if __name__=='__main__':
 # n=1 Lane-Emden radius for P=K*rho^2, K=v2/REF.
 test,v,finite=integrate(1e8,2,20000,None,1e-9)
 expected=np.pi*np.sqrt((20000/REF)/(2*np.pi*G));assert finite and abs(np.exp(test.t[-1])/expected-1)<1e-6
 expected_mass=4*np.pi**2*(expected/np.pi)**3*1e8;assert abs(v/expected_mass-1)<1e-6
 galaxies={b:Galaxy(b) for b in targets};exts={b:external(g,128) for b,g in galaxies.items()}
 out={'scope':'Pressure-fitted effective EOS, finite spherical mass-normalized equilibria; no formation or stability claim','models':[]}
 for gamma in [4/3,5/3,2.]:
  v2,fits=pressure_fit(gamma);rv2,rfits=pressure_fit(gamma,True);model=dict(gamma=gamma,v0_kms=float(np.sqrt(v2)),pressure_fit= fits,extra_pressure_fit=dict(v0_kms=float(np.sqrt(rv2)),rows=rfits),equilibria=[],rejected_candidates=[],scans=[])
  for b,g in galaxies.items():
   target=g.md[-1];logs=np.linspace(2,13,48);scan=[]
   for lr in logs:
    sol,m,finite=integrate(10**lr,gamma,v2,exts[b]);scan.append(dict(log10_central_density=float(lr),mass=float(m),finite=finite,radius=float(np.exp(sol.t[-1]))))
   model['scans'].append(dict(baryons=b,points=scan));roots=[]
   for left,right in zip(scan[:-1],scan[1:]):
    if left['finite'] and right['finite'] and (left['mass']-target)*(right['mass']-target)<0:
     root=brentq(lambda lr:integrate(10**lr,gamma,v2,exts[b])[1]/target-1,left['log10_central_density'],right['log10_central_density'],xtol=1e-9);roots.append(root)
   for root in roots:
    sol,m,finite=integrate(10**root,gamma,v2,exts[b]);radius=np.exp(sol.t[-1])
    if not finite:
     model['rejected_candidates'].append(dict(baryons=b,central_density=float(10**root),reason='Mass match reaches computational boundary without zero-pressure surface',radius_kpc=float(radius)));continue
    enc=sol.sol(np.log(np.minimum(Robs,radius)))[1]*MSCALE
    vb2=Robs*g.field(Robs,0)[1]-G*g.massD(np.log(Robs))/Robs;pred=np.sqrt(vb2+G*enc/Robs)
    fine,mfine,_=integrate(10**root,gamma,v2,exts[b],1e-9);rfine=np.exp(fine.t[-1]);efine=fine.sol(np.log(np.minimum(Robs,rfine)))[1]*MSCALE;vfine=np.sqrt(vb2+G*efine/Robs)
    assert abs(m/target-1)<1e-6
    model['equilibria'].append(dict(baryons=b,central_density_Msun_kpc3=float(10**root),radius_kpc=float(radius),mass_relative_error=float(m/target-1),refined_speed_max_change=float(np.max(np.abs(vfine-pred))),scores=scores(pred),predicted_kms=pred.tolist()))
  out['models'].append(model);print(json.dumps({k:v for k,v in model.items() if k!='scans'}),flush=True)
 (P/'interaction-support-results.json').write_text(json.dumps(out,allow_nan=False),encoding='utf-8',newline='\n')

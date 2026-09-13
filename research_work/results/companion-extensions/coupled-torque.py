"""Mass- and scalar-angular-momentum-conserving circular endpoint solver."""
from settling import *
selected=json.loads((P/'settling-results.json').read_text())['selected']

def ordinary_field(g):
 # Evaluate positive ordinary-matter sources directly, avoiding subtraction
 # of nearly equal total and companion forces below the legacy radial clamp.
 r=np.geomspace(1e-6,1e6,32768);mu,ww=leggauss(256);mu=(mu+1)/2;ww=ww/2
 R=r[:,None]*np.sqrt(1-mu**2);z=r[:,None]*mu
 lr=np.log(np.clip(r,g.r[1],g.r[-1]));mg=g.massG(lr)
 pb=-g.bindG(lr);pb[r>g.r[-1]]=-G*mg[r>g.r[-1]]/r[r>g.r[-1]]
 gb=G*mg/r**2
 for M,a,b in components[g.variant]:
  B=np.sqrt(z*z+b*b);d=R*R+(a+B)**2
  pb-=G*M*(1/np.sqrt(d))@ww
  gb+=G*M*((R*np.sqrt(1-mu**2)+(a+B)*z/B*mu)/d**1.5)@ww
 assert np.all(np.diff(r*r*gb/G)>=0)
 return r,pb,gb

def knots(r,w):
 ix=np.argsort(r);rr=r[ix];ww=w[ix];unique,first=np.unique(rr,return_index=True);mass=np.add.reduceat(ww,first)
 return unique,np.cumsum(mass)-.5*mass

def circular_curve(r,w,ext,mesh=16384):
 # Smooth packet self-force over a fixed logarithmic interpolation mesh.
 # This avoids artificial angular-momentum gaps at nearly coincident thin shells.
 rr=np.geomspace(1e-6,1e6,mesh)
 coord=(np.log(r)-np.log(rr[0]))/np.log(rr[1]/rr[0]);ix=np.floor(coord).astype(int);t=coord-ix
 assert np.all((ix>=0)&(ix<mesh-1))
 nodes=np.bincount(ix,weights=w*(1-t),minlength=mesh)+np.bincount(ix+1,weights=w*t,minlength=mesh)
 mm=np.cumsum(nodes)-.5*nodes
 mb=np.interp(rr,ext[0],ext[0]**2*ext[2]/G)
 jb=np.sqrt(G*rr*(mb+mm))
 assert np.all(np.diff(jb)>0)
 return rr,mm,jb

def prepare(b,n):
 g=Galaxy(b);ext=ordinary_field(g);data=setup(g,n,256);edges,r,dm,sigma,rho,vb=data
 m=selected['mass_eV']*electron_volt/c**2;nn=rho*MSUN/KPC**3/m;rat=m*m*sigma*1e6/(2*np.pi*hbar*hbar)*(zeta(1.5)/nn)**(2/3);f=np.maximum(0,1-np.minimum(rat,1)**1.5)
 pos=np.r_[r,r];weight=np.r_[dm*(1-f),dm*f];donor=np.r_[np.zeros(n,dtype=bool),np.ones(n,dtype=bool)];valid=weight>dm.sum()*1e-18
 return g,ext,data,pos[valid],weight[valid],donor[valid]

def solve(inp,lo,hi,s,seed=0.,mesh=16384):
 g,ext,data,r0,w,donor=inp;recv=(~donor)&(r0>=lo)&(r0<=hi);rr,mm,jcurve=circular_curve(r0,w,ext,mesh);j0=np.interp(r0,rr,jcurve)
 r=r0.copy();r[donor]*=s;r[recv]*=1+seed;res=1.;delta=0.
 for it in range(2000):
  rr,mm,jcurve=circular_curve(r,w,ext,mesh);jc=np.interp(r,rr,jcurve);loss=np.sum(w[donor]*(j0[donor]-jc[donor]));delta=loss/w[recv].sum()
  target=j0.copy();target[recv]+=delta;target[donor]=jc[donor]
  new=np.interp(target,jcurve,rr);new[donor]=s*r0[donor]
  res=float(np.max(np.abs(np.log(new/r))))
  if res<1e-8:break
  r=np.exp(.75*np.log(r)+.25*np.log(new))
 rr,mm,jcurve=circular_curve(r,w,ext,mesh);jc=np.interp(r,rr,jcurve);loss=np.sum(w[donor]*(j0[donor]-jc[donor]));gain=np.sum(w[recv]*(jc[recv]-j0[recv]));untouched=(~donor)&(~recv)
 def speeds(pos):
  x,middle,_=circular_curve(pos,w,ext,mesh);enclosed=np.interp(Robs,x,middle);return np.sqrt(data[-1]+G*enclosed/Robs)
 initial=mechanical(r0,w,ext);final=mechanical(r,w,ext)
 # Independent kinetic energy from assigned circular angular momenta.
 # Retain the virial estimate to quantify shell/grid force inconsistency.
 virialQ=initial['total']-final['total']
 for en,rad,j in [(initial,r0,j0),(final,r,jc)]:
  en['virial_kinetic']=en['support_kinetic'];en['support_kinetic']=float(.5*np.sum(w*(j/rad)**2))
  en['total']=en['self_energy']+en['external_energy']+en['support_kinetic']
 Q=initial['total']-final['total']
 imposed=r0.copy();imposed[donor]*=s
 out=dict(receiver_region_kpc=[lo,hi],converged=bool(res<1e-8),iterations=it+1,radius_residual=res,angular_balance_relative_to_transfer=float((gain-loss)/max(abs(loss),1)),untorqued_max_relative_j_change=float(np.max(np.abs(jc[untouched]/j0[untouched]-1))),receiver_mean_radius_change=float(np.average(r[recv]/r0[recv]-1,weights=w[recv])),untorqued_mean_radius_change=float(np.average(r[untouched]/r0[untouched]-1,weights=w[untouched])),released_energy_J=float(Q*MSUN*1e6),energy_ledger_relative_residual=float((final['total']+Q-initial['total'])/abs(initial['total'])),positive_release=bool(Q>=0),scores=scores(speeds(r)),control_scores=scores(speeds(r0)),imposed_scores=scores(speeds(imposed)),predicted_kms=speeds(r).tolist(),initial_energy=initial,final_energy=final)
 out['virial_vs_circular_release_relative_difference']=float((virialQ-Q)/max(abs(Q),1))
 out['mass_Msun']=float(w.sum())
 if out['converged']:
  assert abs(out['angular_balance_relative_to_transfer'])<1e-5
  assert out['untorqued_max_relative_j_change']<1e-6
 return out,r

if __name__=='__main__':
 out={'scope':'Self-consistent circular endpoint construction, orientation-matched scalar transfer; not evolution or stability','rows':[],'checks':[]}
 for b in ['I','II']:
  inp=prepare(b,8192);fineinp=prepare(b,16384)
  zero,zr=solve(inp,30,60,1.);assert np.max(abs(zr/inp[3]-1))<1e-7
  assert abs(zero['released_energy_J'])<1e-10*MSUN*1e6*abs(zero['initial_energy']['total'])
  for lo,hi in [(15,30),(30,60),(60,120)]:
   q,r=solve(inp,lo,hi,selected['s']);q['baryons']=b;out['rows'].append(q)
   alt,ar=solve(inp,lo,hi,selected['s'],.05);fine,fr=solve(fineinp,lo,hi,selected['s'],mesh=32768)
   assert all(row['converged'] for row in [q,alt,fine])
   assert all(row['positive_release'] for row in [q,alt,fine])
   assert np.max(abs(np.array(q['predicted_kms'])-fine['predicted_kms']))<.1
   out['checks'].append(dict(baryons=b,region=[lo,hi],alternate_seed_converged=alt['converged'],alternate_seed_max_speed_difference=float(np.max(abs(np.array(q['predicted_kms'])-alt['predicted_kms']))),refined_converged=fine['converged'],refined_max_speed_difference=float(np.max(abs(np.array(q['predicted_kms'])-fine['predicted_kms']))),refined_energy_relative_change=float(fine['released_energy_J']/q['released_energy_J']-1),refined_scores=fine['scores']))
   print(json.dumps({k:v for k,v in q.items() if k not in ['predicted_kms','initial_energy','final_energy']}),flush=True)
 out['configuration']=dict(source_bins=8192,refined_source_bins=16384,field_nodes=16384,refined_field_nodes=32768,mass_eV=selected['mass_eV'],s=selected['s'],zero_contraction_checks_passed=True)
 (P/'coupled-torque-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

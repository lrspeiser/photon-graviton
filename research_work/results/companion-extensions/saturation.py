"""Mass-conserving ordered donor packing with a shared saturation threshold."""
from stopping import ct,np,json,P
from scipy.optimize import brentq

def pack(inp,cap):
 _,_,data,r0,w,donor=inp;edges,r=data[:2]
 ix=np.searchsorted(r,r0[donor]);mass=np.zeros(len(r));mass[ix]=w[donor]
 volume=4*np.pi/3*np.diff(edges**3)
 shrunk=np.where(mass>0,np.minimum(volume,mass/cap),volume)
 before=np.r_[0,np.cumsum(shrunk[:-1])]+4*np.pi/3*edges[0]**3
 frac=(r**3-edges[:-1]**3)/np.diff(edges**3)
 final=(3/(4*np.pi)*(before+frac*shrunk))**(1/3)
 ratio=final[ix]/r0[donor]
 assert np.all(ratio<=1+1e-10) and np.all(np.diff(final)>0)
 expected=np.maximum(mass/volume,cap)
 occupied=mass>0
 assert np.max(abs((mass[occupied]/shrunk[occupied])/expected[occupied]-1))<1e-12
 return ratio,dict(mean_contraction=float(np.average(ratio,weights=w[donor])),min_contraction=float(ratio.min()),max_contraction=float(ratio.max()),donor_mass_Msun=float(mass.sum()),initially_above_threshold_mass_fraction=float(mass[mass/volume>=cap].sum()/mass.sum()))

if __name__=='__main__':
 inp=ct.prepare('I',8192);target=ct.selected['s']
 cap=10**brentq(lambda logcap:pack(inp,10**logcap)[1]['mean_contraction']-target,-6,14,xtol=1e-10)
 out=dict(saturation_density_Msun_kpc3=cap,calibration='Baseline I mean donor radius ratio; previous fitted contraction target, no new speed fit',target=target,rows=[],checks=[])
 for b in ['I','II']:
  if b!='I':inp=ct.prepare(b,8192)
  fineinp=ct.prepare(b,16384);s,summary=pack(inp,cap);fs,fsummary=pack(fineinp,cap)
  for lo,hi in [(15,30),(30,60),(60,120)]:
   q,r=ct.solve(inp,lo,hi,s);q.update(baryons=b,packing=summary)
   fq,fr=ct.solve(fineinp,lo,hi,fs,mesh=32768)
   out['rows'].append(q)
   out['checks'].append(dict(baryons=b,receiver_band=[lo,hi],refined_converged=fq['converged'],max_speed_difference_kms=float(np.max(abs(np.array(q['predicted_kms'])-fq['predicted_kms']))),release_relative_difference=float(fq['released_energy_J']/q['released_energy_J']-1),refined_packing=fsummary))
   print(json.dumps({k:v for k,v in q.items() if k not in ['predicted_kms','initial_energy','final_energy']}),flush=True)
 (P/'saturation-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

"""Conditional direct-ray and internal-torque energy requirements."""
from settling import *
import importlib.util
spec=importlib.util.spec_from_file_location('orbital',P/'orbital-support.py');orb=importlib.util.module_from_spec(spec);spec.loader.exec_module(orb)

def run(b,n=4096):
 g=Galaxy(b);ext=external(g,256);inp=json.loads((P/'local-support-results.json').read_text())['rows'];fields={}
 for settled in [False,True]:
  q=next(q for q in inp if q['baryons']==b and q['settled']==settled);_,fields[settled]=orb.analyze(q,g,ext)
 edges,r,dm,sigma,rho,vb=setup(g,n,256);mass=orb.local.selected['mass_eV'];s=orb.local.selected['s'];m=mass*electron_volt/c**2
 nn=rho*MSUN/KPC**3/m;rat=m*m*sigma*1e6/(2*np.pi*hbar*hbar)*(zeta(1.5)/nn)**(2/3);f=np.maximum(0,1-np.minimum(rat,1)**1.5)
 interp=lambda field,k,x:np.interp(x,field[0],field[k])
 gi=interp(fields[False],1,r);gf=interp(fields[True],1,s*r);ji=np.sqrt(r**3*gi);jf=np.sqrt((s*r)**3*gf)
 mask=(r>=.1)&(r<=30)&(f>0);dJ=dm[mask]*f[mask]*(ji[mask]-jf[mask]);assert np.all(dJ>0)
 J=float(dJ.sum());emin=float(np.sum(c/1000*dJ/r[mask]));cool=ledger(setup(g,4096,256),ext,mass,s)['cooling_energy_Msun_kms2']
 receivers=[];rf,af,psif=fields[True];jcurve=np.sqrt(rf**3*af);assert np.all(np.diff(jcurve)>0)
 for lo,hi in [(15,30),(30,60),(60,120)]:
  use=(r>=lo)&(r<=hi);w=dm[use]*(1-f[use]);rr=r[use];mr=w.sum();dj=J/mr
  a=interp(fields[True],1,rr);j=np.sqrt(rr**3*a);dest=np.interp(j+dj,jcurve,rf);anew=interp(fields[True],1,dest)
  before=.5*rr*a-interp(fields[True],2,rr);after=.5*dest*anew-interp(fields[True],2,dest)
  first=float(np.sum(w*np.sqrt(a/rr)*dj));finite=float(np.sum(w*(after-before)))
  receivers.append(dict(region_kpc=[lo,hi],receiver_mass_Msun=float(mr),added_specific_angular_momentum_kpc_kms=float(dj),mean_fractional_j_increase=float(np.average(dj/j,weights=w)),mean_fractional_radius_increase=float(np.average(dest/rr-1,weights=w)),first_order_work_over_cooling=first/cool,frozen_potential_energy_gain_over_cooling=finite/cool,angular_momentum_relative_residual=float((np.sum(w*dj)-J)/J)))
 return dict(baryons=b,moved_source_region_kpc=[.1,30],sum_angular_momentum_magnitudes_Msun_kpc_kms=J,minimum_direct_ray_energy_J=emin*MSUN*1e6,minimum_direct_ray_energy_over_cooling=emin/cool,minimum_direct_ray_energy_over_total_deposit_rest_energy=emin/(g.md[-1]*(c/1000)**2),cooling_energy_J=cool*MSUN*1e6,receivers=receivers)

if __name__=='__main__':
 out={'scope':'Direct orbital-ray emission bound and orientation-matched frozen-potential receiver scenarios; not full evolution','rows':[],'refinement':[]}
 for b in ['I','II']:
  q=run(b);fine=run(b,8192);out['rows'].append(q);out['refinement'].append(dict(baryons=b,direct_energy_relative_change=fine['minimum_direct_ray_energy_J']/q['minimum_direct_ray_energy_J']-1,max_receiver_energy_ratio_change=max(abs(a['frozen_potential_energy_gain_over_cooling']-z['frozen_potential_energy_gain_over_cooling']) for a,z in zip(q['receivers'],fine['receivers']))))
  assert max(abs(x['angular_momentum_relative_residual']) for x in q['receivers'])<1e-12
 (P/'torque-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n');print(json.dumps(out,indent=2))

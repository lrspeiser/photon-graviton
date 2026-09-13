"""Isotropic necessary-positivity test and circular-shell existence construction."""
from settling import *
import importlib.util
spec=importlib.util.spec_from_file_location('local_support',P/'local-support.py');local=importlib.util.module_from_spec(spec);spec.loader.exec_module(local)

def analyze(q,g,ext):
 p=q['profile'];r=np.array(p['r_kpc']);rho=np.array(p['rho_Msun_kpc3'])
 mass=4*np.pi*cumulative_trapezoid(r*r*rho,r,initial=0);acc=np.interp(r,ext[0],ext[2])+G*mass/r**2
 outer=-cumulative_trapezoid((4*np.pi*r*rho)[::-1],r[::-1],initial=0)[::-1]
 psi=-np.interp(r,ext[0],ext[1])+G*(mass/r+outer)
 slope=np.gradient(np.log(rho),np.log(r));epicycle=np.gradient(acc,r,edge_order=2)+3*acc/r
 circ2=r*acc;energy=.5*circ2-psi;mask=(r>=.1)&(r<=30);up=mask&(slope>1e-3)
 starts=np.flatnonzero(up&~np.r_[False,up[:-1]]);ends=np.flatnonzero(up&~np.r_[up[1:],False])
 intervals=[dict(inner_kpc=float(r[i]),outer_kpc=float(r[j]),max_log_slope=float(slope[i:j+1].max())) for i,j in zip(starts,ends)]
 # Spherical Jeans with zero radial stress and total tangential dispersion v_c^2.
 assert np.max(np.abs(circ2/r-acc)/acc)<1e-12
 samples=[dict(r_kpc=x,orbital_speed_kms=float(np.sqrt(np.interp(x,r,circ2))),angular_momentum_kpc_kms=float(x*np.sqrt(np.interp(x,r,circ2))),epicyclic_over_orbital_frequency_squared=float(np.interp(x,r,epicycle/(acc/r)) )) for x in [1.,3.,5.,8.,12.,20.,30.]]
 return dict(baryons=q['baryons'],settled=q['settled'],isotropic_necessary_condition_passes=not bool(np.any(up)),positive_density_slope_intervals=intervals,min_epicyclic_frequency_squared=float(epicycle[mask].min()),max_specific_orbital_energy_kms2=float(energy[mask].max()),all_tested_circular_orbits_bound=bool(np.all(energy[mask]<0)),all_tested_circular_orbits_radially_stable=bool(np.all(epicycle[mask]>0)),samples=samples), (r,acc,psi)

if __name__=='__main__':
 inp=json.loads((P/'local-support-results.json').read_text())['rows'];out={'scope':'Necessary isotropic positivity and singular circular-shell construction; not collective stability or formation','rows':[],'refinement':[],'migration':[]}
 for b in ['I','II']:
  g=Galaxy(b);ext=external(g,256);fields={}
  for settled in [False,True]:
   q=next(v for v in inp if v['baryons']==b and v['settled']==settled);result,field=analyze(q,g,ext);out['rows'].append(result);fields[settled]=field
   original_setup=local.setup
   local.setup=lambda galaxy,n=2048,nmu=192:original_setup(galaxy,8192,384)
   fine,_=analyze(local.run(b,settled,16384),g,ext)
   local.setup=original_setup
   out['refinement'].append(dict(baryons=b,settled=settled,necessary_condition_same=fine['isotropic_necessary_condition_passes']==result['isotropic_necessary_condition_passes'],refined_intervals=fine['positive_density_slope_intervals'],max_sample_speed_difference=float(max(abs(x['orbital_speed_kms']-y['orbital_speed_kms']) for x,y in zip(result['samples'],fine['samples'])))))
  edges,rs,dm,sigma,rho,vb=setup(g,4096,256);m=local.selected['mass_eV']*electron_volt/c**2;nn=rho*MSUN/KPC**3/m
  rat=m*m*sigma*1e6/(2*np.pi*hbar*hbar)*(zeta(1.5)/nn)**(2/3);f=np.maximum(0,1-np.minimum(rat,1)**1.5);s=local.selected['s'];dest=s*rs
  gi=np.interp(rs,fields[False][0],fields[False][1]);gf=np.interp(dest,fields[True][0],fields[True][1]);ji=np.sqrt(rs**3*gi);jf=np.sqrt(dest**3*gf)
  mask=(rs>=.1)&(rs<=30)&(f>0);w=dm[mask]*f[mask];loss=1-jf[mask]/ji[mask]
  ei=.5*rs*gi-np.interp(rs,fields[False][0],fields[False][2]);ef=.5*dest*gf-np.interp(dest,fields[True][0],fields[True][2])
  out['migration'].append(dict(baryons=b,source_region_kpc=[.1,30],mobile_mass_weighted_fractional_specific_angular_momentum_reduction=float(np.average(loss,weights=w)),minimum_reduction=float(loss.min()),maximum_reduction=float(loss.max()),fraction_mobile_mass_needing_reduction=float(np.sum(w[loss>0])/w.sum()),single_orbit_specific_energy_decrease_kms2=float(np.average((ei-ef)[mask],weights=w)),energy_note='Single-orbit energies in different self-consistent potentials; not an additive global cooling budget'))
 # Point mass has kappa^2/Omega^2=1; uniform-density core has ratio 4.
 rt=np.geomspace(.1,30,10000)
 for force,expected in [(rt**-2,1.),(rt,4.)]:
  ratio=(np.gradient(force,rt,edge_order=2)+3*force/rt)/(force/rt)
  assert np.max(abs(ratio[2:-2]-expected))<1e-5
 (P/'orbital-support-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps(out,indent=2))

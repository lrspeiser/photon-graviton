"""Local isotropic support and ideal-Bose pressure ceiling; no speed refit."""
from settling import *
selected=json.loads((P/'settling-results.json').read_text())['selected']

def run(b,settled,n=8192,boundary=1e5):
 g=Galaxy(b);data=setup(g,4096,256);edges,rs,dm,sigma2,rhos,vb=data
 mev=selected['mass_eV'];s=selected['s'] if settled else 1.;m=mev*electron_volt/c**2
 ns=rhos*MSUN/KPC**3/m
 rat=m*m*sigma2*1e6/(2*np.pi*hbar*hbar)*(zeta(1.5)/ns)**(2/3)
 fs=np.maximum(0,1-np.minimum(rat,1)**1.5)
 r=np.geomspace(1e-4,boundary,n)
 def rho0(q):return np.exp(g.logrho(np.log(np.clip(q,g.r[1],g.r[-1]))))
 def f0(q):return np.interp(q,rs,fs,left=fs[0],right=0)
 comp=f0(r/s)*rho0(r/s)/s**3 if settled else np.zeros(n)
 rho=(1-f0(r))*rho0(r)+comp if settled else rho0(r)
 mass=4*np.pi*cumulative_trapezoid(rho*r*r,r,initial=0)
 ext=external(g,128);acc=np.interp(r,ext[0],ext[2])+G*mass/r**2
 integrand=rho*acc
 pressure=-cumulative_trapezoid(integrand[::-1],r[::-1],initial=0)[::-1]
 rhoSI=rho*MSUN/KPC**3;nnum=rhoSI/m
 kTc=2*np.pi*hbar*hbar/m*(nnum/zeta(1.5))**(2/3)
 pcSI=zeta(2.5)/zeta(1.5)*nnum*kTc
 ratio=pressure*MSUN/KPC**3*1e6/pcSI
 fthermal=np.maximum(0,1-ratio**.6)
 assigned=comp/rho
 weight=4*np.pi*r*r*rho;mask=(r>=.1)&(r<=30)
 avg=lambda a:float(np.trapezoid((a*weight)[mask],r[mask])/np.trapezoid(weight[mask],r[mask]))
 deriv=np.gradient(pressure,r,edge_order=2);res=np.abs(deriv+integrand)/integrand
 samples=[]
 for radius in [1.,3.,5.,8.,12.,20.,30.]:
  samples.append(dict(r_kpc=radius,required_sigma_kms=float(np.sqrt(np.interp(radius,r,pressure/rho))),pressure_over_condensation_ceiling=float(np.interp(radius,r,ratio)),assigned_compact_fraction=float(np.interp(radius,r,assigned)),pressure_implied_condensed_fraction=float(np.interp(radius,r,fthermal))))
  entry=samples[-1];pressure_fraction=(1-entry['assigned_compact_fraction'])**(5/3)/entry['pressure_over_condensation_ceiling']
  entry['extra_support_fraction_at_assigned_phase']=max(0.,1-pressure_fraction)
  entry['effective_mass_eV_for_assigned_fraction_at_fixed_profile']=mev*pressure_fraction**(3/8)
 return dict(baryons=b,settled=settled,mass_eV=mev,s=s,region_kpc=[.1,30],mass_weighted_fraction_above_ceiling=avg((ratio>1).astype(float)),assigned_fraction=avg(assigned),pressure_implied_fraction=avg(fthermal),mean_absolute_fraction_difference=avg(np.abs(assigned-fthermal)),max_hydrostatic_relative_residual=float(res[mask].max()),integrated_mass_Msun=float(mass[-1]),mass_relative_to_original=float(mass[-1]/g.md[-1]),samples=samples,profile=dict(r_kpc=r.tolist(),rho_Msun_kpc3=rho.tolist(),pressure_Msun_kpc3_kms2=pressure.tolist(),pressure_ratio=ratio.tolist(),assigned_fraction=assigned.tolist(),implied_fraction=fthermal.tolist()))

if __name__=='__main__':
 # Independent Plummer hydrostatic integral, in G=M=a=1 units.
 rt=np.geomspace(1e-4,1e5,16384);rh=3/(4*np.pi)*(1+rt*rt)**(-2.5);gt=rt/(1+rt*rt)**1.5
 pt=-cumulative_trapezoid((rh*gt)[::-1],rt[::-1],initial=0)[::-1]
 exact=(1+rt*rt)**(-3)/(8*np.pi);testmask=(rt>.01)&(rt<100)
 assert np.max(np.abs(pt[testmask]/exact[testmask]-1))<1e-5
 out={'scope':'Ideal noninteracting Bose support ceiling with spherical force average; not a phase equilibrium or stability solution','rows':[],'refinement':[]}
 for b in ['I','II']:
  for settled in [False,True]:
   q=run(b,settled);out['rows'].append(q)
   if settled:
    for n,boundary,kind in [(16384,1e5,'radial'),(8192,5e4,'boundary')]:
     fine=run(b,settled,n,boundary)
     out['refinement'].append(dict(baryons=b,kind=kind,max_sample_pressure_ratio_change=float(max(abs(x['pressure_over_condensation_ceiling']/y['pressure_over_condensation_ceiling']-1) for x,y in zip(fine['samples'],q['samples']))),fraction_above_ceiling_change=fine['mass_weighted_fraction_above_ceiling']-q['mass_weighted_fraction_above_ceiling']))
   assert q['max_hydrostatic_relative_residual']<.01
   print(json.dumps({k:v for k,v in q.items() if k!='profile'}),flush=True)
 (P/'local-support-results.json').write_text(json.dumps(out,allow_nan=False),encoding='utf-8',newline='\n')

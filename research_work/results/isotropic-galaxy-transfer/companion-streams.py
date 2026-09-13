"""Near/far companion streams: geometry, cached catalog, and finite supply shells."""
import json,hashlib,ast
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2]
alpha=.0002488993286382367
pilot=json.loads((P.parent/'cluster-catalog-pilot/results.json').read_text())
source=P.parent/'cluster-catalog-pilot/run.py';tree=ast.parse(source.read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='read'],type_ignores=[]),str(source),'exec'))
def frac(d):return -np.expm1(-alpha*np.asarray(d))
def incident(d):return frac(d)/np.asarray(d)**2
geometry=[]
for R in [.01,.1,1.]:
    for q in [1.1,2,5,10,20,100,1000,10000]:
        D=q*R
        geometry.append(dict(receiver_radius_Mpc=R,source_distance_Mpc=D,D_over_R=q,converted_fraction=float(frac(D)),near_far_incident_ratio=float(incident(D-R)/incident(D+R)),angular_half_width_degrees=float(np.degrees(np.arcsin(R/D)))))
catalogs=[]
for run in pilot['runs']:
    file=ROOT/f"companion_causal_test/data/{run['model']}.dat"
    assert hashlib.sha256(file.read_bytes()).hexdigest()==run['input_sha256']
    sources,_=read(file);center=next(v for v in sources if v['name']=='NGC4486')['xyz']
    rows=[]
    for v in sources:
        vec=v['xyz']-center;D=float(np.linalg.norm(vec))
        if D<=1:continue
        direction=vec/D;flux=v['lum']*incident(D)/(4*np.pi)
        old=next(z for z in run['sources'] if z['name']==v['name'])
        rows.append(dict(name=v['name'],distance_Mpc=D,direction=direction.tolist(),converted_fraction=float(frac(D)),incident_companion_flux_Lsun_Mpc2=float(flux),prior_absorbed_power_Lsun=old['deposited_power_Lsun']))
    total=sum(r['incident_companion_flux_Lsun_Mpc2'] for r in rows);power=sum(r['prior_absorbed_power_Lsun'] for r in rows)
    weights=np.array([v['incident_companion_flux_Lsun_Mpc2']/total for v in rows]);dirs=np.array([v['direction'] for v in rows]);vector=np.sum(weights[:,None]*dirs,axis=0);tensor=np.einsum('i,ij,ik->jk',weights,dirs,dirs)
    bins=[]
    for lo,hi in [(1,2),(2,5),(5,10),(10,20),(20,50),(50,100),(100,1e9)]:
        group=[r for r in rows if lo<r['distance_Mpc']<=hi]
        bins.append(dict(lower_Mpc=lo,upper_Mpc=hi if hi<1e9 else None,n=len(group),incident_share=sum(r['incident_companion_flux_Lsun_Mpc2'] for r in group)/total,prior_absorbed_share=sum(r['prior_absorbed_power_Lsun'] for r in group)/power))
    ordered=sorted(rows,key=lambda r:r['incident_companion_flux_Lsun_Mpc2'],reverse=True)
    catalogs.append(dict(model=run['model'],external_sources=len(rows),incident_total_Lsun_Mpc2=total,prior_absorbed_total_Lsun=power,net_direction_fraction=float(np.linalg.norm(vector)),angular_second_moment_eigenvalues=np.linalg.eigvalsh(tensor).tolist(),top10_incident_share=sum(r['incident_companion_flux_Lsun_Mpc2'] for r in ordered[:10])/total,bins=bins,top_sources=ordered[:20],maximum_catalog_distance_Mpc=max(r['distance_Mpc'] for r in rows)))
# Homogeneous hypothetical emissivity shell integral, per unit luminosity density.
# Exact integral of conversion fraction from zero to D, no fixed universe cutoff.
shells=[]
for D in [1,10,100,1000,3000,10000]:
    integral=D+np.expm1(-alpha*D)/alpha
    shells.append(dict(Dmax_Mpc=D,integrated_incident_per_emissivity_Mpc=float(integral),converted_fraction_at_edge=float(frac(D))))
assert all(abs(sum(b['incident_share'] for b in c['bins'])-1)<1e-12 for c in catalogs)
assert all(abs(sum(c['angular_second_moment_eigenvalues'])-1)<1e-12 for c in catalogs)
opposite=np.array([[1.,0,0],[-1.,0,0]])
assert np.linalg.norm(opposite.sum(axis=0))==0 and np.linalg.norm(opposite,axis=1).sum()==2
focus=[]
for v in json.loads((P/'third-retention-optics-results.json').read_text())['rows']:
    if v['model']!='attenuated_Chabrier':continue
    bend=v['lens_catalog_arcsec']/v['geometry']['Dls_over_Ds']
    focus.append(dict(Name=v['Name'],required_total_bend_arcsec=bend,required_total_bend_radians=bend/206264.80624709636))
out=dict(alpha_per_Mpc=alpha,scope='Source-flow diagnostic; no measured historical supply or predicted halo',geometry=geometry,catalogs=catalogs,homogeneous_shells=shells,lightlike_deflection_scale=focus,input_sha256={str(file.relative_to(ROOT)):hashlib.sha256(file.read_bytes()).hexdigest() for file in [source,P.parent/'cluster-catalog-pilot/results.json']})
(P/'companion-streams-results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
q=np.geomspace(1.05,1000,400)
for R in [.01,.1,1]:ax[0].loglog(q,incident((q-1)*R)/incident((q+1)*R),label=f'R={R:g} Mpc')
ax[0].set(xlabel='Source distance / receiver radius',ylabel='Companion intensity: near side / far side',title='Nearby streams vary across the receiver');ax[0].legend();ax[0].grid(alpha=.2)
labels=['1–2','2–5','5–10','10–20','20–50','50–100','>100']
for j,c in enumerate(catalogs):ax[1].bar(np.arange(7)+(j-.5)*.36,[100*b['incident_share'] for b in c['bins']],width=.36,label=c['model'])
ax[1].set(xticks=np.arange(7),xticklabels=labels,xlabel='Source distance from M87 (Mpc)',ylabel='Share of catalog external companion intensity (%)',title='Cached nearby-galaxy catalog only');ax[1].legend()
fig.savefig(P/'companion-streams.png',dpi=160);plt.close(fig)
for c in catalogs:print(json.dumps({k:v for k,v in c.items() if k not in ['top_sources']},indent=2))

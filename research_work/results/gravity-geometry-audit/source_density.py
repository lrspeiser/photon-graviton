"""Conditional Newtonian source reconstruction, not a photon energy budget."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator
H=Path(__file__).resolve().parent;ROOT=H.parents[2]
sys.path.insert(0,str(H.parent/'rotating-bar-orbits'))
from run import DiskGrid,Field
from field import bar,nuclei
G=4.30091727003628e-6;KPC=3.085677581491367e19
fit=H.parent/'joint-galaxy-audit/results.json'
pars=json.loads(fit.read_text())['sparc']['parameters']
A,p,astar=pars['A'],pars['p'],pars['a_star_m_s2']*KPC/1e6
disk=DiskGrid(256,129,None)
angles=np.arange(64)*2*np.pi/64
def mean_field(field,R,z):
    R,z=np.broadcast_arrays(np.atleast_1d(R),np.atleast_1d(z))
    xyz=np.stack([R[:,None]*np.cos(angles),R[:,None]*np.sin(angles),np.broadcast_to(z[:,None],(len(R),len(angles)))],axis=-1)
    pot,a=field.evaluate(xyz.reshape(-1,3));a=a.reshape(-1,len(angles),3)
    return pot.reshape(-1,len(angles)).mean(1),np.c_[(a[:,:,0]*np.cos(angles)+a[:,:,1]*np.sin(angles)).mean(1),a[:,:,2].mean(1)]
rr,zz=np.meshgrid([.5,1,1.5,2,3,5,8,12,20],[0,.1,.3,.7,1.1,2,4],indexing='ij')
rr=rr.ravel();zz=zz.ravel();rows=[]
xyz=np.stack([rr[:,None]*np.cos(angles),rr[:,None]*np.sin(angles),np.broadcast_to(zz[:,None],(len(rr),len(angles)))],axis=-1)
declared=(bar(xyz.reshape(-1,3))+nuclei(xyz.reshape(-1,3))).reshape(-1,len(angles)).mean(1)
for sigma,rd,h,rhole,kind in [(1.332e9,2.,.3,2.7,'exp'),(8.97e8,2.8,.9,2.7,'exp'),(5.81e7,7.,.085,4.,'sech2'),(2.68e9,1.5,.045,12.,'sech2')]:
    sig=sigma*np.exp(-rhole/rr-rr/rd)
    vertical=np.exp(-abs(zz)/h)/(2*h) if kind=='exp' else np.exp(-2*np.logaddexp(zz/(2*h),-zz/(2*h))+2*np.log(2))/(4*h)
    declared+=sig*vertical
# Existing central point_potential uses a Plummer sphere with 0.001-kpc softening.
declared+=3*4.1e6*.001**2/(4*np.pi*(rr*rr+zz*zz+.001**2)**2.5)
declared/=1e9
for order in [40,64]:
    f=Field(disk,order=order)
    rgrid=np.geomspace(.2,32,960);potential,acc=mean_field(f,rgrid,0);gb=-acc[:,0]
    assert np.all(gb>0) and np.all(np.diff(potential)>0)
    gc=A*astar*(gb/astar)**p
    radial=PchipInterpolator(rgrid,gc,extrapolate=False)
    ratio=PchipInterpolator(potential,gc/gb,extrapolate=False)
    ph,a=mean_field(f,rr,zz);rad=np.hypot(rr,zz)
    spherical=(2*radial(rad)/rad+radial.derivative()(rad))/(4*np.pi*G)/1e9
    for h in [.004,.002,.001]:
        pr,ar=mean_field(f,rr+h,zz);mr,br=mean_field(f,rr-h,zz)
        pz,az=mean_field(f,rr,zz+h);mz,bz=mean_field(f,rr,zz-h)
        rhob=-((ar[:,0]-br[:,0])/(2*h)+a[:,0]/rr+(az[:,1]-bz[:,1])/(2*h))/(4*np.pi*G)/1e9
        # Exact chain rule for F(Phi_b); numerical baryonic divergence only.
        shape=ratio(ph)*rhob+ratio.derivative()(ph)*np.sum(a*a,axis=1)/(4*np.pi*G)/1e9
        shape_source_substitution=shape+ratio(ph)*(declared-rhob)
        # Independent divergence of the constructed companion acceleration.
        acr=ratio(pr)[:,None]*ar;bcr=ratio(mr)[:,None]*br
        acz=ratio(pz)[:,None]*az;bcz=ratio(mz)[:,None]*bz
        direct=-((acr[:,0]-bcr[:,0])/(2*h)+ratio(ph)*a[:,0]/rr+(acz[:,1]-bcz[:,1])/(2*h))/(4*np.pi*G)/1e9
        assert np.isfinite(shape).all() and np.isfinite(spherical).all()
        for i in range(len(rr)):
            rows.append(dict(bar_order=order,h_kpc=h,R_kpc=rr[i],z_kpc=zz[i],baryon_density_Msun_pc3=rhob[i],
                spherical_density_Msun_pc3=spherical[i],shape_density_Msun_pc3=shape[i],shape_divergence_density_Msun_pc3=direct[i]))
            rows[-1].update(declared_baryon_density_Msun_pc3=declared[i],shape_with_declared_source_Msun_pc3=shape_source_substitution[i])
    print('Completed source reconstruction L'+str(order),flush=True)
final=[r for r in rows if r['bar_order']==64 and r['h_kpc']==.001]
summary={'probe_count':len(final),'parameters':pars,'models':{}}
for model in ['baryon','spherical','shape']:
    key=model+'_density_Msun_pc3'
    vals=np.array([r[key] for r in final])
    summary['models'][model]={'minimum_density_Msun_pc3':float(vals.min()),'maximum_density_Msun_pc3':float(vals.max()),
        'negative_probe_count':int((vals<0).sum()),'minimum_probe':final[int(vals.argmin())]}
summary['asymptotic_if_extended_around_finite_baryonic_mass']={'g_exponent':-2*p,'enclosed_source_mass_exponent':2-2*p,'potential_growth_exponent':1-2*p,
    'finite_potential_at_infinity':bool(p>.5),'finite_total_equivalent_mass':False}
summary['scope']='Conditional ordinary positive-source Poisson interpretation. No global positivity proof, photon/graviton identity, capture law, energy normalization or observational fit. No heldout score evaluated.'
summary['fit_sha256']=hashlib.sha256(fit.read_bytes()).hexdigest()
(H/'source-density-rows.json').write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n',encoding='utf-8')
(H/'source-density-results.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps(summary,indent=2))

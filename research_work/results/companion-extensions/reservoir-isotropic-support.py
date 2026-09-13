"""Necessary monotonicity check for isotropic bound deposited particles."""
from pathlib import Path
import numpy as np
import json,hashlib
from scipy.special import logsumexp
from scipy.integrate import cumulative_trapezoid
P=Path(__file__).resolve().parent;L=P.parent/'isotropic-galaxy-transfer'
files=[P/'free-companion-results.json',P/'reservoir-scale-scan-results.json',L/'capacity-reference-optics-results.json',L/'third-radiation-retention-results.json']
free,scan,ref,fit=[json.loads(f.read_text()) for f in files];k=fit['models']['attenuated']['k0_per_kpc']
cases=[dict(label=v['Name']+' free',ac=v['capture_scale_kpc']) for v in free['rows']]
cases += [dict(label='J1621 scale '+str(v['capture_scale_Re']),ac=v['capture_scale_kpc']) for v in scan['rows'] if v['density_Msun_kpc3']>0]
cases += [dict(label=v['Name']+' transferred',ac=v['capture_scale_kpc']) for v in ref['rows'] if v['retention_mapping']['population']=='Chabrier']
def check(ac,n,order):
 x=np.geomspace(1e-4,1000,n);mu,w=np.polynomial.legendre.leggauss(order)
 t=x[:,None]*mu;b2=1+x[:,None]**2*(1-mu*mu);b=np.sqrt(b2)
 tau=k*ac*(t/(2*b2*(b2+t*t))+(np.arctan(t/b)+np.pi/2)/(2*b**3))
 logr=logsumexp(-tau+np.log(w/2),axis=1)-2*np.log1p(x*x)
 slope=np.gradient(logr,np.log(x),edge_order=2);density=np.exp(logr-max(logr));valid=(density>1e-12);valid[:2]=False;valid[-2:]=False
 rising=valid&(slope>1e-3);peak=int(np.argmax(logr));mass=cumulative_trapezoid(x*x*density,x,initial=0)
 return dict(violates=bool(np.any(rising)),maximum_reliable_log_slope=float(max(slope[valid])),peak_radius_kpc=float(x[peak]*ac),peak_x=float(x[peak]),central_sample_density_over_peak=float(density[0]),sampled_mass_fraction_inside_density_peak=float(mass[peak]/mass[-1]),rising_interval_kpc=[float(x[rising][0]*ac),float(x[rising][-1]*ac)] if any(rising) else None)
threshold=(-1.5*np.pi+np.sqrt((1.5*np.pi)**2+48))/2
out={'central_hollow_threshold_k0_ac':float(threshold),'central_hollow_threshold_ac_kpc':float(threshold/k),'scope':'Necessary test for spherical bound collisionless isotropic deposited particles; no full distribution-function construction','rows':[]}
for row in cases:
 coarse=check(row['ac'],4097,192);fine=check(row['ac'],8193,384);assert coarse['violates']==fine['violates']
 out['rows'].append(dict(**row,coarse=coarse,refined=fine));print(row['label'],fine,flush=True)
out['input_sha256']={f.as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in files+[Path(__file__),P/'reservoir-isotropic-support-protocol.md']}
(P/'reservoir-isotropic-support-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

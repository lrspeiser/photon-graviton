"""Local-well response diagnostic plus extraction of 57 bulge field summaries.

First run requires --html pointing to saved https://arxiv.org/html/2509.06846.
Later runs can use the saved bulge-inputs.json; optional HTML extraction needs bs4.
Field coordinates decoded from names are rounded angular labels, not stellar distances.
"""
from pathlib import Path
import argparse, hashlib, json, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy.polynomial.legendre import leggauss

H=Path(__file__).resolve().parent
save=lambda f,v:(H/f).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8',newline='\n')
p=argparse.ArgumentParser(); p.add_argument('--html',type=Path); a=p.parse_args()
if a.html:
    from bs4 import BeautifulSoup
    raw=a.html.read_bytes(); text=BeautifulSoup(raw,'html.parser').get_text(' ',strip=True)
    start=text.index('Fields characteristics'); stop=text.index('Notes.',start); text=text[start:stop]
    pattern=r'([pm]\d+(?:\.\d+)?[pm]\d+(?:\.\d+)?)\s+(\d+)\s+(-?\d+)\s+±\s+\\pm\s+(\d+)\s+(\d+\.?)\s+±\s+\\pm\s+(\d+)\s+(\d+|-)'
    fields=[]
    for m in re.finditer(pattern,text):
        name,N,mean,em,disp,ed,_=m.groups()
        coord=re.fullmatch(r'([pm])(\d+(?:\.\d+)?)([pm])(\d+(?:\.\d+)?)',name)
        signl,l,signb,b=coord.groups(); prefix=text[:m.start()]
        survey=max(['MUSE-inner','MUSE-outer','MUSE-V18','APOGEE','GIBS'],key=lambda s:prefix.rfind(s))
        fields.append(dict(survey=survey,field=name,N=int(N),l_label_deg=float(l)*(1 if signl=='p' else -1),b_label_deg=float(b)*(1 if signb=='p' else -1),
                           mean_heliocentric_RV_kms=float(mean),mean_error_kms=float(em),dispersion_kms=float(disp),dispersion_error_kms=float(ed)))
    assert len(fields)==57,len(fields)
    save('bulge-inputs.json',dict(source='https://arxiv.org/html/2509.06846',table=2,retrieved='2026-09-09',source_html_sha256=hashlib.sha256(raw).hexdigest(),
         status='Published field moments, not raw stellar phase-space data; approximate coordinates decoded from labels',fields=fields))
fields=json.loads((H/'bulge-inputs.json').read_text())['fields']

# Proposed local response, using the KNOWN softened/Plummer kernel.
# All configurations have the same total response amplitude 1000 (km/s)^2 kpc.
# This is not inferred companion energy, source mass or a Milky Way fit.
phi=np.arange(64)*2*np.pi/64
ring=np.c_[1.5*np.cos(phi),1.5*np.sin(phi),np.zeros(64)]
caps=np.vstack([np.c_[.5*np.cos(phi),.5*np.sin(phi),np.full(64,z)] for z in [-1.5,1.5]])
mu,w=leggauss(32)
shell=np.concatenate([np.c_[1.5*np.sqrt(1-u*u)*np.cos(phi),1.5*np.sqrt(1-u*u)*np.sin(phi),np.full(64,1.5*u)] for u in mu])
configs={'equatorial_deposits':(ring,np.ones(64)/64),'upper_lower_deposits':(caps,np.ones(128)/128),'whole_bulge_shell':(shell,np.repeat(w/2/64,64))}
amplitude=1000.; eps=.3
def response(R,z,centers,weights):
    delta=np.array([R,0,z])-centers; d2=np.sum(delta*delta,axis=1)+eps*eps
    potential=-amplitude*np.sum(weights/np.sqrt(d2))
    acc=-amplitude*np.sum(weights[:,None]*delta/d2[:,None]**1.5,axis=0)
    hzz=amplitude*np.sum(weights*(d2**(-1.5)-3*delta[:,2]**2*d2**(-2.5)))
    return float(potential),float(acc[0]),float(acc[2]),float(hzz)

rows=[]; checks={}; zz=np.linspace(-3,3,241)
for name,(centers,weights) in configs.items():
    assert abs(weights.sum()-1)<1e-12
    for R in [0,1,2,4,8]:
        for z in zz:
            potential,ar,az,hzz=response(R,z,centers,weights)
            rows.append(dict(configuration=name,R_kpc=R,z_kpc=float(z),delta_potential_kms2=potential,extra_aR_kms2_per_kpc=ar,extra_az_kms2_per_kpc=az,extra_vertical_curvature_kms2_per_kpc2=hzz))
    potential,ar,az,hzz=response(0,0,centers,weights)
    checks[name]=dict(center_depth=potential,center_acceleration_magnitude=float(np.hypot(ar,az)),center_vertical_curvature=hzz)
    assert np.hypot(ar,az)<1e-9
    # Independent finite difference of potential versus analytic acceleration.
    step=1e-5
    grad=(response(1,.7+step,centers,weights)[0]-response(1,.7-step,centers,weights)[0])/(2*step)
    assert abs(grad+response(1,.7,centers,weights)[2])<1e-5
save('bulge-response-predictions.json',rows)
save('bulge-results.json',dict(status='57 observed field summaries plus 3615 synthetic local-response positions; not a bulge dynamical fit',checks=checks,
    amplitude_kms2_kpc=amplitude,softening_kpc=eps,field_rows=len(fields),sum_field_star_counts=sum(x['N'] for x in fields),
    no_direct_force_inference='Line-of-sight dispersion is not vertical acceleration or a circular speed. Distances, proper motions, populations, selection and a bar model are required.'))
plt.rcParams.update({'svg.hashsalt':'bulge-local-depth-2026-09-09'})
fig,ax=plt.subplots(1,3,figsize=(16,4.8),layout='constrained')
labels={'equatorial_deposits':'Near the bulge equator','upper_lower_deposits':'Above and below the bulge','whole_bulge_shell':'Around the entire bulge'}
for name in configs:
    sub=[x for x in rows if x['configuration']==name and x['R_kpc']==1]
    ax[0].plot(zz,[x['delta_potential_kms2'] for x in sub],label=labels[name])
    ax[1].plot(zz,[x['extra_az_kms2_per_kpc'] for x in sub],label=labels[name])
ax[0].set(title='Same deposited response, different locations',xlabel='Height above/below disk (kpc)',ylabel='Added well depth ΔΦ, (km/s)²'); ax[0].legend(fontsize=8)
ax[1].axhline(0,c='gray',ls=':'); ax[1].set(title='Motion follows the resulting slope',xlabel='Height above/below disk (kpc)',ylabel='Extra vertical acceleration, (km/s)²/kpc')
for survey in ['MUSE-inner','MUSE-outer','MUSE-V18','APOGEE','GIBS']:
    sub=[x for x in fields if x['survey']==survey]
    ax[2].errorbar([abs(x['b_label_deg']) for x in sub],[x['dispersion_kms'] for x in sub],yerr=[x['dispersion_error_kms'] for x in sub],fmt='.',label=survey,alpha=.65)
ax[2].set(title='Observed bulge fields: velocities vary with location',xlabel='Absolute Galactic latitude (degrees)',ylabel='Line-of-sight velocity dispersion (km/s)'); ax[2].legend(fontsize=7)
fig.suptitle('Left/centre: synthetic response at R = 1 kpc, not fitted to the observations at right',fontsize=11)
fig.savefig(H/'bulge-summary.png',dpi=150)
print(json.dumps(checks,indent=2))

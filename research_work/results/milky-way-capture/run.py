"""Reproduce a finite Milky Way capture diagnostic, not a fitted galaxy theory.

Run: python research_work/results/milky-way-capture/run.py
Only numpy and matplotlib are required. Inputs retain published table provenance.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
G=4.30091727003628e-6  # kpc (km/s)^2 / solar mass
KPC_M=3.085677581491367e19
raw=HERE.joinpath('inputs.json').read_bytes()
data=json.loads(raw)

def angles(n,lo=-1,hi=1):
    x,w=leggauss(n)
    mu=lo+(hi-lo)*(x+1)/2
    w=w*(hi-lo)/2
    phi=(np.arange(2*n)+0.5)*np.pi/n
    nx=np.sqrt(1-mu[:,None]**2)*np.cos(phi)[None,:]
    nz=np.broadcast_to(mu[:,None],nx.shape)
    return nx,nz,w[:,None]/(4*n)

def intensity(R,z,beta,a=30.,h=3.,n=64,lo=-1.,hi=1.):
    """Integral I/I_boundary dOmega/(4pi), interior oblate ellipsoid.

    Backtrace x - s*n to the surface: A*s^2 - 2*B*s + C = 0.
    No focusing, scattering, internal sources, or varying absorption.
    """
    nx,nz,w=angles(n,lo,hi)
    aa=(1-nz*nz)/(a*a)+nz*nz/(h*h)
    bb=R*nx/(a*a)+z*nz/(h*h)
    cc=R*R/(a*a)+z*z/(h*h)-1
    assert cc<=1e-12
    s=(bb+np.sqrt(np.maximum(0,bb*bb-aa*cc)))/aa
    return float(np.sum(w*np.exp(-beta*s)))

rows=[]
for row in data['user_table']['rows']:
    r,p=row['R_kpc'],row['speed_excess_percent']
    observed=next(x for x in data['eilers']['rows'] if x['R_kpc']==r)
    v=observed['vc_kms']; vs=v/(1+p/100)
    extra=v*v-vs*vs
    rows.append(dict(R_kpc=r,quoted_excess_percent=p,vc_kms=v,
        reconstructed_stellar_speed_kms=vs,extra_speed_quadrature_kms=float(np.sqrt(extra)),
        extra_over_stellar_gravity=extra/(vs*vs),extra_fraction_total=extra/(v*v),
        extra_acceleration_m_s2=extra*1e6/(r*KPC_M),
        quoted_line_excess_percent=5.78+2.36*r))
r=np.array([x['R_kpc'] for x in rows]); p=np.array([x['quoted_excess_percent'] for x in rows])
coeff=np.polyfit(r,p,1); residual=p-np.polyval(coeff,r)
fit=dict(n=6,intercept_percent=float(coeff[1]),slope_percent_per_kpc=float(coeff[0]),
         R_squared=float(1-np.sum(residual**2)/np.sum((p-p.mean())**2)),
         RMS_percentage_points=float(np.sqrt(np.mean(residual**2))))

capture=[]; convergence=[]
for beta in [0.001,0.1,1.0]:
    for R in [0.,5.27,8.19,12.25,15.22,20.27,24.82,28.5]:
        for fraction in [0.,0.9]:
            z=fraction*3*np.sqrt(1-(R/30)**2)
            low=intensity(R,z,beta,n=64)
            high=intensity(R,z,beta,n=128)
            polar=intensity(R,z,beta,n=128,lo=.7,hi=1)+intensity(R,z,beta,n=128,lo=-1,hi=-.7)
            convergence.append(abs(low/high-1))
            capture.append(dict(beta_per_kpc=beta,R_kpc=R,z_kpc=float(z),
                fraction_local_halfheight=fraction,q_over_4pi_beta_I=high,
                polar_direction_fraction=polar/high))

checks={}
checks['max_64_vs_128_relative_difference']=max(convergence)
assert max(convergence)<1e-4
checks['transparent_limit']=intensity(8,1,0)
assert abs(checks['transparent_limit']-1)<1e-13
checks['sphere_center_error']=abs(intensity(0,0,.1,a=30,h=30)-np.exp(-3))
assert checks['sphere_center_error']<1e-13
checks['reflection_error']=abs(intensity(8,1,.1)-intensity(8,-1,.1))
assert checks['reflection_error']<1e-13
assert all(0<x['q_over_4pi_beta_I']<=1 for x in capture)
assert all(0<x['polar_direction_fraction']<1 for x in capture)
assert all(x['extra_acceleration_m_s2']>0 for x in rows)
assert len(data['eilers']['rows'])==38 and len(data['bovy']['rows'])==43
assert all(abs(x['R_kpc']+x['R0_minus_R_kpc']-8)<1e-12 for x in data['bovy']['rows'])

# Known flattened logarithmic potential: an identifiability demonstration.
# Phi=.5*v0^2*ln[(rc^2+R^2+z^2/q^2)/rc^2]. All q give identical
# radial forces at z=0. This is NOT a density produced by the capture toy.
R0=8.19; rc=1.; z=1.1
sun=rows[1]; vextra2=sun['extra_speed_quadrature_kms']**2
v02=vextra2*(R0*R0+rc*rc)/(R0*R0)
shape=[]
for q in [.75,1.,1.5]:
    for R in [5.27,8.19,12.25,15.22,20.27,24.82]:
        v2=v02*R*R/(R*R+rc*rc)
        kz=v02*z/(q*q*(rc*rc+R*R+z*z/(q*q)))
        shape.append(dict(q=q,R_kpc=R,z_kpc=z,extra_vc_kms=float(np.sqrt(v2)),
                          Kz_over_2piG_Msun_pc2=kz/(2*np.pi*G*1e6)))
checks['same_radial_different_vertical']=max(x['extra_vc_kms'] for x in shape if x['R_kpc']==R0)-min(x['extra_vc_kms'] for x in shape if x['R_kpc']==R0)
assert checks['same_radial_different_vertical']==0
assert abs(rows[1]['extra_over_stellar_gravity']-.5625)<1e-12
# Differentiate the potential numerically as a check on the force expressions.
force_errors=[]
for item in shape:
    R=item['R_kpc']; q=item['q']; z=item['z_kpc']; step=1e-4
    def phi(rr,zz):
        return .5*v02*np.log((rc*rc+rr*rr+zz*zz/(q*q))/(rc*rc))
    numerical_kz=(phi(R,z+step)-phi(R,z-step))/(2*step)
    numerical_v2=R*(phi(R+step,0)-phi(R-step,0))/(2*step)
    force_errors.extend([abs(numerical_kz/(item['Kz_over_2piG_Msun_pc2']*2*np.pi*G*1e6)-1),abs(numerical_v2/item['extra_vc_kms']**2-1)])
checks['potential_finite_difference_relative_error']=max(force_errors)
assert max(force_errors)<1e-7
# Galactic tidal order of magnitude, not an ephemeris residual prediction.
omega2=(sun['vc_kms']*1000/(R0*KPC_M))**2
au=149597870700.; solar_gm=1.32712440018e20
tidal=[dict(radius_AU=x,tidal_over_solar_acceleration=omega2*(x*au)**3/solar_gm) for x in [1,100]]
results=dict(input_sha256=hashlib.sha256(raw).hexdigest(),scope='Diagnostic and data ingestion; no joint Milky Way fit or source-budget validation',
             user_rows_recalculated=rows,six_point_fit=fit,all_angle_capture=capture,
             potential_shape_demo=shape,solar_tidal_order_of_magnitude=tidal,checks=checks)
HERE.joinpath('results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8',newline='\n')

plt.rcParams.update({'font.size':10,'svg.hashsalt':'milky-way-capture-2026-09-09'})
fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
ax[0,0].plot(r,[x['extra_acceleration_m_s2']/1e-11 for x in rows],'o-')
ax[0,0].set(xlabel='Galactic radius (kpc)',ylabel='Extra inward acceleration (10^-11 m/s²)',title='Your percentages + published circular speeds')
ax[0,0].text(.03,.06,'Stellar-only residual; gas still included in “extra”',transform=ax[0,0].transAxes,fontsize=9)
for beta in [.001,.1,1.]:
    subset=[x for x in capture if x['beta_per_kpc']==beta and x['fraction_local_halfheight']==0]
    central=subset[0]['q_over_4pi_beta_I']
    ax[0,1].plot([x['R_kpc'] for x in subset],[x['q_over_4pi_beta_I']/central for x in subset],'o-',label=f'capture length {1/beta:g} kpc')
ax[0,1].set(xlabel='Radius (kpc)',ylabel='Midplane deposition / central deposition',title='All-direction bath: illustrative 30 × 3 kpc ellipsoid')
ax[0,1].legend(fontsize=8)
br=data['bovy']['rows']
ax[1,0].errorbar([x['R_kpc'] for x in br],[x['Kz_over_2piG_Msun_pc2'] for x in br],yerr=[x['Kz_error'] for x in br],fmt='.',alpha=.7)
ax[1,0].set(xlabel='Radius (kpc), published R₀ = 8 kpc frame',ylabel='|Kz| / (2πG), M☉/pc²',title='43 published vertical-force estimates at |z| = 1.1 kpc')
ax[1,0].text(.03,.06,'Bovy & Rix 2013; total force, model-dependent inference',transform=ax[1,0].transAxes,fontsize=8)
for q in [.75,1,1.5]:
    subset=[x for x in shape if x['q']==q]
    ax[1,1].plot([x['R_kpc'] for x in subset],[x['Kz_over_2piG_Msun_pc2'] for x in subset],label=f'potential flattening q={q:g}')
ax[1,1].set(xlabel='Radius (kpc)',ylabel='Extra |Kz| / (2πG), M☉/pc²',title='Same midplane rotation, different vertical pull')
ax[1,1].legend(fontsize=8)
fig.savefig(HERE/'summary.png',dpi=160)
fig.savefig(HERE/'summary.svg',metadata={'Date':None})
svg=HERE/'summary.svg'; svg.write_text('\n'.join(x.rstrip() for x in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'fit':fit,'checks':checks,'user_rows':rows,'solar_shape':[x for x in shape if x['R_kpc']==R0],'tidal':tidal},indent=2))

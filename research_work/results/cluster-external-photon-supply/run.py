"""Exact external point-source photon/companion interception by a sphere."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre
from scipy.linalg import expm
OUT=Path(__file__).resolve().parent

def ray(y,d,a,k):
    entry=np.sqrt(d*d-1+y*y)-y
    length=2*y
    p0=np.exp(-a*entry); c0=-np.expm1(-a*entry)
    pg=np.exp(-a*length)
    cg=a*length*pg if abs(k-a)<1e-12 else a*(pg-np.exp(-k*length))/(k-a)
    p=p0*pg
    c=c0*np.exp(-k*length)+p0*cg
    dep=c0*(-np.expm1(-k*length))+p0*(1-pg-cg)
    return np.array([p,c,dep])

def source(d,a,k,n=128):
    if d<=1: raise ValueError('External source must have D/R > 1')
    hit=1/(2*d*d*(1+np.sqrt(1-1/d**2)))
    x,w=roots_legendre(n);y=(x+1)/2
    weight=w/2*y/(2*d*np.sqrt(d*d-1+y*y))
    energy=np.sum(ray(y,d,a,k)*weight,axis=1)
    return hit,energy

rows=[]; max_matrix=0
for d in [1.01,2,10,100,1000,10000]:
    for k in [.1,1,10]:
        a=.0002488993265191759
        hit,f=source(d,a,k); _,fine=source(d,a,k,256)
        reference=np.array([quad(lambda y:ray(y,d,a,k)[j]*y/(2*d*np.sqrt(d*d-1+y*y)),0,1,epsabs=1e-15,epsrel=1e-10)[0] for j in range(3)])
        err=float(np.max(np.abs(fine-reference))/hit)
        assert err<1e-9 and np.max(np.abs(f-fine))/hit<1e-9
        assert abs(np.sum(fine)/hit-1)<1e-10
        assert fine[2]>=0 and fine[2]<=hit
        for y in [.01,.5,1]:
            entry=np.sqrt(d*d-1+y*y)-y
            initial=np.array([np.exp(-a*entry),-np.expm1(-a*entry),0])
            matrix=np.array([[-a,0,0],[a,-k,0],[0,k,0]])
            e=float(np.max(np.abs(expm(matrix*2*y)@initial-ray(y,d,a,k))))
            max_matrix=max(max_matrix,e);assert e<1e-12
        rows.append(dict(distance_over_radius=d,alpha_R=a,kappa_R=k,intercepted_fraction=hit,
            missed_fraction=1-hit,escaping_photon_fraction=fine[0],escaping_companion_fraction=fine[1],
            deposited_fraction=fine[2],deposited_share_of_intercepted=fine[2]/hit,
            point_area_approximation_relative_error=abs((1/(4*d*d))/hit-1),
            quadrature_difference_per_intercepted=err))
# Degenerate limits: absent conversion or capture cannot deposit photon energy.
for a,k in [(0,1),(.1,0),(1,1)]:
    hit,f=source(2,a,k)
    assert abs(sum(f)/hit-1)<1e-10
    if a*k==0: assert abs(f[2])<1e-14
result=dict(scope='Isotropic external stellar pulse, no intervening capture, straight common-c propagation, fixed spherical receiver, constant conversion/capture. Pulse followed to completion; no astronomical catalog or chosen cosmic age.',matrix_max_error=max_matrix,rows=rows)
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
for r in rows:
    if r['kappa_R']==10: print(json.dumps(r))

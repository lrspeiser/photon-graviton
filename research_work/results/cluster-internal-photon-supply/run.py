"""Photon-funded internal capture in a stipulated sphere, not a cluster fit."""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
from scipy.special import roots_legendre
OUT=Path(__file__).resolve().parent

def fractions(a,k,s):
    photon=np.exp(-a*s)
    if abs(k-a)<1e-12:
        companion=a*s*photon
    else:
        # Integral of conversion at x followed by companion survival to s.
        companion=a*(np.exp(-a*s)-np.exp(-k*s))/(k-a)
    return np.array([photon,companion,1-photon-companion])

def avg(a,k):
    return np.array([quad(lambda l:.75*(1-l*l/4)*fractions(a,k,l)[j],0,2,epsabs=1e-12)[0] for j in range(3)])

def spatial(a,k,n):
    x,w=roots_legendre(n); r=(x+1)/2; wr=w/2
    mu,wm=roots_legendre(2*n)
    ell=-r[:,None]*mu+np.sqrt(1-r[:,None]**2+r[:,None]**2*mu**2)
    f=fractions(a,k,ell)
    return np.array([1.5*np.sum(wr*r*r*(v@wm)) for v in f])

rows=[]; max_matrix_error=0
for a in [.0002488993265191759,.01,1.0]:
    for k in [.1,1.,10.]:
        f=avg(a,k); f1=spatial(a,k,128); f2=spatial(a,k,256)
        matrix=np.array([[-a,0,0],[a,-k,0],[0,k,0]])
        for length in [0,.01,.5,2]:
            error=float(np.max(np.abs(expm(matrix*length)@np.array([1,0,0])-fractions(a,k,length))))
            max_matrix_error=max(max_matrix_error,error);assert error<1e-12
        assert abs(sum(f)-1)<1e-12 and min(f)>-1e-12
        assert np.max(np.abs(f2-f))<1e-7
        converted=1-f[0]
        rows.append(dict(alpha_R=a,kappa_R=k,escaped_photon_fraction=f[0],escaped_companion_fraction=f[1],
            deposited_fraction=f[2],maximum_deposit_with_perfect_capture=converted,
            captured_share_of_converted=f[2]/converted,spatial_path_max_difference=float(np.max(np.abs(f2-f))),
            refinement_difference=float(np.max(np.abs(f2-f1)))))
result=dict(scope='Uniform isotropic internal photon emission, constant alpha and kappa, same-direction c propagation, fixed sphere and permanent capture; full pulse followed until exit or deposit. No external supply, focusing, recycling, self-gravity or observed cluster fit.',matrix_crosscheck_max_error=max_matrix_error,rows=rows)
(OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))

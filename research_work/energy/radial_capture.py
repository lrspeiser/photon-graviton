"""Conditional steady radial transport from a point source, with no secondary loss.

Units: alpha and beta are inverse length; r uses the matching length unit.
Returned luminosities are fractions of the central photon luminosity.
Deposits are assumed to stay at their formation radii. No binding law is supplied.
"""
import numpy as np
from scipy.integrate import quad


def fractions(radius, alpha, beta):
    r=np.asarray(radius,dtype=float)
    if not np.isfinite(r).all() or np.any(r<0):
        raise ValueError('Radius must be finite and nonnegative.')
    if any(not np.isfinite(v) or v<0 for v in [alpha,beta]):
        raise ValueError('Rates must be finite and nonnegative.')
    p=np.exp(-alpha*r)
    if alpha==0:
        c=np.zeros_like(r)
    elif alpha==beta:
        c=alpha*r*p
    elif beta>alpha:
        c=alpha*np.exp(-alpha*r)*(-np.expm1(-(beta-alpha)*r))/(beta-alpha)
    else:
        c=alpha*np.exp(-beta*r)*(-np.expm1(-(alpha-beta)*r))/(alpha-beta)
    d=-np.expm1(-alpha*r)-c
    # For very short paths, subtracting two first-order terms loses relative
    # accuracy in the second-order deposited fraction. Integrate the positive
    # conversion-and-subsequent-capture kernel in those cases.
    flat=np.asarray(d).reshape(-1).copy()
    rr=r.reshape(-1)
    for i,x in enumerate(rr):
        if max(alpha,beta)*x<.05:
            flat[i]=quad(lambda s: alpha*np.exp(-alpha*s)*(-np.expm1(-beta*(x-s))),0,float(x),epsabs=1e-25,epsrel=1e-11)[0]
    return p,c,flat.reshape(r.shape)


def deposited_velocity_km_s(radius_kpc, luminosity_Lsun, duration_years, alpha_per_kpc, beta_per_kpc):
    """Deposited contribution only, under stipulated Newtonian cold-mass response."""
    r=np.asarray(radius_kpc,dtype=float)
    if np.any(r<=0):
        raise ValueError('Velocity requires strictly positive radii.')
    if any(not np.isfinite(v) or v<0 for v in [luminosity_Lsun,duration_years]):
        raise ValueError('Luminosity and duration must be finite and nonnegative.')
    d=fractions(r,alpha_per_kpc,beta_per_kpc)[2]
    mass=luminosity_Lsun*3.828e26*duration_years*31557600*d/(299792458.**2*1.98847e30)
    return np.sqrt(4.30091727003628e-6*mass/r)

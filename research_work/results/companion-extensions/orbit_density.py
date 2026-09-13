"""Shared Abel deprojection and tracer slope; no fitting or file I/O."""
import numpy as np
from numpy.polynomial.legendre import leggauss

def deproject_slope(profile,x,order):
    nodes,weights=leggauss(order);total=np.zeros_like(x);slope_weight=np.zeros_like(x)
    for comp in profile['components']:
        Re=comp['R_arcsec']/profile['computed_equal_area_half_light_arcsec'];n=comp['n'];bn=comp['bn'];amp=comp['amp_at_R']
        upper=np.arccosh(np.maximum(1.,Re*(1+100/bn)**n/x))
        u=upper[:,None]*(nodes+1)/2;R=x[:,None]*np.cosh(u)
        q=(R/Re)**(1/n)
        integrand=amp*np.exp(-bn*(q-1))*bn/(n*Re)*(R/Re)**(1/n-1)
        density=np.sum(integrand*weights,axis=1)*upper/(2*np.pi)
        numerator=np.sum(integrand*(1-1/n+bn*q/n)*weights,axis=1)*upper/(2*np.pi)
        total+=density;slope_weight+=numerator
    assert np.all(total>0)
    return total,slope_weight/total

"""Gas-only McMillan profile via galpy DiskSCF; no dark halo imported.

Source: McMillan 2017 Table 1, https://arxiv.org/html/1608.00971
Dependency: galpy (tested version recorded in generated metadata).
"""
import numpy as np
from galpy.potential import DiskSCFPotential

RO=8.; VO=220.; G=4.30091727003628e-6
PARAMS=[dict(Sigma0=53.1,Rd=7.,Rm=4.,h=.085),dict(Sigma0=2180.,Rd=1.5,Rm=12.,h=.045)]

def make_gas(order=30):
    # Natural mass unit is VO^2*RO/G; surface unit is VO^2/(G*RO).
    sigmaunit=VO**2/(G*RO)/1e6
    def density(R,z):
        R=np.asarray(R); safe=np.maximum(R,1e-12)
        return sum((p['Sigma0']/sigmaunit)/(4*p['h']/RO)*np.exp(-p['Rm']/RO/safe-safe/(p['Rd']/RO)-2*np.logaddexp(z/(2*p['h']/RO),-z/(2*p['h']/RO))+2*np.log(2)) for p in PARAMS)
    return DiskSCFPotential(dens=density,
        Sigma=[dict(type='exp',h=p['Rd']/RO,amp=p['Sigma0']/sigmaunit,Rhole=p['Rm']/RO) for p in PARAMS],
        hz=[dict(type='sech2',h=p['h']/RO) for p in PARAMS],a=2.5,N=order,L=order,ro=RO,vo=VO)

def force(pot,R,z):
    # Positive inward radial / downward vertical force in (km/s)^2/kpc.
    return -pot.Rforce(R/RO,z/RO,use_physical=False)*VO**2/RO,-pot.zforce(R/RO,z/RO,use_physical=False)*VO**2/RO

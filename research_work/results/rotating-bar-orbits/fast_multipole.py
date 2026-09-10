"""Evaluate the same harmonic expansion using shared Legendre recurrences."""
import numpy as np
from scipy.special import sph_legendre_p_all
from field import Multipole


class FastMultipole(Multipole):
    def evaluate(self,xyz):
        xyz=np.atleast_2d(xyz).astype(float);x,y,z=xyz.T
        R=np.hypot(x,y);r=np.sqrt(R*R+z*z)
        if np.any((r<self.r[0])|(r>self.r[-1])|(R==0)):
            raise ValueError('Outside numerical domain or exactly on polar axis')
        theta=np.arctan2(R,z);phi=np.arctan2(y,x)
        # Known Y_lm=P_lm(theta)*exp(i*m*phi). Compute all P once instead
        # of recomputing the lower-order recurrences for every requested Y.
        leg=sph_legendre_p_all(int(self.l.max()),int(self.m.max()),theta,diff_n=1)
        norm=np.where(self.m[:,None]>0,np.sqrt(2),1)
        co=np.cos(self.m[:,None]*phi);si=np.sin(self.m[:,None]*phi)
        Y=leg[0,self.l,self.m]*co*norm
        dtheta=leg[1,self.l,self.m]*co*norm
        dphi=-self.m[:,None]*leg[0,self.l,self.m]*si*norm
        coef=self.spline(np.log(r));dr=self.spline(np.log(r),1)/r[:,None]
        potential=np.sum(coef*Y.T,axis=1)
        gr=np.sum(dr*Y.T,axis=1)
        gt=np.sum(coef*dtheta.T,axis=1)
        gp=np.sum(coef*dphi.T,axis=1)
        aR=-gr*R/r-gt*z/r**2;az=-gr*z/r+gt*R/r**2;ap=-gp/R
        return potential,np.c_[aR*x/R-ap*y/R,aR*y/R+ap*x/R,az]

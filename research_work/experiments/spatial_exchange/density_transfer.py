"""Scalar-gradient extension with exact reciprocal Hamiltonian reaction."""
from joint_transfer import JointTransfer
from local_rotor import derivative


class DensityTransfer(JointTransfer):
    def __init__(self,zeta=.5,**kwargs):
        super().__init__(**kwargs);self.zeta=zeta

    def evaluate(self,y,flow=True):
        mapped=y.copy();f,m,_,_=self.unpack(mapped)
        for j in range(3):m[4+j]-=self.zeta*derivative(f[0],j,self.h)
        result=super().evaluate(mapped,flow)
        if flow:
            fd,pd,_,_=self.unpack(result)
            pd[0]-=self.zeta*sum(derivative(fd[4+j],j,self.h) for j in range(3))
        return result

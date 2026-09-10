"""Independent analytic limits and potential-gradient checks, without fitting."""
import json
import numpy as np
from scipy.integrate import quad
from run import HERE, G, A, P, ASTAR, Completion, force

class Analytic:
    def __init__(self,kuzmin=False):self.kuzmin=kuzmin
    def evaluate(self,r,mu):
        R=r*np.sqrt(1-mu**2);z=r*mu
        mass=5e10;scale=3.
        q=R*R+(scale+abs(z))**2 if self.kuzmin else r*r+scale*scale
        phi=-G*mass/np.sqrt(q)
        gR=G*mass*R/q**1.5
        gz=G*mass*(scale+abs(z))*np.sign(z)/q**1.5 if self.kuzmin else G*mass*z/q**1.5
        return phi,gR*R/r+gz*z/r,gR*z/r-gz*R/r

R=np.array([1.,3.,5.,8.,12.,20.])
z=np.array([.5,1.1,2.,1.1,2.,3.])
checks={}
for name,source in [('spherical_plummer',Analytic()),('kuzmin_disk',Analytic(True))]:
    # A razor-thin disk has a force jump, unlike the finite-thickness MW disks.
    # Use explicitly higher angular order for this singular analytic limit.
    out=Completion(source,True,400,ellmax=512,angular_nodes=2048)
    ordinary=force(source,R,z)
    exact=ordinary*(A*(np.linalg.norm(ordinary,axis=1)/ASTAR)**(P-1))[:,None]
    calculated=force(out,R,z)
    relative=np.linalg.norm(calculated-exact,axis=1)/np.linalg.norm(exact,axis=1)
    # In each half-space Kuzmin's potential depends on distance to one shifted
    # point source. The algebraic extra force is therefore already conservative.
    gradient_errors=[]
    for rr,zz in zip(R,z):
        def potential(rr,zz):
            rad=np.hypot(rr,zz)
            return out.evaluate(np.array([rad]),np.array([zz/rad]))[0][0]
        eps=1e-4
        numerical=-np.array([(potential(rr+eps,zz)-potential(rr-eps,zz))/(2*eps),
                             (potential(rr,zz+eps)-potential(rr,zz-eps))/(2*eps)])
        actual=force(out,rr,zz)[0]
        gradient_errors.append(float(np.linalg.norm(numerical-actual)/np.linalg.norm(actual)))
    def work(order):
        x,w=np.polynomial.legendre.leggauss(order);total=0.
        for start,end in [([2,.4],[8,.4]),([8,.4],[8,2]),([8,2],[2,2]),([2,2],[2,.4])]:
            start=np.array(start);end=np.array(end)
            points=(start+end)/2+x[:,None]*(end-start)/2
            total+=np.sum(w*(force(out,points[:,0],points[:,1])@((end-start)/2)))
        return float(total)
    loop=work(256);loop_low=work(128)
    checks[name]=dict(max_analytic_extra_force_relative_error=float(relative.max()),
        relative_errors=relative.tolist(),max_potential_gradient_relative_error=max(gradient_errors),
        closed_loop_work_kms2=loop,closed_loop_quadrature_change_kms2=abs(loop-loop_low))
    assert relative.max()<.01,checks[name]
    assert max(gradient_errors)<1e-6,checks[name]
    assert abs(loop)<.1,checks[name]

# Disk mass independently from the radial surface density, including holes.
mass=0.
for sigma,rd,hole in [(1.332e9,2.,2.7),(8.97e8,2.8,2.7),(5.81e7,7.,4.),(2.68e9,1.5,12.)]:
    mass+=2*np.pi*sigma*quad(lambda R:R*np.exp(-hole/R-R/rd),0,np.inf,epsabs=1e-9)[0]
checks['disk_mass_independent_Msun']=mass
(HERE/'analytic-verification.json').write_text(json.dumps(checks,indent=2)+'\n',newline='\n')
print(json.dumps(checks,indent=2))

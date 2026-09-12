from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent

class Model:
    def __init__(self,R,n,nq):
        self.R=R;self.n=n
        self.edges=e=np.r_[0.,np.geomspace(1e-3,R,n)]
        self.a=a=e[:-1];self.b=b=e[1:];self.r=r=(a+b)/2
        self.V=4*np.pi/3*(b**3-a**3)
        z,w=roots_legendre(nq)
        u=(a*a)[:,None]+(b*b-a*a)[:,None]*(1+z)/2
        self.area=(np.pi*(b*b-a*a)[:,None]*w/2).ravel()
        self.length=np.sqrt(np.maximum(b[None,:]**2-u.ravel()[:,None],0))-np.sqrt(np.maximum(a[None,:]**2-u.ravel()[:,None],0))
        self.wb=1/np.sqrt(1+r*r)
    def gravity(self,rho):
        mass=rho*self.V;before=np.cumsum(mass)-mass
        shell=2*np.pi*rho*(self.b**2-self.a**2)
        outside=np.cumsum(shell[::-1])[::-1]-shell
        depth=(before+4*np.pi/3*rho*(self.r**3-self.a**3))/self.r+2*np.pi*rho*(self.b**2-self.r**2)+outside
        return mass,depth,float(shell.sum())
    def absorption(self,kappa):
        t=self.length*kappa[None,:]
        lower=np.cumsum(t,axis=1)-t
        half=t.sum(axis=1)
        higher=np.maximum(half[:,None]-lower-t,0)
        absorbed=(-np.expm1(-t))*(np.exp(-higher)+np.exp(-half[:,None]-lower))
        shell=self.area@absorbed
        total=float(self.area@(-np.expm1(-2*half)))
        transmitted=float(self.area@np.exp(-2*half))
        return shell,total,transmitted
    def run(self,p,C,k0,fine):
        B=C/k0
        def rhs(s,y):
            _,wd,_=self.gravity(y[:self.n])
            W=self.wb+wd;k=k0*(W/(1+W))**p
            shell,total,out=self.absorption(k)
            return np.r_[B*shell/self.V,B*total,B*out]
        sol=solve_ivp(rhs,(0,1),np.zeros(self.n+2),rtol=1e-8 if fine else 2e-6,atol=1e-12 if fine else 1e-10)
        assert sol.success and sol.y.min()>=-1e-10
        rho=sol.y[:self.n,-1];mass,wd,central=self.gravity(rho)
        cumulative=np.r_[0.,np.cumsum(mass)]
        v=[float(np.interp(x**3,self.edges**3,cumulative)/x) for x in (1.,3.,10.)]
        total=float(mass.sum());absorbed=float(sol.y[-2,-1]);out=float(sol.y[-1,-1]);incoming=B*np.pi*self.R**2
        ledger=abs(total+out-incoming)/incoming
        chord=abs(total-absorbed)/max(total,1e-100)
        assert ledger<1e-7 and chord<1e-7
        return dict(p=p,C=C,kappa0=k0,B=B,R=self.R,shells=self.n,evaluations=sol.nfev,
                    mass=total,central_depth=central,vc2=v,incoming=incoming,transmitted=out,
                    capture_fraction=total/incoming,energy_ledger_error=ledger,chord_mass_error=chord,
                    edge_density=float(rho[-1]),diagnostics=[total,central]+v)

def main():
    results=[];checks=[];analytic=[]
    for R in (30.,100.):
        coarse=Model(R,128,4);fine=Model(R,256,8)
        for p,C in ((4,1.),(6,1.),(6,100.)):
            for k in (.1,10.,1000.):
                c=coarse.run(p,C,k,False);f=fine.run(p,C,k,True)
                change=float(max(abs(np.array(c['diagnostics'])/f['diagnostics']-1)))
                checks.append(dict(R=R,p=p,C=C,kappa0=k,max_relative_change=change,passes=change<.01))
                results.append(dict(coarse=c,refined=f))
                print(R,p,C,k,f['mass'],f['capture_fraction'],change,flush=True)
        for k in (0.,.001,1.):
            shells,total,out=fine.absorption(np.full(fine.n,k))
            exact=quad(lambda u:np.pi*(-np.expm1(-2*k*np.sqrt(max(R*R-u,0)))),0,R*R,epsabs=1e-9)[0]
            error=abs(total-exact)/max(exact,1.)
            analytic.append(dict(R=R,kappa=k,absorption_relative_error=error,passes=error<.001))
        density=3/(4*np.pi)*(1+fine.r*fine.r)**-2.5
        _,wd,_=fine.gravity(density)
        exact=1/np.sqrt(1+fine.r*fine.r)-(1+R*R)**-1.5
        error=float(max(abs(wd/exact-1)))
        analytic.append(dict(R=R,gravity_relative_error=error,passes=error<.005))
    out=dict(results=results,refinement_checks=checks,analytic_checks=analytic)
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8',newline='\n')
    assert all(c['passes'] for c in checks+analytic),'Failed gate retained; refine affected cases'

if __name__ == "__main__":
    main()

"""Axisymmetric QUMOND-style completion of the frozen empirical extra force.

Known Poisson/Legendre mathematics; no fitted parameters or photon derivation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy.special import eval_legendre, roots_legendre
from scipy.interpolate import CubicSpline
from scipy.integrate import cumulative_trapezoid

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
G=4.30091727003628e-6
FIT=HERE.parent/'joint-galaxy-audit/results.json'
PARAM=json.loads(FIT.read_text())['sparc']['parameters']
A,P=PARAM['A'],PARAM['p']
ASTAR=PARAM['a_star_m_s2']*3.085677581491367e19/1e6

def basis(ells,mu):
    v=eval_legendre(ells[:,None],mu[None,:]).T
    prev=eval_legendre(np.maximum(ells-1,0)[:,None],mu[None,:]).T
    d=ells[None,:]*(prev-mu[:,None]*v)/(1-mu[:,None]**2)
    return v,d

def scaled_integrals(r,inner,outer,ells):
    """Exact power-kernel integration of a log-linearly interpolated source.

    Treating large-l kernels with endpoint trapezoids produces appreciable
    error even for a smooth source; integrate their exponential factors.
    """
    a=np.zeros_like(inner);b=np.zeros_like(outer)
    for i in range(1,len(r)):
        h=np.log(r[i]/r[i-1]);power=ells+2;x=power*h
        q=np.exp(-h*(ells+1))
        k0=-np.expm1(-x)/power
        k1=(np.expm1(-x)+x*np.exp(-x))/power**2
        a[i]=q*a[i-1]+r[i]*(inner[i]*k0+(inner[i]-inner[i-1])*k1/h)
    for i in range(len(r)-2,-1,-1):
        h=np.log(r[i+1]/r[i]);power=1-ells;x=power*h
        q=np.exp(-h*ells)
        k0=np.expm1(x)/power
        k1=(np.expm1(x)*(x-1)+x)/power**2
        b[i]=q*b[i+1]+r[i]*(outer[i]*k0+(outer[i+1]-outer[i])*k1/h)
    return a,b

class Expansion:
    def __init__(self,r,ells,phi,grad):
        self.r=r;self.ells=ells
        self.phi=CubicSpline(np.log(r),phi,axis=0)
        self.grad=CubicSpline(np.log(r),grad,axis=0)
    def evaluate(self,r,mu):
        r,mu=np.broadcast_arrays(np.atleast_1d(r),np.atleast_1d(mu))
        b,d=basis(self.ells,mu)
        # Differentiate the very same interpolated potential used for angular
        # forces. The independent integral-derived radial field remains stored
        # in self.grad for checking, but is not a separate force prescription.
        ph=self.phi(np.log(r));gr=self.phi(np.log(r),1)/r[:,None]
        return np.sum(ph*b,axis=1),np.sum(gr*b,axis=1),-np.sqrt(1-mu**2)*np.sum(ph*d,axis=1)/r

class Disk(Expansion):
    def __init__(self,refined):
        r=np.geomspace(1e-5,500,3072 if refined else 1536)
        ells=np.arange(0,257 if refined else 129,2)
        mu,w=roots_legendre(1024 if refined else 512)
        leg,_=basis(ells,mu)
        angular=leg*(w[:,None]*(2*ells+1)[None,:]/2)
        coeff=np.empty((len(r),len(ells)))
        params=[(1.332e9,2.,.3,2.7,'exp'),(8.97e8,2.8,.9,2.7,'exp'),
                (5.81e7,7.,.085,4.,'sech2'),(2.68e9,1.5,.045,12.,'sech2')]
        for start in range(0,len(r),32):
            rr=r[start:start+32,None];R=rr*np.sqrt(1-mu**2);z=rr*mu
            rho=np.zeros_like(R)
            for sigma,rd,h,hole,kind in params:
                surface=sigma*np.exp(-hole/R-R/rd)
                if kind=='exp':rho+=surface*np.exp(-abs(z)/h)/(2*h)
                else:rho+=surface*np.exp(-2*np.logaddexp(z/(2*h),-z/(2*h))+2*np.log(2))/(4*h)
            coeff[start:start+len(rr)]=rho@angular
        inn,out=scaled_integrals(r,coeff*r[:,None],coeff*r[:,None],ells)
        fac=4*np.pi*G/(2*ells+1)
        super().__init__(r,ells,-fac*(inn+out),fac*((ells+1)*inn-ells*out)/r[:,None])
        self.mass=float(4*np.pi*inn[-1,0]*r[-1])

class CachedAxisymmetric:
    def __init__(self,path):
        f=np.load(path);mask=f['m']==0
        self.ells=f['l'][mask]
        coeff=f['coeff'][:,mask]*np.sqrt((2*self.ells+1)/(4*np.pi))
        self.spline=CubicSpline(np.log(f['r']),coeff,axis=0)
        self.edge=100.
    def evaluate(self,r,mu):
        b,d=basis(self.ells,mu);s=np.minimum(r,self.edge)
        ph=self.spline(np.log(s));gr=self.spline(np.log(s),1)/s[:,None]
        outside=r>self.edge
        if outside.any():
            ph[outside]*=(self.edge/r[outside,None])**(self.ells+1)
            gr[outside]=-(self.ells+1)*ph[outside]/r[outside,None]
        return np.sum(ph*b,axis=1),np.sum(gr*b,axis=1),-np.sqrt(1-mu**2)*np.sum(ph*d,axis=1)/r

class Baryons:
    def __init__(self,refined):
        self.disk=Disk(refined)
        cache=ROOT/'research_work/data-cache/bar-field'
        self.bar=CachedAxisymmetric(cache/'bar-L64.npz')
        self.nuclei=CachedAxisymmetric(cache/'nuclei-L16.npz')
    def evaluate(self,r,mu):
        r,mu=np.broadcast_arrays(np.atleast_1d(r),np.atleast_1d(mu))
        out=np.array(self.disk.evaluate(r,mu))
        for component in [self.bar,self.nuclei]:out+=component.evaluate(r,mu)
        q=r*r+.001**2
        out[0]-=G*4.1e6/np.sqrt(q);out[1]+=G*4.1e6*r/q**1.5
        return out

class Completion(Expansion):
    def __init__(self,baryons,refined,rmax=200,ellmax=None,angular_nodes=None):
        r=np.geomspace(1e-4,rmax,1200 if refined else 600)
        ells=np.arange(0,(128 if refined else 64) + 1 if ellmax is None else ellmax+1,2)
        mu,w=roots_legendre((512 if refined else 256) if angular_nodes is None else angular_nodes)
        leg,dleg=basis(ells,mu)
        angular=leg*w[:,None]*(2*ells+1)[None,:]/2
        tangent=dleg*(w*np.sqrt(1-mu**2))[:,None]*(2*ells+1)[None,:]/2
        br=np.empty((len(r),len(ells)));tt=np.empty_like(br)
        for start in range(0,len(r),8):
            rr=np.broadcast_to(r[start:start+8,None],(len(r[start:start+8]),len(mu)))
            mm=np.broadcast_to(mu,rr.shape)
            _,gr,gt=baryons.evaluate(rr.ravel(),mm.ravel())
            gr=gr.reshape(rr.shape);gt=gt.reshape(rr.shape)
            strength=np.hypot(gr,gt)
            mult=A*(strength/ASTAR)**(P-1)
            br[start:start+len(rr)]=(mult*gr)@angular
            tt[start:start+len(rr)]=(mult*gt)@tangent
        # Integration by parts projects divergence without noisy differentiation.
        inn,out=scaled_integrals(r,tt-ells*br,tt+(ells+1)*br,ells)
        phi=-(inn+out)/(2*ells+1)
        grad=br+((ells+1)*inn-ells*out)/((2*ells+1)*r[:,None])
        # Monopole has a divergent absolute potential for this uncut power law.
        # Fix a gauge at inner radius; its force is exactly the radial average.
        phi[:,0]=cumulative_trapezoid(br[:,0],r,initial=0)
        grad[:,0]=br[:,0]
        super().__init__(r,ells,phi,grad)
        self.radial_average=br[:,0]

def force(component,R,z):
    R,z=np.broadcast_arrays(np.atleast_1d(R),np.atleast_1d(z))
    r=np.hypot(R,z);_,gr,gt=component.evaluate(r,z/r)
    return np.c_[-gr*R/r-gt*z/r,-gr*z/r+gt*R/r]

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--refined',action='store_true')
    parser.add_argument('--outer',type=float,default=200.)
    args=parser.parse_args();suffix=('-refined' if args.refined else '')+(f'-outer{args.outer:g}' if args.outer!=200 else '')
    print('Building global ordinary-matter field',flush=True)
    baryons=Baryons(args.refined)
    print('Projecting frozen empirical response',flush=True)
    completion=Completion(baryons,args.refined,args.outer)
    data=json.loads((HERE.parent/'milky-way-capture/inputs.json').read_text())
    rows=[];scores={}
    for observable,source in [('rotation',data['eilers']['rows']),('vertical',data['bovy']['rows'])]:
        R=np.array([row['R_kpc'] for row in source]);z=np.zeros(len(R)) if observable=='rotation' else np.full(len(R),1.1)
        ordinary=force(baryons,R,z);additional=force(completion,R,z)
        observed=np.array([row['vc_kms'] if observable=='rotation' else row['Kz_over_2piG_Msun_pc2'] for row in source])
        for label,a in [('ordinary',ordinary),('conservative_completion',ordinary+additional)]:
            predicted=np.sqrt(-R*a[:,0]) if observable=='rotation' else abs(a[:,1])/(2*np.pi*G*1e6)
            scores[f'{observable}/{label}']=dict(rows=len(R),rmse=float(np.sqrt(np.mean((predicted-observed)**2))),bias=float(np.mean(predicted-observed)))
            for i in range(len(R)):
                rows.append(dict(observable=observable,model=label,R_kpc=float(R[i]),z_kpc=float(z[i]),observed=float(observed[i]),predicted=float(predicted[i])))
    result=dict(classification='Known quasi-linear conservative completion of frozen empirical fit; not a photon mechanism or fresh holdout',
        parameters=PARAM,refined=args.refined,outer_radius_kpc=args.outer,disk_mass_Msun=baryons.disk.mass,scores=scores,
        fit_sha256=hashlib.sha256(FIT.read_bytes()).hexdigest(),holdout_scores_opened=False,new_parameters_fitted=False)
    for name,obj in [('results',result),('predictions',rows)]:
        (HERE/f'{name}{suffix}.json').write_text(json.dumps(obj,indent=2)+'\n',newline='\n')
    print(json.dumps(result,indent=2),flush=True)

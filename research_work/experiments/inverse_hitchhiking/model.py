"""Deterministic phenomenological response and static optical projection."""
import numpy as np
from scipy.special import expit
from scipy.optimize import least_squares
from data import CODE_TO_SI,G,C
FAMILIES={"L":("c","q"),"R":("c","q","s"),"M":("c","q","m"),
          "RM":("c","q","s","m"),"T":("c","q","s","m","d")}
CAPS=(3e-9,1e-8,3e-8,1e-7)
def features(gb,r,mass,family):
    gb,r=np.broadcast_arrays(np.asarray(gb,float),np.asarray(r,float))
    mass=np.broadcast_to(mass,gb.shape)
    columns=dict(c=np.ones_like(gb),q=np.log(np.maximum(gb*CODE_TO_SI,1e-30)/1e-10),
        s=np.log(r/10),m=np.log(mass/1e11),d=np.tanh(np.log(mass/1e12)/2))
    return np.stack([columns[k] for k in FAMILIES[family]],axis=-1)
def extra(theta,capacity,F):
    return capacity/CODE_TO_SI*expit(F@theta)
def bounds(family):
    limits=dict(c=(-25,5),q=(-2,3),s=(-2,2),m=(-2,2),d=(-10,10))
    return np.array([limits[k][0] for k in FAMILIES[family]]),np.array([limits[k][1] for k in FAMILIES[family]])
def alpha(g,b,order=128,breaks=()):
    b=float(b);edges=[0.]+sorted(np.arccos(b/r) for r in breaks if r>b)+[np.pi/2]
    nodes,weights=np.polynomial.legendre.leggauss(order)
    total=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        angle=lo+(nodes+1)*(hi-lo)/2
        r=b/np.cos(angle)
        total+=float(np.sum(weights*g(r)*r)*(hi-lo)/2)
    return 4/C**2*total
def shear_shape(g,r,order=128,h=.002,breaks=()):
    a=np.array([alpha(g,b,order,breaks) for b in r])
    derivative=np.array([(alpha(g,b*(1+h),order,breaks)-alpha(g,b*(1-h),order,breaks))/(2*h*b) for b in r])
    return dict(alpha=a,shape=a/r-derivative)
def optical_force(source,fit,reach,eta=1):
    def law(r):
        r=np.asarray(r)
        mb=np.interp(r,source["r"],source["mb"],left=source["mb"][0],right=source["mass"])
        gb=G*mb/r**2
        rc=np.minimum(r,reach)
        mc=np.interp(rc,source["r"],source["mb"],left=source["mb"][0],right=source["mass"])
        base=G*mc/rc**2
        ah=extra(np.asarray(fit["theta"]),fit["capacity"],features(base,rc,source["mass"],fit["family"]))
        ah*=np.minimum(1.,(reach/r)**2)
        return gb+eta*ah
    return law
def shape_fit(shape,y,error):
    amp=max(0.,float(np.sum(shape*y/error**2)/np.sum((shape/error)**2)))
    pred=amp*shape
    return dict(amplitude=amp,prediction=pred,chi2=float(np.sum(((pred-y)/error)**2)))
def controls():
    results=[]
    def add(name,value,tol):
        results.append(dict(name=name,value=float(value),limit=tol,passed=bool(value<=tol)))
    f0=.2;feq=.8;ell=3.;steps=100;ds=.01
    f=f0
    for _ in range(steps):f=feq+(f-feq)*np.exp(-ds/ell)
    add("deterministic relaxation exact",abs(f-(feq+(f0-feq)*np.exp(-steps*ds/ell))),1e-12)
    z=np.linspace(-30,30,101);p=expit(z)
    add("positive bounded fraction",float(max(0,-p.min(),p.max()-1)),0)
    x=np.geomspace(.01,100,100);p=x/(1+x)
    add("inverse required odds",float(np.max(abs(p/(1-p)/x-1))),1e-12)
    add("one code acceleration in SI",abs(CODE_TO_SI*3.085677581491367e19/1e6-1),1e-14)
    mass=1e12;rr=np.array([1.,10.,100.])
    pred=shear_shape(lambda r:G*mass/r**2,rr)
    add("point mass deflection",np.max(abs(pred["alpha"]/(4*G*mass/(rr*C*C))-1)),1e-10)
    add("point mass shear finite difference",np.max(abs(pred["shape"]/(8*G*mass/(rr*rr*C*C))-1)),3e-6)
    rr=np.geomspace(.1,100,40);gb=np.geomspace(1e-12,1e-8,40)/CODE_TO_SI
    F=features(gb,rr,1e11,"L");truth=np.array([-4.,.7]);target=extra(truth,1e-8,F)
    result=least_squares(lambda t:(extra(t,1e-8,F)-target)/max(target),[-3.,.5],
                         xtol=1e-13,ftol=1e-13,gtol=1e-13)
    add("synthetic parameter recovery",np.max(abs(result.x-truth)),1e-7)
    if not all(r["passed"] for r in results):raise AssertionError(results)
    return results

"""CMF-2 source-gradient features; CMF-1 source stays immutable."""
import laws as B
from laws import *
BASES={**B.BASES,"G":B.BASES["RM"]+["s","xs","ys","zs","s2"],
       "GH":B.BASES["RMH"]+["s","xs","ys","zs","s2","hs"]}
SPECS=[("G",0.)]+[("GH",e) for e in (.3,1.,3.)]

def design(gb,r,mass,family,ell):
    if family not in ("G","GH"):return B.design(gb,r,mass,family,ell)
    base=B.design(gb,r,mass,"RM" if family=="G" else "RMH",ell)
    logged=np.log(np.maximum(np.asarray(gb)*D.CODE_TO_SI,1e-30))
    derivative=np.gradient(logged,np.log(r)) if len(logged)>1 else np.zeros_like(logged)
    s=np.tanh(derivative/2)
    extra=[s,base[:,1]*s,base[:,4]*s,base[:,7]*s,s*s]
    if family=="GH":extra.append(base[:,11]*s)
    return np.column_stack([base,*extra])

def prepared(galaxies,clusters,spec):
    family,ell=spec
    return [{**g,"F":design(g["gb"],g["r"],g["mass"],family,ell)} for g in galaxies],[
        {**c,"F":design(c["gb"],c["r"],c["mass"],family,ell)} for c in clusters]

def candidate_force(fit):
    return lambda o:acceleration(fit["theta"],design(o["gb"],o["r"],o["mass"],fit["family"],fit["ell"]),o["gb"])

def gradient_controls(gals,cs):
    rows=[]
    for family,ell in SPECS:
        pg,pc=prepared(gals,cs,(family,ell));fun,jac=objective(pg,pc,1,.01)
        n=len(BASES[family]);theta=np.linspace(-.2,.4,n);eps=1e-5
        fd=np.column_stack([(fun(theta+eps*np.eye(n)[i])-fun(theta-eps*np.eye(n)[i]))/(2*eps) for i in range(n)])
        error=float(np.max(abs(fd-jac(theta)))/max(1.,np.max(abs(fd))))
        rows.append(dict(name=f"{family} ell{ell} full Jacobian",relative_error=error,passed=error<1e-7))
    r=np.geomspace(1,100,25);gb=1/r**2
    F=design(gb,r,1e11,"G",0)
    error=float(np.max(abs(F[:,11]-np.tanh(-1))))
    rows.append(dict(name="spherical inverse square slope",error=error,passed=error<1e-12))
    return rows

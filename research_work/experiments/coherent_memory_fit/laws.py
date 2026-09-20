"""CMF-1 effective radial-memory laws. Established tools, hypothetical closure."""
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/"inverse_hitchhiking"))
import data as D
import model as OLD
from campaign import plain,save
BASES={
 "L":["1","x","x2","x3"],
 "R":["1","x","x2","x3","y","xy","y2"],
 "M":["1","x","x2","x3","z","xz","z2"],
 "RM":["1","x","x2","x3","y","xy","y2","z","xz","z2","yz"],
 "RMH":["1","x","x2","x3","y","xy","y2","z","xz","z2","yz","h","xh","zh"]}
SPECS=[(f,0.) for f in ("L","R","M","RM")]+[("RMH",e) for e in (.3,1.,3.)]

def memory(t,x,ell):
    if ell<=0:return x.copy()
    h=np.empty_like(x);h[0]=x[0]
    for i in range(1,len(x)):
        dt=t[i]-t[i-1]
        if dt<=0:raise ValueError("Memory needs strictly increasing radii")
        one=-np.expm1(-dt/ell);e=1-one
        h[i]=e*h[i-1]+one*x[i-1]+(x[i]-x[i-1])/dt*(dt-ell*one)
    return h

def design(gb,r,mass,family,ell):
    gb=np.asarray(gb);r=np.asarray(r)
    x=np.tanh(-np.log(np.maximum(gb*D.CODE_TO_SI,1e-30)/1e-10)/4)
    y=np.tanh(np.log(r/10)/4);z=np.full_like(x,np.tanh(np.log(mass/1e11)/4))
    h=memory(np.log(r),x,ell)-x
    c={"1":np.ones_like(x),"x":x,"x2":x*x,"x3":x**3,
       "y":y,"xy":x*y,"y2":y*y,"z":z,"xz":x*z,"z2":z*z,
       "yz":y*z,"h":h,"xh":x*h,"zh":z*h}
    return np.stack([c[k] for k in BASES[family]],axis=1)

def acceleration(theta,F,gb,jac=False):
    u=np.tanh(F@np.asarray(theta)/8)
    release=1/(1+(np.asarray(gb)*D.CODE_TO_SI/1e-7)**2)
    g=gb*np.exp(8*release*u)
    if not jac:return g
    return g,g[:,None]*release[:,None]*(1-u*u)[:,None]*F

def pressure_matrix(c):
    # Each column is one independent grid-force coefficient in the trapezoid.
    conv=.6*D.MP*1e12/D.KEV_CM3_PA
    cum=cumulative_trapezoid(np.diag(c["ne"]*conv),c["r"],axis=0,initial=0)
    tail=cum[-1]-cum
    index=np.clip(np.searchsorted(c["r"],c["rp"])-1,0,len(c["r"])-2)
    frac=(c["rp"]-c["r"][index])/(c["r"][index+1]-c["r"][index])
    return (1-frac[:,None])*tail[index]+frac[:,None]*tail[index+1]

def pressure_prediction(c,g,dg=None):
    p=c["Pmat"]@g;w=1/c["error"]**2;w/=w.sum()
    boundary=max(0.,float(w@(c["y"]-p)))
    if dg is None:return p+boundary,boundary
    dp=c["Pmat"]@dg
    if boundary>0:dp=dp-(w@dp)[None,:]
    return p+boundary,boundary,dp

def prepared(galaxies,clusters,spec):
    family,ell=spec
    return [{**g,"F":design(g["gb"],g["r"],g["mass"],family,ell)} for g in galaxies],[
        {**c,"F":design(c["gb"],c["r"],c["mass"],family,ell)} for c in clusters]

def objective(gals,clusters,weight,ridge):
    gal=[g for g in gals if g["split"]=="train"]
    cs=[c for c in clusters if c["split"]=="train"]
    F=np.concatenate([g["F"] for g in gal])
    gb=np.concatenate([g["gb"] for g in gal]);r=np.concatenate([g["r"] for g in gal])
    y=np.concatenate([g["y"] for g in gal])
    weights=np.concatenate([np.full(len(g["r"]),1/(20*np.sqrt(len(gal)*len(g["r"])))) for g in gal])
    cnorm=np.sqrt(weight/(10*sum(len(c["y"]) for c in cs)))
    n=F.shape[1];pen=np.sqrt(ridge/n)
    cache={}
    def calc(theta):
        if "theta" in cache and np.array_equal(cache["theta"],theta):return cache["res"],cache["jac"]
        force,df=acceleration(theta,F,gb,True)
        v=np.sqrt(r*force)
        # derivative vanishes at g_b=0; no division by zero.
        dv=np.zeros_like(df);positive=force>0
        dv[positive]=.5*v[positive,None]*df[positive]/force[positive,None]
        res=[(v-y)*weights];jac=[dv*weights[:,None]]
        if weight>0:
            for c in cs:
                force,df=acceleration(theta,c["F"],c["gb"],True)
                p,b,dp=pressure_prediction(c,force,df)
                res.append((p-c["y"])/c["error"]*cnorm)
                jac.append(dp/c["error"][:,None]*cnorm)
        res.append(pen*theta);jac.append(pen*np.eye(n))
        rr=np.concatenate(res);jj=np.concatenate(jac)
        cache.update(theta=theta.copy(),res=rr,jac=jj)
        return rr,jj
    return lambda t:calc(t)[0],lambda t:calc(t)[1]

def evaluate(gals,clusters,force,splits=("train","validation","test"),details=False):
    out={}
    for split in splits:
        gr=[];cr=[]
        for g in gals:
            if g["split"]!=split:continue
            pred=np.sqrt(g["r"]*force(g))
            row=dict(name=g["name"],n=len(pred),mse=float(np.mean((pred-g["y"])**2)),
                chi2=float(np.sum(((pred-g["y"])/g["error"])**2)))
            if details:row.update(radius_kpc=g["r"],observed=g["y"],predicted=pred,error=g["error"],
                distance_Mpc=g["catalog"]["distance"])
            gr.append(row)
        for c in clusters:
            if c["split"]!=split:continue
            pred,b=pressure_prediction(c,force(c))
            row=dict(name=c["name"],n=len(pred),chi2=float(np.sum(((pred-c["y"])/c["error"])**2)),boundary_pressure=b)
            if details:row.update(radius_kpc=c["rp"],observed=c["y"],predicted=pred,error=c["error"])
            cr.append(row)
        gm=np.sqrt(np.mean([r["mse"] for r in gr]));cp=sum(r["chi2"] for r in cr)/sum(r["n"] for r in cr)
        out[split]=dict(galaxy_rmse=gm,galaxy_chi2_per_point=sum(r["chi2"] for r in gr)/sum(r["n"] for r in gr),
            cluster_chi2_per_point=cp,joint_objective=(gm/20)**2+cp/10,galaxies=gr,clusters=cr)
    return out

def candidate_force(fit):
    return lambda obj:acceleration(fit["theta"],design(obj["gb"],obj["r"],obj["mass"],fit["family"],fit["ell"]),obj["gb"])

def optical_law(source,force,reach=9000,n=4096,rmin=.1):
    r=np.geomspace(rmin,reach,n)
    mb=np.interp(r,source["r"],source["mb"],left=source["mb"][0],right=source["mass"])
    inner=r<source["r"][0];mb[inner]=source["mb"][0]*(r[inner]/source["r"][0])**3
    gb=D.G*mb/r**2
    total=force(dict(r=r,gb=gb,mass=source["mass"]))
    ah=total-gb
    def law(q):
        q=np.asarray(q)
        mass=np.interp(q,source["r"],source["mb"],left=source["mb"][0],right=source["mass"])
        base=D.G*mass/q**2
        value=np.interp(np.log(np.minimum(q,reach)),np.log(r),ah)
        return base+value*np.minimum(1.,(reach/q)**2)
    return law

def optics(source,force,reach=9000,n=4096,rmin=.1,order=128,h=.002):
    law=optical_law(source,force,reach,n,rmin)
    obs=D.coma_data()
    result=OLD.shear_shape(law,obs["r"],order,h,tuple(sorted(set((3000.,reach)))))
    unit=50000*result["shape"] # D_l=100,000 kpc, beta=1
    w=1/obs["error"]**2
    beta=float(np.sum(unit*obs["y"]*w)/np.sum(unit*unit*w))
    def score(b):
        p=b*unit
        return dict(beta=b,predicted_shear=p,chi2=float(np.sum(((p-obs["y"])/obs["error"])**2)))
    return dict(reach_kpc=reach,rmin_kpc=rmin,radial_grid=n,alpha=result["alpha"],shape=result["shape"],
        unit_beta_shear=unit,beta_unbounded=beta,
        bounded=score(float(np.clip(beta,0,1))),fixed=[score(b) for b in (.5,.9,1.)],
        shape_only=score(max(0.,beta)))

def controls(gals,clusters):
    rows=[]
    def add(name,value,tolerance):
        rows.append(dict(name=name,value=float(value),limit=tolerance,passed=bool(value<=tolerance)))
    add("same galaxy objects",abs(len(gals)-149),0)
    add("same positive radius rows",abs(sum(len(g["r"]) for g in gals)-3152),0)
    add("same cluster pressure rows",abs(sum(len(c["y"]) for c in clusters)-246),0)
    pg,pc=prepared(gals,clusters,("RMH",1))
    fun,jac=objective(pg,pc,1,.01)
    t=np.linspace(-.2,.4,14);eps=1e-5
    fd=np.column_stack([(fun(t+eps*np.eye(14)[i])-fun(t-eps*np.eye(14)[i]))/(2*eps) for i in range(14)])
    add("full analytic residual Jacobian",np.max(abs(fd-jac(t)))/max(1,np.max(abs(fd))),1e-7)
    for c in clusters:
        g=c["gb"]*(1+.2*np.sin(np.log(c["r"])))
        a=c["Pmat"]@g;b=D.pressure(c,g)
        add("pressure linear map "+c["name"],np.max(abs(a-b))/max(abs(b)),1e-10)
    t=np.linspace(-3,3,51);x=np.sin(t)*.7;ell=.7
    exact=memory(t,x,ell)
    ode=solve_ivp(lambda u,h:(np.interp(u,t,x)-h)/ell,(t[0],t[-1]),[x[0]],t_eval=t,rtol=1e-10,atol=1e-12,max_step=.005)
    add("memory independent ODE",np.max(abs(exact-ode.y[0])),1e-7)
    add("constant memory",np.max(abs(memory(t,np.full_like(t,.4),ell)-.4)),1e-14)
    refine_t=np.sort(np.r_[t,(t[:-1]+t[1:])/2])
    add("piecewise linear memory refinement",np.max(abs(memory(refine_t,np.interp(refine_t,t,x),ell)[::2]-exact)),1e-12)
    theta=np.ones(14)*8;gb=np.full(3,1./D.CODE_TO_SI);r=np.array([1.,10.,100.])
    gg=acceleration(theta,design(gb,r,1e11,"RMH",1),gb)
    add("large acceleration ordinary limit",np.max(abs(gg/gb-1)),1e-12)
    zero=acceleration(np.ones(4),np.ones((1,4)),np.zeros(1))
    add("zero source no generated acceleration",zero[0],0)
    rows.extend(OLD.controls())
    return rows

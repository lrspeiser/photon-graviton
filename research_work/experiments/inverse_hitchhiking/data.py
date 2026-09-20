"""Fixed observed inputs and ordinary-matter reductions for IH-1."""
import io,json,zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
G=4.30091727003628e-6
KPC_M=3.085677581491367e19
CODE_TO_SI=1e6/KPC_M
C=299792.458
MP=1.67262192e-27
MSUN=1.98847e30
KEV_CM3_PA=1.602176634e-10
FILES=[
 ROOT/"companion_wave_test/data/sparc_frozen.json",
 ROOT/"temporal_candidate_audit/data/SPARC_Lelli2016c.mrt",
 ROOT/"temporal_candidate_audit/data/Rotmod_LTG.zip",
 ROOT/"research_work/results/path-memory/cl2-inputs-xcop-profiles.json",
 ROOT/"research_work/results/cluster-observation-readiness/kubo-figure-data.json"]
def read(p):return json.loads(p.read_text(encoding="utf8"))
def galaxies():
    splits=read(FILES[0])["split"]
    catalog={}
    for line in FILES[1].read_text(encoding="utf8").splitlines():
        f=line.split()
        if len(f)!=19:continue
        try:catalog[f[0]]=dict(distance=float(f[2]),luminosity=float(f[7])*1e9,
                              rd=float(f[11]),MHI=float(f[13])*1e9)
        except ValueError:pass
    out=[]
    with zipfile.ZipFile(FILES[2]) as z:
        for split,names in splits.items():
            for name in names:
                a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(name+"_rotmod.dat"))))
                a=a[a[:,0]>0]
                cat=catalog[name];r=a[:,0]
                gb=np.maximum(a[:,3]*np.abs(a[:,3])+.5*a[:,4]**2+.7*a[:,5]**2,0)/r
                out.append(dict(name=name,split=split,r=r,y=a[:,1],error=a[:,2],
                    gb=gb,mass=.5*cat["luminosity"]+1.33*cat["MHI"],catalog=cat))
    assert len(out)==149
    return out
def loginterp(x,xp,yp,inner=None):
    x=np.asarray(x);xp=np.asarray(xp);yp=np.asarray(yp)
    result=np.exp(np.interp(np.log(x),np.log(xp),np.log(yp)))
    if inner is not None:
        mask=x<xp[0];result[mask]=yp[0]*(x[mask]/xp[0])**inner
    return result
def cluster_inputs(n=600):
    raw=read(FILES[3])["clusters"]
    xs=np.geomspace(.01,2,60);fractions=[]
    for c in raw.values():
        if "stellar_mass" in c:
            r500=c["header"]["R500_kpc"];sm=c["stellar_mass"];gm=c["gas_mass"]
            fractions.append(loginterp(xs*r500,sm["radius_kpc"],sm["Mstar"],3)/
                             loginterp(xs,gm["RADIUS"],gm["MGAS"],3))
    frac=np.median(fractions,axis=0)
    out=[]
    for i,name in enumerate(sorted(raw)):
        c=raw[name];r500=c["header"]["R500_kpc"];d=c["density"]
        rm=np.sqrt(np.array(d["r_in_kpc"])*d["r_out_kpc"]);ne=np.array(d["ne_cm3"])
        rout=min(d["r_out_kpc"][-1],c["pressure"]["sz"]["r_over_R500"][-1]*r500)
        r=np.geomspace(1,rout,n)
        lne=np.interp(np.log(r),np.log(rm),np.log(ne),left=np.log(ne[0]))
        slope=np.log(ne[-1]/ne[-2])/np.log(rm[-1]/rm[-2])
        outer=r>rm[-1];lne[outer]=np.log(ne[-1])+slope*np.log(r[outer]/rm[-1])
        ne=np.exp(lne)
        rho=1.14*MP*ne*1e6*KPC_M**3/MSUN
        mint=cumulative_trapezoid(4*np.pi*r*r*rho,r,initial=0)+4*np.pi/3*r[0]**3*rho[0]
        gm=c["gas_mass"];rr=np.array(gm["RADIUS"])*r500;mm=np.array(gm["MGAS"])
        mg=loginterp(r,rr,mm);inner=r<rr[0]
        mg[inner]=mint[inner]*(mm[0]/np.interp(rr[0],r,mint))
        if "stellar_mass" in c:
            sm=c["stellar_mass"];ms=loginterp(r,sm["radius_kpc"],sm["Mstar"],3)
            star_source="measured stellar profile"
        else:
            ms=np.interp(np.log(r/r500),np.log(xs),frac,left=frac[0],right=frac[-1])*mg
            star_source="frozen median fraction of seven measured stellar profiles"
        mb=mg+ms;pts=[];p500=c["pressure"]["P500_keV_cm3"]
        for kind in ("xray","sz"):
            p=c["pressure"][kind]
            for rr,pp,ee in zip(p["r_over_R500"],p["P_over_P500"],p["err"]):
                if r[0]<=rr*r500<=rout:pts.append((rr*r500,pp*p500,ee*p500,kind))
        pts.sort()
        out.append(dict(name=name,split="train" if i%4 in (0,1) else ("validation" if i%4==2 else "test"),
            r=r,ne=ne,mb=mb,gb=G*mb/r**2,mass=mb[-1],rp=np.array([p[0] for p in pts]),
            y=np.array([p[1] for p in pts]),error=np.array([p[2] for p in pts]),
            kinds=[p[3] for p in pts],star_source=star_source,
            gas_integral_check=float(np.interp(r500,r,mint)/np.interp(r500,r,mg))))
    return out
def pressure(c,g):
    # Electron pressure; ne cm^-3, g (km/s)^2/kpc, r kpc.
    conv=.6*MP*1e12/KEV_CM3_PA
    integral=-cumulative_trapezoid((c["ne"]*g*conv)[::-1],c["r"][::-1],initial=0)[::-1]
    return np.interp(c["rp"],c["r"],integral)
def pressure_residual(c,g):
    pred=pressure(c,g)
    w=1/c["error"]**2
    boundary=max(0.,float(np.sum(w*(c["y"]-pred))/w.sum()))
    return (pred+boundary-c["y"])/c["error"],pred+boundary,boundary
def coma_source(ne0,mstar,n=4096):
    r=np.geomspace(.1,3000,n)
    rho=1.17*MP*ne0*1e6*KPC_M**3/MSUN*(1+(r/296)**2)**(-1.125)
    mg=cumulative_trapezoid(4*np.pi*r*r*rho,r,initial=0)+4*np.pi/3*r[0]**3*rho[0]
    mb=mg*(1+mstar/mg[-1])
    return dict(r=r,mb=mb,mass=mb[-1],ne0=ne0,mstar=mstar)
def coma_data():
    d=read(FILES[4]);rows=d["rows"]
    return dict(r=np.array([r["published_radius_h_inverse_Mpc"] for r in rows])*1000/.7,
        y=np.array([r["shear_t"] for r in rows]),error=np.array([r["plotted_sigma_t"] for r in rows]),
        status=d["status"],source=d["source"])

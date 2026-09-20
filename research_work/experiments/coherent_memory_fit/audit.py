"""Independent evidence, formula, geometry and provenance audit."""
import hashlib,json,subprocess,warnings
import numpy as np
from scipy.integrate import quad,solve_ivp
import gradient_laws as L
import run as C
D=L.D;HERE=L.HERE
def read(p):return json.loads(p.read_text(encoding="utf8"))
def main():
    rows=[]
    def check(name,value,limit=0):
        rows.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
    evidence_count=0;pinned=0
    for file in ("evidence-sha256.json","gradient-evidence-sha256.json"):
        for path,expected in read(HERE/file).items():
            check("evidence "+path,0 if C.digest(HERE/path)==expected else 1);evidence_count+=1
    for folder in ("evidence-v1","evidence-gradient-v1"):
        manifest=read(HERE/folder/"manifest.json");head=manifest["git_head"]
        for path,expected in manifest["hashes"].items():
            actual=C.digest(L.ROOT/path)
            check("current input/source "+folder+" "+path,0 if actual==expected else 1)
            original=subprocess.check_output(["git","show",head+":"+path],cwd=L.ROOT)
            check("pinned Git "+folder+" "+path,0 if hashlib.sha256(original).hexdigest()==expected else 1)
            pinned+=1
        controls=read(HERE/folder/"controls.json")
        check("saved preliminary controls "+folder,sum(not r["passed"] for r in controls))
        numerical=read(HERE/folder/"numerical-checks.json")
        check("saved refinements "+folder,sum(not r.get("passed",True) for r in numerical))
    sel=read(HERE/"evidence-gradient-v1"/"selection-before-test.json")["selected"]
    pred=read(HERE/"evidence-gradient-v1"/"predictions.json");gals=D.galaxies();cs=C.clusters()
    # Independently write the selected local polynomial, without design()/acceleration().
    fit=sel["joint"]
    if fit["family"]!="RM":raise AssertionError("Independent polynomial audit needs update for changed selection")
    theta=np.array(fit["theta"])
    def independent(o):
        gb=np.asarray(o["gb"]);r=np.asarray(o["r"])
        x=np.tanh(-np.log(np.maximum(gb*D.CODE_TO_SI,1e-30)/1e-10)/4)
        y=np.tanh(np.log(r/10)/4);z=np.tanh(np.log(o["mass"]/1e11)/4)
        values=[np.ones_like(x),x,x*x,x*x*x,y,x*y,y*y,
                np.full_like(x,z),x*z,np.full_like(x,z*z),y*z]
        p=sum(a*b for a,b in zip(theta,values))
        return gb*np.exp(8*np.tanh(p/8)/(1+(gb*D.CODE_TO_SI/1e-7)**2))
    max_velocity=0.;max_pressure=0.
    for g in gals:
        row=next(r for r in pred["joint"][g["split"]]["galaxies"] if r["name"]==g["name"])
        check("fixed distance "+g["name"],abs(row["distance_Mpc"]-g["catalog"]["distance"]))
        max_velocity=max(max_velocity,float(np.max(abs(np.sqrt(g["r"]*independent(g))-row["predicted"]))))
    check("independent selected velocities km/s",max_velocity,1e-10)
    for c in cs:
        p=D.pressure(c,independent(c));w=1/c["error"]**2
        p+=max(0.,float(np.sum(w*(c["y"]-p))/w.sum()))
        row=next(r for r in pred["joint"][c["split"]]["clusters"] if r["name"]==c["name"])
        max_pressure=max(max_pressure,float(np.max(abs(p-row["predicted"]))))
    check("independent selected pressure",max_pressure,1e-12)
    # An independent ODE solver evaluates the winning galaxy-only memory.
    memory_fit=sel["galaxy"];ell=memory_fit["ell"];memory_error=0.
    for g in [gals[i] for i in (0,20,60,100,148)]:
        t=np.log(g["r"]);x=np.tanh(-np.log(np.maximum(g["gb"]*D.CODE_TO_SI,1e-30)/1e-10)/4)
        ode=solve_ivp(lambda u,h:(np.interp(u,t,x)-h)/ell,(t[0],t[-1]),[x[0]],t_eval=t,rtol=1e-11,atol=1e-12,max_step=.01)
        error=float(np.max(abs(ode.y[0]-L.memory(t,x,ell))))
        memory_error=max(memory_error,error)
    check("selected memory independent ODE",memory_error,1e-7)
    lens=read(HERE/"evidence-gradient-v1"/"lensing.json")
    alpha_error=0.;integration_warnings=[]
    for row in [r for r in lens["models"] if r["model"]=="joint" and r["primary"]]:
        src=D.coma_source(row["ne0"],row["mstar"])
        law=L.optical_law(src,independent)
        for b,archived in zip(D.coma_data()["r"],row["alpha"]):
            edges=[0.]+sorted(np.arccos(b/r) for r in (3000.,9000.) if r>b)+[np.pi/2]
            value=0.
            for lo,hi in zip(edges[:-1],edges[1:]):
                with warnings.catch_warnings(record=True) as captured:
                    warnings.simplefilter("always")
                    v,_=quad(lambda angle:float(law(np.array([b/np.cos(angle)]))[0])*b/np.cos(angle),
                             lo,hi,epsabs=1e-5,epsrel=2e-7,limit=500)
                integration_warnings.extend(str(w.message) for w in captured)
                value+=v
            value*=4/D.C**2
            alpha_error=max(alpha_error,abs(value/archived-1))
        unit=np.array(row["shape"])*100000/2
        check("fixed static lens scale "+str(row["ne0"]),np.max(abs(unit-np.array(row["unit_beta_shear"]))),1e-16)
        b=float(np.clip(row["beta_unbounded"],0,1))
        chi=np.sum(((b*unit-D.coma_data()["y"])/D.coma_data()["error"])**2)
        check("bounded lens score "+str(row["ne0"]),abs(chi-row["bounded"]["chi2"]),1e-10)
    check("adaptive deflection quadrature",alpha_error,1e-4)
    refs=read(HERE/"evidence-v1"/"references.json")
    a0=refs["MOND_a0_SI"]
    mf=lambda o:.5*(o["gb"]+np.sqrt(o["gb"]**2+4*a0/D.CODE_TO_SI*o["gb"]))
    for split in ("train","validation","test"):
        group=[g for g in gals if g["split"]==split]
        rmse=np.sqrt(np.mean([np.mean((np.sqrt(g["r"]*mf(g))-g["y"])**2) for g in group]))
        check("matched MOND "+split,abs(rmse-refs["scores"]["MOND"][split]["galaxy_rmse"]),1e-12)
    totals=[read(HERE/f/"complete.json") for f in ("evidence-v1","evidence-gradient-v1")]
    check("132 fits",abs(sum(t["fit_count"] for t in totals)-132))
    check("396 starts",abs(sum(t["attempts"] for t in totals)-396))
    result=dict(passed=all(r["passed"] for r in rows),evidence_files=evidence_count,pinned_source_input_checks=pinned,
        unchanged_distances=len(gals),independent_velocity_error_kms=max_velocity,independent_pressure_error=max_pressure,
        independent_memory_error=memory_error,independent_optics_relative_error=alpha_error,
        adaptive_integration_warnings=integration_warnings,
        warning_scope="Requested adaptive tolerance may not be achieved at piecewise interpolation knots; compare independently at the declared 1e-4 audit threshold.",checks=rows)
    L.save(HERE/"audit.json",result)
    print(json.dumps({k:v for k,v in result.items() if k!="checks"},indent=2))
    if not result["passed"]:raise AssertionError("Evidence audit failed")
if __name__=="__main__":main()

"""Fit the declared IH-1 response families; preserve every attempt."""
import hashlib,json,platform,subprocess,time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import least_squares
import data as D
import model as M
HERE=Path(__file__).resolve().parent
def plain(x):
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [plain(v) for v in x]
    if isinstance(x,(np.bool_,bool)):return bool(x)
    if isinstance(x,(np.integer,int)):return int(x)
    if isinstance(x,(np.floating,float)):
        if not np.isfinite(x):raise ValueError("nonfinite result")
        return float(x)
    return x
def save(path,x):path.write_text(json.dumps(plain(x),indent=2,allow_nan=False)+"\n",encoding="utf8",newline="\n")
def pack(gals,split):
    chosen=[g for g in gals if g["split"]==split]
    return dict(r=np.concatenate([g["r"] for g in chosen]),gb=np.concatenate([g["gb"] for g in chosen]),
        y=np.concatenate([g["y"] for g in chosen]),mass=np.concatenate([np.full(len(g["r"]),g["mass"]) for g in chosen]),
        weight=np.concatenate([np.full(len(g["r"]),1/(20*np.sqrt(len(chosen)*len(g["r"])))) for g in chosen]))
def galaxy_score(gals,fit,split,details=False):
    records=[]
    for g in gals:
        if g["split"]!=split:continue
        ah=0 if fit is None else M.response(fit,g["gb"],g["r"],g["mass"])
        pred=np.sqrt(g["r"]*(g["gb"]+ah))
        row=dict(name=g["name"],rmse=float(np.sqrt(np.mean((pred-g["y"])**2))),n=len(pred),
            chi2=float(np.sum(((pred-g["y"])/g["error"])**2)))
        if details:row.update(radius_kpc=g["r"],observed_kms=g["y"],error_kms=g["error"],predicted_kms=pred,
            baryon_kms=np.sqrt(g["r"]*g["gb"]),distance_Mpc=g["catalog"]["distance"],source_mass_proxy_Msun=g["mass"])
        records.append(row)
    return dict(rmse=float(np.sqrt(np.mean([r["rmse"]**2 for r in records]))),
        galaxies=len(records),points=sum(r["n"] for r in records),
        chi2_per_point=sum(r["chi2"] for r in records)/sum(r["n"] for r in records),rows=records)
def cluster_score(clusters,fit,split,details=False):
    records=[]
    for c in clusters:
        if c["split"]!=split:continue
        ah=0 if fit is None else M.response(fit,c["gb"],c["r"],c["mass"])
        resid,pred,bound=D.pressure_residual(c,c["gb"]+ah)
        row=dict(name=c["name"],n=len(resid),chi2=float(np.sum(resid**2)),boundary_pressure=bound)
        if details:row.update(radius_kpc=c["rp"],observed_pressure=c["y"],error=c["error"],prediction=pred,
            source_mass_proxy_Msun=c["mass"],star_source=c["star_source"])
        records.append(row)
    total=sum(r["chi2"] for r in records);points=sum(r["n"] for r in records)
    return dict(chi2_per_point=total/points,chi2=total,points=points,clusters=len(records),rows=records)
def scores(gals,clusters,fit,splits,details=False):
    out={}
    for split in splits:
        g=galaxy_score(gals,fit,split,details);c=cluster_score(clusters,fit,split,details)
        out[split]=dict(galaxy=g,cluster=c,joint_objective=(g["rmse"]/20)**2+c["chi2_per_point"]/10)
    return out
def fit_one(gpack,ctrain,family,capacity,mode,rng):
    gf=M.features(gpack["gb"],gpack["r"],gpack["mass"],family)
    cf=[M.features(c["gb"],c["r"],c["mass"],family) for c in ctrain]
    cnorm=np.sqrt(sum(len(c["y"]) for c in ctrain)*10)
    def residual(theta):
        gp=np.sqrt(gpack["r"]*(gpack["gb"]+M.extra(theta,capacity,gf)))
        r=[(gp-gpack["y"])*gpack["weight"]]
        if mode=="joint":
            for c,F in zip(ctrain,cf):
                rr,_,_=D.pressure_residual(c,c["gb"]+M.extra(theta,capacity,F))
                r.append(rr/cnorm)
        return np.concatenate(r)
    lo,hi=M.bounds(family)
    x=np.zeros(len(lo));x[0]=np.log(1e-10/capacity);x[1]=.5
    starts=[x]+[np.clip(x+rng.normal(0,.5,len(x)),lo+1e-6,hi-1e-6) for _ in range(2)]
    attempts=[]
    for start in starts:
        t=time.monotonic()
        res=least_squares(residual,start,bounds=(lo,hi),xtol=1e-9,ftol=1e-9,gtol=1e-9,max_nfev=500)
        attempts.append(dict(initial=start,theta=res.x,objective=float(np.sum(res.fun**2)),success=res.success,
            status=res.status,message=res.message,nfev=res.nfev,optimality=res.optimality,
            bound_parameters=[k for k,v,a,b in zip(M.FAMILIES[family],res.x,lo,hi) if min(v-a,b-v)<1e-4],
            seconds=time.monotonic()-t))
    best=min(attempts,key=lambda a:a["objective"])
    return dict(id=f"{mode}-{family}-cap{capacity:.0e}",mode=mode,family=family,capacity=capacity,
        theta=best["theta"],parameters=dict(zip(M.FAMILIES[family],best["theta"])),objective=best["objective"],
        success=best["success"],bound_parameters=best["bound_parameters"],attempts=attempts)
def main():
    root=HERE/"evidence-v1";root.mkdir(exist_ok=False)
    files=D.FILES+list(HERE.glob("*.py"))+list(HERE.glob("*.md"))
    save(root/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=D.ROOT,text=True).strip(),
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        inputs_and_sources={p.relative_to(D.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}))
    save(root/"controls.json",M.controls())
    gals=D.galaxies();clusters=D.cluster_inputs()
    save(root/"data-inventory.json",dict(galaxies=[dict(name=g["name"],split=g["split"],rows=len(g["r"]),distance_Mpc=g["catalog"]["distance"]) for g in gals],
        clusters=[dict(name=c["name"],split=c["split"],rows=len(c["y"]),star_source=c["star_source"],gas_integral_check=c["gas_integral_check"]) for c in clusters],
        historical_exposure=True,no_dark_matter=True,no_expansion=True,no_distance_fitting=True))
    gp=pack(gals,"train");ct=[c for c in clusters if c["split"]=="train"]
    fits=[];rng=np.random.default_rng(190921)
    for family in M.FAMILIES:
        for cap in M.CAPS:
            for mode in ("galaxy","joint"):
                fit=fit_one(gp,ct,family,cap,mode,rng)
                fit["scores"]=scores(gals,clusters,fit,("train","validation"))
                fits.append(fit);save(root/"fits.json",fits)
                v=fit["scores"]["validation"]
                print(f"{fit['id']}: val galaxy {v['galaxy']['rmse']:.3f} km/s, cluster {v['cluster']['chi2_per_point']:.3f}/point, success {fit['success']}",flush=True)
    selected={}
    for mode in ("galaxy","joint"):
        eligible=[f for f in fits if f["mode"]==mode and f["success"]]
        if not eligible:raise RuntimeError("No optimizer-success candidate in "+mode)
        def key(f):
            v=f["scores"]["validation"]
            score=v["joint_objective"] if mode=="joint" else v["galaxy"]["rmse"]
            return (score,len(f["theta"]),f["id"])
        selected[mode]=min(eligible,key=key)
    # Persist selection before any test-block or shear calculation.
    save(root/"selection-before-test.json",{k:dict(id=f["id"],family=f["family"],capacity=f["capacity"],theta=f["theta"],
        criterion="validation only",test_not_used=True,coma_not_used=True) for k,f in selected.items()})
    save(root/"selected-predictions.json",{mode:scores(gals,clusters,fit,("train","validation","test"),True) for mode,fit in selected.items()})
    save(root/"baryons-only.json",scores(gals,clusters,None,("train","validation","test"),True))
    inverse=[]
    for g in gals:
        required=(g["y"]**2/g["r"]-g["gb"])*D.CODE_TO_SI
        inverse.append(dict(name=g["name"],split=g["split"],radius_kpc=g["r"],required_extra_acceleration_SI=required,
            required_fraction={f"{cap:.0e}":required/cap for cap in M.CAPS},
            inadmissible_negative=int(np.sum(required<0)),
            capped_counts={f"{cap:.0e}":int(np.sum(required>=cap)) for cap in M.CAPS}))
    save(root/"inverse-requirements.json",inverse)
    fit=selected["joint"];data=D.coma_data();light=[]
    for ne0,mstar in ((2.5e-3,.5e13),(4.5e-3,2e13)):
        src=D.coma_source(ne0,mstar)
        for reach in (3000.,9000.,30000.):
            for eta in (0.,.25,1.,4.):
                g=M.optical_force(src,fit,reach,eta)
                projected=M.shear_shape(g,data["r"],128,.002,(3000.,reach) if reach!=3000 else (reach,))
                result=M.shape_fit(projected["shape"],data["y"],data["error"])
                light.append(dict(ne0=ne0,mstar=mstar,reach_kpc=reach,eta=eta,primary=reach==9000 and eta==1,
                    alpha=projected["alpha"],shear_shape=projected["shape"],fit=result))
    save(root/"coma.json",dict(data=data,models=light,selected_law=fit["id"],
        status="Exposed six-bin shape-only fit; amplitude degenerate with unknown source geometry; no absolute lensing validation."))
    fine=D.cluster_inputs(1200);numerical=[]
    for old,new in zip(clusters,fine):
        g0=old["gb"]+M.extra(fit["theta"],fit["capacity"],M.features(old["gb"],old["r"],old["mass"],fit["family"]))
        g1=new["gb"]+M.extra(fit["theta"],fit["capacity"],M.features(new["gb"],new["r"],new["mass"],fit["family"]))
        p0=D.pressure(old,g0);p1=D.pressure(new,g1)
        err=float(np.max(abs(p0-p1))/np.max(abs(p1)))
        numerical.append(dict(test="pressure grid",name=old["name"],relative_max_error=err,passed=err<.005))
    for row in [r for r in light if r["primary"]]:
        src=D.coma_source(row["ne0"],row["mstar"]);g=M.optical_force(src,fit,9000,1)
        ref=M.shear_shape(g,data["r"],256,.001,(3000.,9000.))
        error=float(np.max(abs(np.array(row["shear_shape"])-ref["shape"]))/np.max(abs(ref["shape"])))
        numerical.append(dict(test="light quadrature/derivative",ne0=row["ne0"],relative_max_error=error,passed=error<.005))
    save(root/"numerical-checks.json",numerical)
    save(root/"complete.json",dict(complete=True,fits=len(fits),attempts=sum(len(f["attempts"]) for f in fits),
        numerical_passed=all(r["passed"] for r in numerical),
        operational_validation_target=fit["scores"]["validation"]["galaxy"]["rmse"]<=22 and fit["scores"]["validation"]["cluster"]["chi2_per_point"]<=10))
    save(HERE/"evidence-sha256.json",{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob("*") if p.is_file()})
    print("Completed: "+str({m:f["id"] for m,f in selected.items()}),flush=True)
if __name__=="__main__":main()

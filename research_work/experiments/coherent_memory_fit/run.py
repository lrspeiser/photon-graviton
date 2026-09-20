"""Execute the predeclared CMF-1 campaign; never overwrite evidence."""
import hashlib,json,platform,subprocess,time
from datetime import datetime,timezone
import numpy as np
import scipy
from scipy.optimize import least_squares,minimize_scalar
import laws as L
D=L.D;HERE=L.HERE;save=L.save

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clusters(n=600):
    cs=D.cluster_inputs(n)
    for c in cs:c["Pmat"]=L.pressure_matrix(c)
    return cs

def fit_one(gals,cs,family,ell,weight,ridge,rng):
    residual,jac=L.objective(gals,cs,weight,ridge)
    n=len(L.BASES[family]);start=np.zeros(n);start[:2]=[.6,1.]
    starts=[start]+[start+rng.normal(0,.4,n) for _ in range(2)]
    attempts=[]
    for st in starts:
        t=time.monotonic()
        res=least_squares(residual,st,jac=jac,bounds=(-8.,8.),max_nfev=400,
            ftol=1e-9,xtol=1e-9,gtol=1e-9)
        attempts.append(dict(initial=st,theta=res.x,objective=float(res.fun@res.fun),success=res.success,
            status=res.status,message=res.message,nfev=res.nfev,optimality=res.optimality,seconds=time.monotonic()-t))
    successful=[a for a in attempts if a["success"]]
    best=min(successful or attempts,key=lambda a:a["objective"])
    return dict(id=f"{family}-ell{ell:g}-w{weight:g}-ridge{ridge:g}",family=family,ell=ell,
        weight=weight,ridge=ridge,n_parameters=n,theta=best["theta"],basis=L.BASES[family],
        objective=best["objective"],success=bool(successful),attempts=attempts,
        bound_parameters=[name for name,val in zip(L.BASES[family],best["theta"]) if abs(val)>7.9999])

def bootstrap(a,b):
    ra=a["test"]["galaxies"];rb=b["test"]["galaxies"]
    assert [r["name"] for r in ra]==[r["name"] for r in rb]
    av=np.array([r["mse"] for r in ra]);bv=np.array([r["mse"] for r in rb])
    rng=np.random.default_rng(200926);idx=rng.integers(0,len(av),(10000,len(av)))
    delta=np.sqrt(av[idx].mean(axis=1))-np.sqrt(bv[idx].mean(axis=1))
    return dict(delta_rmse_kms=float(np.sqrt(av.mean())-np.sqrt(bv.mean())),
        percentile95=np.quantile(delta,[.025,.975]),resample_fraction_lower=float(np.mean(delta<0)),
        galaxies=len(av),resamples=len(idx),
        interpretation="Paired object bootstrap conditional on exposed data and fixed selected models; not post-selection or fresh-test evidence.")

def main():
    out=HERE/"evidence-v1";out.mkdir(exist_ok=False)
    files=D.FILES+list(HERE.glob("*.py"))+list(HERE.glob("*.md"))+[
        HERE.parent/"inverse_hitchhiking"/f for f in ("data.py","model.py","campaign.py","summary.json")]
    head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=L.ROOT,text=True).strip()
    save(out/"manifest.json",dict(git_head=head,utc=datetime.now(timezone.utc).isoformat(),
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        hashes={p.relative_to(L.ROOT).as_posix():digest(p) for p in files}))
    gals=D.galaxies();cs=clusters()
    save(out/"inventory.json",dict(galaxies=[dict(name=g["name"],split=g["split"],n=len(g["r"]),
        distance_Mpc=g["catalog"]["distance"]) for g in gals],
        clusters=[dict(name=c["name"],split=c["split"],n=len(c["y"]),star_source=c["star_source"]) for c in cs],
        historical_exposure=True,dark_matter=False,expansion=False,distance_fitting=False))
    controls=L.controls(gals,cs);save(out/"controls.json",controls)
    if not all(r["passed"] for r in controls):raise AssertionError("See preserved controls.json")
    print(f"All {len(controls)} controls pass.",flush=True)
    def mond(a0):
        return lambda o:.5*(o["gb"]+np.sqrt(o["gb"]**2+4*a0/D.CODE_TO_SI*o["gb"]))
    train=[g for g in gals if g["split"]=="train"]
    def loss(loga):
        force=mond(10.**loga)
        return np.mean([np.mean((np.sqrt(g["r"]*force(g))-g["y"])**2) for g in train])
    opt=minimize_scalar(loss,bounds=(-12.,-8.),method="bounded",options={"xatol":1e-11})
    if not opt.success:raise RuntimeError("MOND benchmark optimizer failed")
    a0=10.**opt.x
    old=D.read(HERE.parent/"inverse_hitchhiking"/"summary.json")["final_formula"]
    refs={"ordinary":lambda o:o["gb"],"MOND":mond(a0),
          "IH1":lambda o:o["gb"]+L.OLD.response(old,o["gb"],o["r"],o["mass"])}
    # At this stage only training and validation targets are scored.
    validation_refs={k:L.evaluate(gals,cs,f,("train","validation")) for k,f in refs.items()}
    save(out/"reference-before-test.json",dict(MOND_a0_SI=a0,n_parameters=1,scores=validation_refs))
    rng=np.random.default_rng(200919);fits=[]
    for family,ell in L.SPECS:
        pg,pc=L.prepared(gals,cs,(family,ell))
        for ridge in (0.,.001,.01):
            for weight in (0.,.1,1.,10.):
                fit=fit_one(pg,pc,family,ell,weight,ridge,rng)
                fit["scores"]=L.evaluate(pg,pc,lambda o:L.acceleration(fit["theta"],o["F"],o["gb"]),("train","validation"))
                fits.append(fit);save(out/"fits.json",fits)
                v=fit["scores"]["validation"]
                print(f"{len(fits)}/84 {fit['id']}: validation galaxy {v['galaxy_rmse']:.3f}, pressure {v['cluster_chi2_per_point']:.3f}, success {fit['success']}",flush=True)
    eligible=[f for f in fits if f["success"]]
    def key(f):return (f["scores"]["validation"]["joint_objective"],f["n_parameters"],f["id"])
    joint=[f for f in eligible if f["weight"]>0]
    selected={"joint":min(joint,key=key),
        "galaxy":min((f for f in eligible if f["weight"]==0),key=lambda f:(f["scores"]["validation"]["galaxy_rmse"],f["n_parameters"],f["id"]))}
    gate=[f for f in joint if f["scores"]["validation"]["galaxy_rmse"]<=validation_refs["MOND"]["validation"]["galaxy_rmse"]
          and f["scores"]["validation"]["cluster_chi2_per_point"]<=10]
    if gate:selected["gated"]=min(gate,key=key)
    save(out/"selection-before-test.json",dict(selected=selected,gate_candidates=len(gate),
        selection_uses_test=False,selection_uses_coma=False))
    predictions={k:L.evaluate(gals,cs,L.candidate_force(f),details=True) for k,f in selected.items()}
    references={k:L.evaluate(gals,cs,f,details=True) for k,f in refs.items()}
    save(out/"predictions.json",predictions);save(out/"references.json",dict(MOND_a0_SI=a0,scores=references))
    save(out/"paired-bootstrap.json",{k:bootstrap(p,references["MOND"]) for k,p in predictions.items()})
    print("Selections saved before test and optics: "+str({k:f["id"] for k,f in selected.items()}),flush=True)
    lights=[]
    forces={**refs,**{k:L.candidate_force(f) for k,f in selected.items()}}
    for name,force in forces.items():
        for ne0,mstar in ((2.5e-3,.5e13),(4.5e-3,2e13)):
            src=D.coma_source(ne0,mstar)
            for reach in (3000.,9000.,30000.):
                row=L.optics(src,force,reach)
                row.update(model=name,ne0=ne0,mstar=mstar,primary=reach==9000.)
                lights.append(row)
    save(out/"lensing.json",dict(distance_Mpc=100.,geometry="fixed static Euclidean approximate cross-catalog geometry",
        source_beta="unknown, physically bounded to [0,1]; not derived from redshift",
        data=D.coma_data(),models=lights,
        limitations="Six figure bins; no full covariance, source weights, bin edges or new angular reduction; conditional envelope, not absolute lens validation."))
    fine=clusters(1200);checks=[];force=L.candidate_force(selected["joint"])
    for c,cf in zip(cs,fine):
        coarse=D.pressure(c,force(c));refined=D.pressure(cf,force(cf))
        err=float(np.max(abs(coarse-refined))/np.max(abs(refined)))
        checks.append(dict(test="pressure radial refinement",name=c["name"],relative_error=err,passed=err<.005))
    for row in [r for r in lights if r["model"]=="joint" and r["primary"]]:
        src=D.coma_source(row["ne0"],row["mstar"],8192)
        refined=L.optics(src,force,9000,8192,.1,256,.001)
        err=float(np.max(abs(np.array(row["unit_beta_shear"])-refined["unit_beta_shear"]))/np.max(abs(refined["unit_beta_shear"])))
        checks.append(dict(test="combined source-memory-optical refinement",ne0=row["ne0"],relative_error=err,passed=err<.005))
        for rmin in (.01,1.):
            r=L.optics(src,force,9000,8192,rmin,256,.001)
            change=float(np.max(abs(r["unit_beta_shear"]-refined["unit_beta_shear"]))/np.max(abs(refined["unit_beta_shear"])))
            checks.append(dict(test="inner memory boundary sensitivity",ne0=row["ne0"],rmin=rmin,relative_change=change,
                note="Sensitivity, not a convergence pass/fail",bounded_chi2=r["bounded"]["chi2"]))
    save(out/"numerical-checks.json",checks)
    main=predictions["joint"];ref=references["MOND"]
    matched=all(main[s]["galaxy_rmse"]<ref[s]["galaxy_rmse"] for s in ("validation","test"))
    bootstrap_primary=bootstrap(main,ref)
    save(out/"complete.json",dict(campaign_executed=True,fit_count=len(fits),
        attempts=sum(len(f["attempts"]) for f in fits),
        optimizer_success_fits=sum(f["success"] for f in fits),
        optimizer_success_attempts=sum(a["success"] for f in fits for a in f["attempts"]),
        numerical_passed=all(r.get("passed",True) for r in checks),
        point_score_beats_MOND=matched,paired95_entirely_better=bool(bootstrap_primary["percentile95"][1]<0),
        cluster_validation_under10=main["validation"]["cluster_chi2_per_point"]<=10,
        complete_solution=False,absolute_lensing_validated=False,derived_energy_funding=False))
    save(HERE/"evidence-sha256.json",{p.relative_to(HERE).as_posix():digest(p) for p in sorted(out.glob("*.json"))})
    print("Completed campaign and evidence digests.",flush=True)
if __name__=="__main__":main()

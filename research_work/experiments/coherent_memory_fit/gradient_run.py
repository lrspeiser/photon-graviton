"""Execute the separately declared CMF-2, retaining prior results."""
import hashlib,json,platform,subprocess,time
from datetime import datetime,timezone
import numpy as np
import scipy
from scipy.optimize import least_squares
import gradient_laws as L
import run as C
D=L.D;HERE=L.HERE;save=L.save
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


def main():
    out=HERE/"evidence-gradient-v1";out.mkdir(exist_ok=False)
    files=D.FILES+list(HERE.glob("*.py"))+list(HERE.glob("*.md"))+[HERE/"evidence-v1"/f for f in
        ("fits.json","references.json","reference-before-test.json","manifest.json")]
    head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=L.ROOT,text=True).strip()
    save(out/"manifest.json",dict(git_head=head,utc=datetime.now(timezone.utc).isoformat(),
        python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        hashes={p.relative_to(L.ROOT).as_posix():C.digest(p) for p in files}))
    gals=D.galaxies();cs=C.clusters()
    controls=L.gradient_controls(gals,cs);save(out/"controls.json",controls)
    if not all(r["passed"] for r in controls):raise AssertionError("Gradient controls failed")
    print(f"All {len(controls)} gradient controls pass.",flush=True)
    rng=np.random.default_rng(200927);fits=[]
    for family,ell in L.SPECS:
        pg,pc=L.prepared(gals,cs,(family,ell))
        for ridge in (0.,.001,.01):
            for weight in (0.,.1,1.,10.):
                fit=fit_one(pg,pc,family,ell,weight,ridge,rng)
                fit["scores"]=L.evaluate(pg,pc,lambda o:L.acceleration(fit["theta"],o["F"],o["gb"]),("train","validation"))
                fits.append(fit);save(out/"fits.json",fits)
                v=fit["scores"]["validation"]
                print(f"{len(fits)}/48 {fit['id']}: validation galaxy {v['galaxy_rmse']:.3f}, pressure {v['cluster_chi2_per_point']:.3f}, success {fit['success']}",flush=True)
    old=D.read(HERE/"evidence-v1"/"fits.json")
    eligible=[f for f in old+fits if f["success"]]
    def key(f):return (f["scores"]["validation"]["joint_objective"],f["n_parameters"],f["id"])
    selected={"joint":min((f for f in eligible if f["weight"]>0),key=key),
        "galaxy":min((f for f in eligible if f["weight"]==0),key=lambda f:(f["scores"]["validation"]["galaxy_rmse"],f["n_parameters"],f["id"]))}
    best_new=min((f for f in fits if f["success"] and f["weight"]>0),key=key)
    refs=D.read(HERE/"evidence-v1"/"references.json")["scores"]
    gate=[f for f in eligible if f["weight"]>0 and f["scores"]["validation"]["galaxy_rmse"]<=refs["MOND"]["validation"]["galaxy_rmse"]
          and f["scores"]["validation"]["cluster_chi2_per_point"]<=10]
    if gate:selected["gated"]=min(gate,key=key)
    save(out/"selection-before-test.json",dict(selected=selected,best_new_joint=best_new,gate_candidates=len(gate),
        selection_uses_test=False,selection_uses_coma=False,historical_test_exposure=True))
    predictions={k:L.evaluate(gals,cs,L.candidate_force(f),details=True) for k,f in selected.items()}
    save(out/"predictions.json",predictions)
    paired={k:C.bootstrap(p,refs["MOND"]) for k,p in predictions.items()};save(out/"paired-bootstrap.json",paired)
    print("Pooled selections: "+str({k:f["id"] for k,f in selected.items()}),flush=True)
    lights=[]
    for name,fit in selected.items():
        force=L.candidate_force(fit)
        for ne0,mstar in ((2.5e-3,.5e13),(4.5e-3,2e13)):
            src=D.coma_source(ne0,mstar)
            for reach in (3000.,9000.,30000.):
                row=L.optics(src,force,reach);row.update(model=name,ne0=ne0,mstar=mstar,primary=reach==9000.)
                lights.append(row)
    save(out/"lensing.json",dict(distance_Mpc=100.,data=D.coma_data(),models=lights,
        scope="Same fixed approximate static geometry and missing-data limitations as CMF-1; no free optical coupling."))
    fine=C.clusters(1200);checks=[];force=L.candidate_force(selected["joint"])
    for c,cf in zip(cs,fine):
        p0=D.pressure(c,force(c));p1=D.pressure(cf,force(cf))
        err=float(np.max(abs(p0-p1))/np.max(abs(p1)))
        checks.append(dict(test="pressure radial refinement",name=c["name"],relative_error=err,passed=err<.005))
    for row in [r for r in lights if r["model"] in ("joint","galaxy") and r["primary"]]:
        src=D.coma_source(row["ne0"],row["mstar"],8192);force=L.candidate_force(selected[row["model"]])
        refined=L.optics(src,force,9000,8192,.1,256,.001)
        err=float(np.max(abs(np.array(row["unit_beta_shear"])-refined["unit_beta_shear"]))/np.max(abs(refined["unit_beta_shear"])))
        checks.append(dict(test="source-memory-optical refinement",model=row["model"],ne0=row["ne0"],relative_error=err,passed=err<.005))
        for rmin in (.01,1.):
            result=L.optics(src,force,9000,8192,rmin,256,.001)
            change=float(np.max(abs(result["unit_beta_shear"]-refined["unit_beta_shear"]))/np.max(abs(refined["unit_beta_shear"])))
            checks.append(dict(test="inner memory boundary sensitivity",model=row["model"],ne0=row["ne0"],rmin=rmin,
                relative_change=change,bounded_chi2=result["bounded"]["chi2"],note="sensitivity, not a convergence assertion"))
    for name,fit in selected.items():
        force=L.candidate_force(fit)
        for g in gals:
            r=g["r"];newr=np.sort(np.r_[r,np.sqrt(r[:-1]*r[1:])])
            # Preserve the input zero exactly at its original row.
            gb=np.exp(np.interp(np.log(newr),np.log(r),np.log(np.maximum(g["gb"],1e-30/D.CODE_TO_SI))))
            gb[::2]=g["gb"]
            gf={**g,"r":newr,"gb":gb};before=np.sqrt(r*force(g));after=np.sqrt(newr*force(gf))[::2]
            checks.append(dict(test="galaxy profile grid sensitivity",model=name,name=g["name"],
                maximum_speed_change_kms=float(np.max(abs(after-before))),
                note="interpolated source-grid sensitivity; no new observed velocity used"))
    save(out/"numerical-checks.json",checks)
    p=predictions["joint"]
    save(out/"complete.json",dict(campaign_executed=True,fit_count=len(fits),attempts=sum(len(f["attempts"]) for f in fits),
        optimizer_success_fits=sum(f["success"] for f in fits),optimizer_success_attempts=sum(a["success"] for f in fits for a in f["attempts"]),
        numerical_passed=all(r.get("passed",True) for r in checks),
        point_score_beats_MOND=all(p[s]["galaxy_rmse"]<refs["MOND"][s]["galaxy_rmse"] for s in ("validation","test")),
        paired95_entirely_better=bool(paired["joint"]["percentile95"][1]<0),
        cluster_validation_under10=p["validation"]["cluster_chi2_per_point"]<=10,
        complete_solution=False,absolute_lensing_validated=False,derived_energy_funding=False))
    save(HERE/"gradient-evidence-sha256.json",{p.relative_to(HERE).as_posix():C.digest(p) for p in sorted(out.glob("*.json"))})
    print("Gradient campaign complete; failures retained.",flush=True)
if __name__=="__main__":main()

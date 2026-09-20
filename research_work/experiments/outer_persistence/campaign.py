"""Execute declared OP-1 and preserve all candidate results."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
import json,platform,subprocess,time
from datetime import datetime,timezone
import numpy as np
import scipy
import persistence as P
import run as C
L=P.L;D=P.D;HERE=P.HERE;save=L.save
def residual_map(gals,force):
    rows=[]
    for g in gals:
        v=np.sqrt(g["r"]*force(g))
        rows.extend(dict(name=g["name"],split=g["split"],r_kpc=float(r),fractional_radius=float(r/g["r"][-1]),
            gb_SI=float(b*D.CODE_TO_SI),mass=g["mass"],observed=float(y),predicted=float(p),
            signed_residual=float(p-y),error=float(e)) for r,b,y,p,e in zip(g["r"],g["gb"],g["y"],v,g["error"]))
    summaries=[]
    for split in ("train","validation","test"):
        for lower,upper in ((0,.33),(.33,.67),(.67,1.01)):
            rr=[r for r in rows if r["split"]==split and lower<=r["fractional_radius"]<upper]
            summaries.append(dict(split=split,fraction=[lower,upper],rows=len(rr),
                mean_signed_kms=np.mean([r["signed_residual"] for r in rr]),
                pooled_rmse_kms=np.sqrt(np.mean([r["signed_residual"]**2 for r in rr]))))
    return dict(rows=rows,bins=summaries,interpretation="Row-weighted radial diagnostics; selection uses equal-object RMSE.")
def safe_scores(gals,cs,fn,splits):
    for obj in gals+cs:
        if obj["split"] not in splits:continue
        a=fn(obj)
        if not np.all(np.isfinite(a)) or np.any(a<0):
            return None,dict(reason="nonfinite or negative force",object=obj["name"],minimum=float(np.nanmin(a)))
    return L.evaluate(gals,cs,fn,splits),None
def main():
    out=HERE/"evidence-v1";out.mkdir(exist_ok=False)
    files=D.FILES+list(HERE.glob("*.py"))+list(HERE.glob("*.md"))
    files+=list(L.HERE.glob("*.py"))+list((HERE.parent/"inverse_hitchhiking").glob("*.py"))
    files += [L.HERE/"summary.json",L.HERE/"evidence-v1"/"references.json"]
    save(out/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=P.ROOT,text=True).strip(),
        utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        hashes={p.relative_to(P.ROOT).as_posix():P.digest(p) for p in files}))
    gals=D.galaxies();cs=C.clusters()
    save(out/"inventory.json",dict(galaxies=[dict(name=g["name"],split=g["split"],n=len(g["r"]),distance_Mpc=g["catalog"]["distance"]) for g in gals],
        clusters=[dict(name=c["name"],split=c["split"],n=len(c["y"])) for c in cs],
        dark_matter=False,expansion=False,distance_fitting=False,historical_exposure=True))
    controls=P.controls(gals,cs);save(out/"controls.json",controls)
    if not all(c["passed"] for c in controls):raise AssertionError("Preserved control failure")
    print(f"{len(controls)} controls passed",flush=True)
    save(out/"baseline-residual-map.json",residual_map(gals,P.BASE_FORCE))
    objs=gals+cs;prepared={id(o):P.prepare(o) for o in objs}
    baseline={id(o):P.BASE_FORCE(o) for o in objs}
    base_scores=L.evaluate(gals,cs,P.BASE_FORCE,("train","validation"))
    cases=[dict(spec=P.BASELINE,scores=base_scores,valid=True)]
    current=None
    for i,spec in enumerate(P.specs(),1):
        key=tuple(spec[k] for k in ("delta","L0","eta"))
        if key!=current:
            corr={id(o):P.correction(prepared[id(o)],spec) for o in objs};current=key
        values={id(o):baseline[id(o)]+spec["A"]*corr[id(o)] for o in objs}
        scores,error=safe_scores(gals,cs,lambda o:values[id(o)],("train","validation"))
        cases.append(dict(spec=spec,scores=scores,valid=error is None,error=error))
        if i%25==0:
            save(out/"candidates.json",cases)
            print(f"{i}/225 candidate settings scored",flush=True)
    valid=[c for c in cases if c["valid"]]
    eligible=[c for c in valid if all(c["scores"][s]["cluster_chi2_per_point"]<=1.05*base_scores[s]["cluster_chi2_per_point"] for s in ("train","validation"))]
    protected=min(eligible,key=lambda c:(c["scores"]["validation"]["galaxy_rmse"],c["scores"]["train"]["galaxy_rmse"],c["spec"]["id"]))
    joint=min(valid,key=lambda c:(c["scores"]["validation"]["joint_objective"],c["spec"]["id"]))
    selected=dict(protected=protected["spec"],joint=joint["spec"])
    save(out/"selection-before-test.json",dict(selected=selected,eligible_pressure_protected=len(eligible),
        selection_uses_test=False,selection_uses_coma=False,historical_exposure=True))
    print("Frozen selections: "+str(selected),flush=True)
    refs=D.read(L.HERE/"evidence-v1"/"references.json")
    a0=refs["MOND_a0_SI"]
    mond=lambda o:.5*(o["gb"]+np.sqrt(o["gb"]**2+4*a0/D.CODE_TO_SI*o["gb"]))
    forces=dict(baseline=P.BASE_FORCE,MOND=mond,**{k:P.force(s) for k,s in selected.items()})
    predictions={};invalid=[]
    for name,fn in forces.items():
        sc,err=safe_scores(gals,cs,fn,("train","validation","test"))
        if err:invalid.append(dict(model=name,error=err))
        else:predictions[name]=L.evaluate(gals,cs,fn,details=True)
    save(out/"predictions.json",predictions);save(out/"invalid-final-selections.json",invalid)
    save(out/"paired-bootstrap.json",{k:{r:C.bootstrap(p,predictions[r]) for r in ("baseline","MOND")} for k,p in predictions.items() if k in selected})
    save(out/"selected-residual-maps.json",{k:residual_map(gals,forces[k]) for k in selected if k in predictions})
    lights=[];light_errors=[]
    for name in ("baseline","protected","joint"):
        if name not in predictions:continue
        for ne,stars in ((.0025,.5e13),(.0045,2e13)):
            src=D.coma_source(ne,stars)
            for reach in (3000.,9000.,30000.):
                def guarded(obj):
                    g=forces[name](obj)
                    if np.any(g<0) or not np.all(np.isfinite(g)):raise ValueError("invalid optical-domain force")
                    return g
                try:row=L.optics(src,guarded,reach)
                except ValueError as exc:
                    light_errors.append(dict(model=name,ne0=ne,reach=reach,error=str(exc)))
                    continue
                row.update(model=name,ne0=ne,mstar=stars,primary=reach==9000)
                lights.append(row)
    save(out/"lensing.json",dict(models=lights,failures=light_errors,data=D.coma_data(),distance_Mpc=100.,
        limitations="Exposed reconstructed bins, missing covariance/source distances; bounded beta is conditional, not absolute validation."))
    checks=[];sens=[]
    def check(test,name,value,limit=.005):checks.append(dict(test=test,name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
    fine=C.clusters(1200)
    for name,spec in selected.items():
        if name not in predictions:continue
        fn=forces[name];ref=P.force(spec,1024)
        for g in gals:
            before=fn(g);after=ref(g)
            check("galaxy internal grid refinement",name+" "+g["name"],np.max(abs(before-after))/max(np.max(abs(after)),1e-30))
            start=max(1,len(g["r"])//10)
            gg={**g,"r":g["r"][start:],"gb":g["gb"][start:]}
            changed=fn(gg)
            sens.append(dict(model=name,name=g["name"],removed_rows=start,
                max_speed_change_kms=float(np.max(abs(np.sqrt(gg["r"]*changed)-np.sqrt(gg["r"]*before[start:])))),
                interpretation="Missing interior-domain sensitivity; no refit."))
        for c,cf in zip(cs,fine):
            p0=D.pressure(c,fn(c));p1=D.pressure(cf,ref(cf))
            check("cluster source/internal grid refinement",name+" "+c["name"],np.max(abs(p0-p1))/np.max(abs(p1)))
        for row in [r for r in lights if r["model"]==name and r["primary"]]:
            source=D.coma_source(row["ne0"],row["mstar"],8192)
            high=L.optics(source,ref,9000,8192,.1,256,.001)
            err=np.max(abs(np.array(row["unit_beta_shear"])-high["unit_beta_shear"]))/max(abs(high["unit_beta_shear"]))
            check("combined optical refinement",name+" "+str(row["ne0"]),err)
    save(out/"numerical-checks.json",checks);save(out/"inner-domain-sensitivity.json",sens)
    results={}
    for name in selected:
        if name not in predictions:continue
        pp=predictions[name];bb=predictions["baseline"]
        improvement=all(pp[s]["galaxy_rmse"]<bb[s]["galaxy_rmse"] for s in ("validation","test"))
        pressure=all(pp[s]["cluster_chi2_per_point"]<=1.05*bb[s]["cluster_chi2_per_point"] for s in ("train","validation","test"))
        numerical=all(c["passed"] for c in checks if c["name"].startswith(name+" "))
        results[name]=dict(galaxy_improved=improvement,pressure_protected=pressure,numerical_passed=numerical,
            descriptive_gate_passed=improvement and pressure and numerical)
    save(out/"complete.json",dict(candidate_settings=225,baseline_cases=1,valid_train_validation_cases=len(valid),
        invalid_train_validation_cases=len(cases)-len(valid),results=results,
        complete_physical_solution=False,energy_funding_derived=False,vector_swirl_derived=False,
        microscopic_attachment_derived=False,absolute_lensing_validated=False))
    save(HERE/"evidence-sha256.json",{p.relative_to(HERE).as_posix():P.digest(p) for p in sorted(out.glob("*.json"))})
    print("Campaign complete: "+str(results),flush=True)
if __name__=="__main__":main()

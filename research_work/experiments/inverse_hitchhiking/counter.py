"""Declared sign-changing coherent response; first campaign remains intact."""
import hashlib,subprocess,time
from pathlib import Path
import numpy as np
from scipy.special import expit
from scipy.optimize import least_squares
import data as D
import model as M
from campaign import save,pack,scores
HERE=Path(__file__).resolve().parent

def fit_one(gp,ct,base,steepness,release_acceleration=None):
    Fg=M.features(gp["gb"],gp["r"],gp["mass"],"T")
    Fc=[M.features(c["gb"],c["r"],c["mass"],"T") for c in ct]
    norm=np.sqrt(sum(len(c["y"]) for c in ct)*10);cap=base["capacity"]
    def signed(gb,F,t):
        opposite=t[5]*gb*expit(steepness*(F[:,1]-t[6]))
        if release_acceleration is not None:
            opposite*=expit(2*(np.log(release_acceleration/1e-10)-F[:,1]))
        return M.extra(t[:5],cap,F)-opposite
    def residual(t):
        pred=np.sqrt(gp["r"]*(gp["gb"]+signed(gp["gb"],Fg,t)))
        chunks=[(pred-gp["y"])*gp["weight"]]
        if base["mode"]=="joint":
            for c,F in zip(ct,Fc):
                resid,_,_=D.pressure_residual(c,c["gb"]+signed(c["gb"],F,t))
                chunks.append(resid/norm)
        return np.concatenate(chunks)
    lo,hi=M.bounds("T");lo=np.r_[lo,0.,-8.];hi=np.r_[hi,.95,8.]
    attempts=[]
    for eps,flip in ((.1,0.),(.5,-2.),(.8,2.)):
        initial=np.r_[base["theta"],eps,flip];start=time.monotonic()
        out=least_squares(residual,initial,bounds=(lo,hi),xtol=1e-9,ftol=1e-9,gtol=1e-9,max_nfev=500)
        names=list(M.FAMILIES["T"])+["epsilon","x_flip"]
        attempts.append(dict(initial=initial,theta=out.x,success=out.success,status=out.status,message=out.message,
            nfev=out.nfev,optimality=out.optimality,objective=float(np.sum(out.fun**2)),seconds=time.monotonic()-start,
            bound_parameters=[n for n,v,a,b in zip(names,out.x,lo,hi) if min(v-a,b-v)<1e-4]))
    best=min(attempts,key=lambda a:a["objective"]);t=best["theta"]
    return dict(id=f"signed-{base['mode']}-cap{cap:.0e}-k{steepness:g}"+(f"-release{release_acceleration:.0e}" if release_acceleration is not None else ""),release_acceleration=release_acceleration,mode=base["mode"],family="T",
        capacity=cap,theta=t[:5],epsilon=t[5],x_flip=t[6],steepness=steepness,success=best["success"],
        parameters=dict(zip(M.FAMILIES["T"],t[:5])),bound_parameters=best["bound_parameters"],attempts=attempts)

def main(recovery=False):
    out=HERE/("evidence-recovery-v1" if recovery else "evidence-counter-v1");out.mkdir(exist_ok=False)
    sources=D.FILES+([HERE/"recovery.py",HERE/"recovery-protocol.md"] if recovery else [])+[HERE/p for p in ("data.py","model.py","campaign.py","counter.py","counter-protocol.md","protocol.md",
                                    "evidence-v1/fits.json")]
    save(out/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=D.ROOT,text=True).strip(),
        inputs_and_sources={p.relative_to(D.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}))
    bases=[f for f in D.read(HERE/"evidence-v1/fits.json") if f["family"]=="T"]
    gals=D.galaxies();clusters=D.cluster_inputs();gp=pack(gals,"train");ct=[c for c in clusters if c["split"]=="train"]
    fit=bases[0];zero=dict(fit,epsilon=0.,x_flip=0.,steepness=1.)
    expected=M.extra(fit["theta"],fit["capacity"],M.features(gp["gb"],gp["r"],gp["mass"],"T"))
    error=float(np.max(abs(M.response(zero,gp["gb"],gp["r"],gp["mass"])-expected)))
    lower=1-.95*expit(np.linspace(-100,100,1000))
    controls=[dict(name="epsilon zero restores original force",value=error,passed=error==0),
              dict(name="signed total acceleration lower bound",minimum=float(lower.min()),passed=bool(lower.min()>=.05-1e-14))]
    if not all(c["passed"] for c in controls):raise AssertionError(controls)
    save(out/"controls.json",controls)
    fits=[]
    cases=([(b,2.,a) for b in bases if b["capacity"]==3e-9 for a in (1e-9,1e-8,1e-7)] if recovery else
           [(b,k,None) for b in bases for k in (.5,1.,2.)])
    for base,k,release in cases:
        fit=fit_one(gp,ct,base,k,release)
        fit["scores"]=scores(gals,clusters,fit,("train","validation"))
        fits.append(fit);save(out/"fits.json",fits)
        v=fit["scores"]["validation"]
        print(f"{fit['id']}: galaxy {v['galaxy']['rmse']:.3f}, cluster {v['cluster']['chi2_per_point']:.3f}; success {fit['success']}",flush=True)
    selected={}
    for mode in ("galaxy","joint"):
        valid=[f for f in fits if f["mode"]==mode and f["success"]]
        selected[mode]=min(valid,key=lambda f:(f["scores"]["validation"]["joint_objective"] if mode=="joint" else
                                            f["scores"]["validation"]["galaxy"]["rmse"],f["id"]))
    save(out/"selection-before-test.json",{m:{k:f[k] for k in ("id","family","capacity","theta","epsilon","x_flip","steepness","release_acceleration")}
                                          for m,f in selected.items()})
    save(out/"selected-predictions.json",{m:scores(gals,clusters,f,("train","validation","test"),True) for m,f in selected.items()})
    fit=selected["joint"]
    if recovery:
        gb=np.array([1e4/D.CODE_TO_SI]);rr=np.array([10.])
        ratio=(gb+M.response(fit,gb,rr,1e11))/gb
        bound=fit["capacity"]/1e4+fit["epsilon"]*(fit["release_acceleration"]/1e4)**2
        controls.append(dict(name="large acceleration ordinary-force limit",error=float(abs(ratio[0]-1)),
            analytic_bound=bound,passed=bool(abs(ratio[0]-1)<1e-10 and abs(ratio[0]-1)<=bound+1e-15)))
        save(out/"controls.json",controls)
    coma=D.coma_data();light=[]
    for ne0,mstar in ((2.5e-3,.5e13),(4.5e-3,2e13)):
        source=D.coma_source(ne0,mstar)
        for reach in (3000.,9000.,30000.):
            for eta in (0.,.25,1.,4.):
                g=M.optical_force(source,fit,reach,eta)
                result=M.shear_shape(g,coma["r"],128,.002,(3000.,reach) if reach!=3000 else (reach,))
                light.append(dict(ne0=ne0,mstar=mstar,reach_kpc=reach,eta=eta,primary=reach==9000 and eta==1,
                    alpha=result["alpha"],shear_shape=result["shape"],fit=M.shape_fit(result["shape"],coma["y"],coma["error"])))
    save(out/"coma.json",dict(data=coma,models=light,selected_law=fit["id"],status="Shape only, exposed data; one normalization per profile."))
    numerical=[]
    for c,fine in zip(clusters,D.cluster_inputs(1200)):
        p=D.pressure(c,c["gb"]+M.response(fit,c["gb"],c["r"],c["mass"]))
        q=D.pressure(fine,fine["gb"]+M.response(fit,fine["gb"],fine["r"],fine["mass"]))
        err=float(np.max(abs(p-q))/np.max(abs(q)))
        numerical.append(dict(test="pressure grid",name=c["name"],relative_max_error=err,passed=err<.005))
    for row in [r for r in light if r["primary"]]:
        source=D.coma_source(row["ne0"],row["mstar"]);g=M.optical_force(source,fit,9000,1.)
        q=M.shear_shape(g,coma["r"],256,.001,(3000.,9000.))
        err=float(np.max(abs(np.array(row["shear_shape"])-q["shape"]))/np.max(abs(q["shape"])))
        numerical.append(dict(test="light projection",ne0=row["ne0"],relative_max_error=err,passed=err<.005))
    save(out/"numerical-checks.json",numerical)
    save(out/"complete.json",dict(complete=True,fits=len(fits),attempts=sum(len(f["attempts"]) for f in fits),
        numerical_passed=all(r["passed"] for r in numerical),
        operational_validation_target=fit["scores"]["validation"]["galaxy"]["rmse"]<=22 and fit["scores"]["validation"]["cluster"]["chi2_per_point"]<=10))
    save(HERE/("recovery-evidence-sha256.json" if recovery else "counter-evidence-sha256.json"),{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for p in out.rglob("*") if p.is_file()})
    print("Completed "+str({m:f["id"] for m,f in selected.items()}),flush=True)
if __name__=="__main__":main()

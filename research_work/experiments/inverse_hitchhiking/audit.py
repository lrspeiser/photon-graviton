"""Independent record, source and observational-prediction verification."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
from scipy.integrate import quad
import data as D
import model as M
from campaign import save
HERE=Path(__file__).resolve().parent
RUNS=("evidence-v1","evidence-counter-v1","evidence-recovery-v1")
CHECKSUMS=("evidence-sha256.json","counter-evidence-sha256.json","recovery-evidence-sha256.json")
def main():
    checks=[];hashes=0;sources=0
    def add(name,passed,detail=None):checks.append(dict(name=name,passed=bool(passed),detail=detail))
    for file in CHECKSUMS:
        bad=[]
        for rel,digest in D.read(HERE/file).items():
            hashes+=1
            if hashlib.sha256((HERE/rel).read_bytes()).hexdigest()!=digest:bad.append(rel)
        add(file+" byte digests",not bad,bad)
    for run in RUNS:
        manifest=D.read(HERE/run/"manifest.json");bad=[]
        for rel,digest in manifest["inputs_and_sources"].items():
            sources+=1
            raw=subprocess.check_output(["git","show",manifest["git_head"]+":"+rel],cwd=D.ROOT,stderr=subprocess.PIPE)
            if hashlib.sha256(raw).hexdigest()!=digest:bad.append(rel)
        add(run+" input/source commit",not bad,bad)
        fits=D.read(HERE/run/"fits.json");selection=D.read(HERE/run/"selection-before-test.json")
        for mode in ("galaxy","joint"):
            candidates=[f for f in fits if f["mode"]==mode and f["success"]]
            key=lambda f:(f["scores"]["validation"]["joint_objective"] if mode=="joint" else
                          f["scores"]["validation"]["galaxy"]["rmse"],len(f["theta"]),f["id"])
            winner=min(candidates,key=key)
            add(run+" "+mode+" validation-only selection",winner["id"]==selection[mode]["id"])
        add(run+" controls",all(c["passed"] for c in D.read(HERE/run/"controls.json")))
        add(run+" quadrature refinements",all(c["passed"] for c in D.read(HERE/run/"numerical-checks.json")))
    gals={g["name"]:g for g in D.galaxies()}
    final=D.read(HERE/"evidence-recovery-v1/selection-before-test.json")["joint"]
    pred=D.read(HERE/"evidence-recovery-v1/selected-predictions.json")["joint"]
    maxima=0.;negative=0;negative3=0;zero=0
    for split,s in pred.items():
        expected_rmse=[];chi=0;n=0
        for row in s["galaxy"]["rows"]:
            g=gals[row["name"]];x=np.log(np.maximum(g["gb"]*D.CODE_TO_SI,1e-30)/1e-10)
            y=np.log(g["r"]/10);z=np.log(g["mass"]/1e11);w=np.tanh(np.log(g["mass"]/1e12)/2)
            c,q,ss,m,d=final["theta"];u=c+q*x+ss*y+m*z+d*w
            fraction=1/(1+np.exp(-u))
            turn=1/(1+np.exp(-final["steepness"]*(x-final["x_flip"])))
            release=1/(1+(g["gb"]*D.CODE_TO_SI/final["release_acceleration"])**2)
            total=g["gb"]*(1-final["epsilon"]*turn*release)+final["capacity"]/D.CODE_TO_SI*fraction
            velocity=np.sqrt(g["r"]*total)
            maxima=max(maxima,float(np.max(abs(velocity-row["predicted_kms"]))))
            expected_rmse.append(np.mean((velocity-g["y"])**2))
            chi+=float(np.sum(((velocity-g["y"])/g["error"])**2));n+=len(velocity)
            vb=np.sqrt(g["r"]*g["gb"]);negative+=int(np.sum(g["y"]<vb))
            negative3+=int(np.sum(g["y"]+3*g["error"]<vb));zero+=int(np.sum(g["gb"]==0))
            add("fixed distance "+g["name"],row["distance_Mpc"]==g["catalog"]["distance"])
        add(split+" galaxy aggregate RMSE",abs(np.sqrt(np.mean(expected_rmse))-s["galaxy"]["rmse"])<1e-9)
        add(split+" galaxy aggregate chi2",abs(chi/n-s["galaxy"]["chi2_per_point"])<1e-8)
    add("independent final formula versus saved velocities",maxima<1e-8,maxima)
    coma=D.read(HERE/"evidence-recovery-v1/coma.json");pr=[r for r in coma["models"] if r["primary"]][0]
    source=D.coma_source(pr["ne0"],pr["mstar"]);force=M.optical_force(source,final,9000,1)
    independent=[]
    for b,target in zip(coma["data"]["r"],pr["alpha"]):
        points=[float(np.arccos(b/r)) for r in (3000,9000) if r>b]
        value=4/D.C**2*quad(lambda t:float(force(np.array([b/np.cos(t)]))[0])*b/np.cos(t),
                           0,np.pi/2,points=points,epsabs=1e-5,epsrel=1e-8,limit=500)[0]
        independent.append(abs(value/target-1))
    add("adaptive versus Gauss projected deflection",max(independent)<.005,max(independent))
    # All input target measurements and adopted distances remain identical to their manifests.
    fitcounts=[len(D.read(HERE/r/"fits.json")) for r in RUNS]
    attempts=[a for r in RUNS for f in D.read(HERE/r/"fits.json") for a in f["attempts"]]
    result=dict(all_passed=all(c["passed"] for c in checks),hashes_checked=hashes,pinned_files_checked=sources,
        fits=sum(fitcounts),attempts=len(attempts),optimizer_failures=sum(not a["success"] for a in attempts),
        ordinary_source_zero_acceleration_rows=zero,inverse_negative_rows=negative,
        inverse_negative_beyond_three_velocity_errors=negative3,
        note="Velocity-only significance omits ordinary-matter model uncertainties. Numerical checks are not physical validation.",
        checks=checks)
    save(HERE/"audit.json",result)
    print(json.dumps({k:v for k,v in result.items() if k!="checks"}))
    if not result["all_passed"]:raise SystemExit("Audit failed")
if __name__=="__main__":main()

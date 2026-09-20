"""Numerical diagnostics declared in numerical-followup.md."""
import hashlib,json,subprocess
from pathlib import Path
import numpy as np
import hitch as h
import curvature as c
from campaign import dump
HERE=Path(__file__).resolve().parent

def compare(a,b):
    aa={(r["id"],r["impact"]):r for r in a["results"]}
    result=[]
    for r in b["results"]:
        q=aa[(r["id"],r["impact"])]
        arrived=q["arrived"] and r["arrived"]
        da=abs(q["angle"]-r["angle"]) if arrived else None
        dy=abs(q["observer_y"]-r["observer_y"]) if arrived else None
        result.append(dict(id=r["id"],impact=r["impact"],both_arrived=arrived,
            both_censored=not q["arrived"] and not r["arrived"],
            arrival_mismatch=q["arrived"]!=r["arrived"],
            angle_difference=da,position_difference=dy,
            observer_pass=bool(da<.01 and dy<.01) if arrived else None))
    return result

def main():
    out=HERE/"evidence-numerical-v1";out.mkdir(exist_ok=False)
    dump(out/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=h.REPO,text=True).strip(),
        sources={p.relative_to(h.REPO).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in [HERE/"numerical.py",HERE/"numerical-followup.md",HERE/"hitch.py",HERE/"curvature.py",h.STREAM/"campaign.py"]}))
    checks=[]
    def check(name,value,limit):
        checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
    # Exact autonomous single-bound heading equation tan(delta/2)=tan(delta0/2)*exp(-k*t).
    y=np.array([[0.,0.,.3]]);heading=np.array([1.]);k=np.array([2.])
    for _ in range(100):y=h.fly(y,.01,np.array([True]),heading,k)
    exact=1+2*np.arctan(np.tan((.3-1)/2)*np.exp(-2))
    check("bound steering versus exact solution",abs(y[0,2]-exact),1e-8)
    p={k:np.array([v]) for k,v in dict(pitch=2.,on=100.,off=.7,kappa=4.).items()}
    for mode,rhs in [("heading",h.mrhs),("curvature",c.rhs)]:
        y=np.zeros((1,9));y[0,:2]=[100.,100.];y[0,3]=.8;y[0,4:6]=[.3,.4]
        for step in range(100):y=h.rk(step*.01,y,.01,p,rhs)
        check(mode+" detached-front decay",np.max(np.abs(y[0,3:6]-np.array([.8,.3,.4])*np.exp(-.7))),1e-8)
    ledger=[]
    for root in ["evidence-v1","evidence-curvature-v1"]:
        for run in ["events-dt001","events-dt0005"]:
            f=HERE/root/run/"trajectories.npz"
            with np.load(f) as z:
                residual=int(np.max(np.abs(z["captures"]-z["releases"]-z["bound"][:,-1])))
                check(root+"/"+run+" event balance",residual,0)
                check(root+"/"+run+" nonnegative events",max(0,-int(z["captures"].min()),-int(z["releases"].min())),0)
                arrived=np.isfinite(z["exit_time"])
                check(root+"/"+run+" observer position",np.max(np.abs(z["state"][arrived,-1,0]-20)),1e-10)
                dist=np.linalg.norm(np.diff(z["state"][:,:,:2],axis=1),axis=2)
                excess=np.max(dist-np.diff(z["time"])[None,:])
                check(root+"/"+run+" unit speed chord bound",max(0,excess),1e-9)
                ledger.append(dict(archive=str(f.relative_to(HERE)),histories=len(arrived),arrivals=int(arrived.sum())))
    suites={}
    for mode,rows,rhs,bound in [("heading",[h.REPS[1]],h.mrhs,1.),
                               ("curvature",h.REFINE+[h.row(2,100,.1,1)],c.rhs,2.)]:
        lo=h.mean_trace(rows,out/(mode+"-dt0005"),.005,fun=rhs,moment_bound=bound)
        hi=h.mean_trace(rows,out/(mode+"-dt00025"),.0025,fun=rhs,moment_bound=bound)
        diffs=compare(lo,hi)
        with np.load(out/(mode+"-dt0005")/"trajectories.npz") as a,np.load(out/(mode+"-dt00025")/"trajectories.npz") as b:
            censored=~np.isfinite(a["exit_time"])&~np.isfinite(b["exit_time"])
            terminal=np.max(np.abs(a["state"][censored,-1]-b["state"][censored,-1]),axis=0).tolist() if censored.any() else []
        suites[mode]=dict(comparisons=diffs,arrived_pass=sum(x["observer_pass"] is True for x in diffs),
            both_arrived=sum(x["both_arrived"] for x in diffs),both_censored=sum(x["both_censored"] for x in diffs),
            arrival_mismatches=sum(x["arrival_mismatch"] for x in diffs),censored_terminal_max_by_component=terminal,
            occupations_pass=lo["occupation_pass"] and hi["occupation_pass"])
        if mode=="curvature":
            reflected=[h.row(-r["pitch"],r["on"],r["off"],r["kappa"]) for r in rows]
            mirror=h.mean_trace(reflected,out/"curvature-mirror-dt0005",.005,fun=rhs,moment_bound=bound)
            with np.load(out/"curvature-dt0005"/"trajectories.npz") as a,np.load(out/"curvature-mirror-dt0005"/"trajectories.npz") as b:
                v=b["state"].reshape(len(rows),10,-1,9)[:,::-1].reshape(a["state"].shape).copy()
                v[:,:,1]*=-1;v[:,:,2]*=-1;v[:,:,5]*=-1
                check("HH2 full reflected state",np.max(np.abs(v-a["state"])),1e-7)
    result=dict(controls=checks,controls_passed=all(x["passed"] for x in checks),event_ledger=ledger,suites=suites)
    dump(out/"summary.json",result)
    dump(out/"complete.json",dict(complete=True))
    dump(HERE/"numerical-evidence-sha256.json",{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
         for p in out.rglob("*") if p.is_file()})
    print(json.dumps({k:{q:v for q,v in a.items() if q not in ("comparisons","censored_terminal_max_by_component")} for k,a in suites.items()}),flush=True)
    print(f"controls pass: {result['controls_passed']}",flush=True)
if __name__=="__main__":main()

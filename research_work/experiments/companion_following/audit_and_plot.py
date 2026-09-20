"""Audit preserved GF-1 artifacts, compare refinements, and render actual trajectories."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from campaign import HERE, REPO, write_json

EVIDENCE = HERE/"evidence"
def read(path): return json.loads(path.read_text(encoding="utf8"))
def mapping(rows): return {r["id"]: r for r in rows}

def main():
    screen = read(EVIDENCE/"screen-completion-v1"/"screen.json")
    refine = read(EVIDENCE/"refine-v1"/"refinement.json")
    selection = read(EVIDENCE/"refine-v1"/"selection.json")
    controls = read(EVIDENCE/"controls-v1"/"controls.json")
    audit = dict(controls_passed=controls["passed"], controls_count=controls["count"],
                 cases={}, family_counts={}, refinement=[], symmetry=[], trajectory_count=0,
                 nonfinite_trajectories=0, checks=[], inputs_are_exposed=True,
                 observational_fit_performed=False)
    def check(name, passed, detail=None):
        audit["checks"].append(dict(name=name, passed=bool(passed), detail=detail))
    for group in screen:
        a=mapping(screen[group]["chain"]); b=mapping(screen[group]["ring"])
        joint=[i for i in a if a[i]["behavioral_pass"] and b[i]["behavioral_pass"]]
        audit["cases"][group] = dict(count=len(a),
            chain_pass=sum(r["behavioral_pass"] for r in a.values()),
            ring_pass=sum(r["behavioral_pass"] for r in b.values()),
            joint_behavior_pass=len(joint), joint_behavior_ids=joint,
            conservative_ring_pass=sum(r["conservation_pass"] for r in b.values()),
            unavailable=sum(r.get("unavailable",False) for r in b.values()))
        for family in sorted(set(i.split("-")[0] for i in a)):
            ids=[i for i in a if i.split("-")[0]==family]
            audit["family_counts"][family] = dict(count=len(ids),
                chain=sum(a[i]["behavioral_pass"] for i in ids),
                ring=sum(b[i]["behavioral_pass"] for i in ids),
                both=sum(a[i]["behavioral_pass"] and b[i]["behavioral_pass"] for i in ids))
        for fixture in ("chain","ring"):
            fine=mapping(refine[f"{group}-{fixture}-dt002"])
            coarse=mapping(screen[group][fixture])
            metric="chain_gain" if fixture=="chain" else "radius_ratio_final"
            for identifier in selection[group]:
                aa,bb=coarse[identifier],fine[identifier]
                change=abs(aa[metric]-bb[metric])/max(1,abs(bb[metric]))
                audit["refinement"].append(dict(id=identifier,fixture=fixture,
                    metric=metric,normalized_change=change,converged=change<.05,
                    coarse_behavioral_pass=aa["behavioral_pass"],
                    fine_behavioral_pass=bb["behavioral_pass"],
                    fine_conservation_pass=bb["conservation_pass"]))
    for identifier in selection["phenomenological"]:
        # Compare complete reflected state histories, with omega a pseudoscalar.
        coarse=np.load(EVIDENCE/"screen-v1"/"phenomenological-chain"/"trajectories.npz")
        rev=np.load(EVIDENCE/"refine-v1"/"reverse-turn"/"trajectories.npz")
        ia=list(coarse["ids"]).index(identifier); ib=list(rev["ids"]).index(identifier)
        expected=coarse["state"][ia].copy()
        expected[..., [1,3,5,7,8]] *= -1
        error=float(np.max(np.abs(expected-rev["state"][ib])))
        audit["symmetry"].append(dict(id=identifier,reflection_state_error=error))
    check("all selected mirror trajectories reflect",all(r["reflection_state_error"]<1e-9 for r in audit["symmetry"]))
    for name in ("zero-input","conservative-zero-input"):
        check(name+" stays straight",all(r["no_turn_transverse_max"]<1e-12 for r in refine[name]))
    # Check source hashes against the exact version named in each run manifest.
    source_checks=[]
    for path in sorted(EVIDENCE.rglob("manifest.json")):
        data=read(path)
        for filename,digest in data["source_sha256"].items():
            rel=(HERE/filename).relative_to(REPO).as_posix()
            committed=subprocess.check_output(["git","show",data["git_head"]+":"+rel],cwd=REPO)
            source_checks.append(hashlib.sha256(committed).hexdigest()==digest)
    check("all archived source hashes match their recorded commits",all(source_checks),len(source_checks))

    max_leader_error=0.
    for summary_path in sorted(EVIDENCE.rglob("summary.json")):
        data=read(summary_path)
        count=len(data["results"])
        audit["trajectory_count"]+=count
        audit["nonfinite_trajectories"]+=sum(r["first_nonfinite_time"] is not None for r in data["results"])
        raw=np.load(summary_path.with_name("trajectories.npz"))
        check(str(summary_path.parent.relative_to(EVIDENCE))+" IDs match",
              list(raw["ids"])==[r["id"] for r in data["results"]])
        if data["fixture"]=="chain":
            max_leader_error=max(max_leader_error,max(r["leader_speed_error"] for r in data["results"]))
        if data["fixture"]=="ring" and str(raw["ids"][0]).startswith("C-"):
            # Independent direct pair-sum H and canonical L, without importing the model.
            x=raw["state"][...,:2]; v=raw["state"][...,2:4]
            dx=x[:,:,:,None,:]-x[:,:,None,:,:]
            dv=v[:,:,:,None,:]-v[:,:,None,:,:]
            rr=np.linalg.norm(dx,axis=-1); s=rr/2.5
            nn=x.shape[2]; mask=1-np.eye(nn)
            U=np.sum((2*np.exp(-rr/.4)-np.exp(-rr/2))*mask,axis=(-2,-1))/32
            H=[]; L=[]
            for i,identifier in enumerate(raw["ids"]):
                kid=int(str(identifier).split("-")[1][1:])
                eta=float(str(identifier).split("-")[2][1:])
                ss=s[i]
                if kid==0: K=np.exp(-ss**2)
                elif kid==1: K=np.exp(-ss)
                elif kid==2: K=(1+ss**2)**-1
                elif kid==3: K=(1+ss**2)**-2
                elif kid==4: K=np.maximum(1-ss,0)**4*(1+4*ss)
                else: K=ss**2*np.exp(1-ss**2)
                K*=mask
                momentum=v[i]+eta/16*np.sum(K[...,None]*dv[i],axis=2)
                H.append(.5*np.sum(v[i]*v[i],axis=(-2,-1))
                     +eta/64*np.sum(K*np.sum(dv[i]**2,axis=-1),axis=(-2,-1))+U[i])
                L.append(np.sum(x[i,...,0]*momentum[...,1]-x[i,...,1]*momentum[...,0],axis=-1))
            check(str(summary_path.parent.relative_to(EVIDENCE))+" independent energy",
                  np.max(np.abs(np.array(H)-raw["energy"]))<1e-10)
            check(str(summary_path.parent.relative_to(EVIDENCE))+" independent canonical angular momentum",
                  np.max(np.abs(np.array(L)-raw["angular"]))<1e-10)
    audit["max_prescribed_leader_speed_error"]=max_leader_error
    check("prescribed leader speed integration",max_leader_error<1e-7)
    audit["all_evidence_checks_passed"]=all(c["passed"] for c in audit["checks"])
    write_json(HERE/"audit.json",audit)

    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    fig,axs=plt.subplots(2,2,figsize=(13,9),layout="constrained")
    families=list(audit["family_counts"])
    ind=np.arange(len(families))
    for j,(metric,label,color) in enumerate([("chain","Turn chain","#247BA0"),("ring","Prepared ring","#E8AA42"),("both","Both","#8C5AA8")]):
        axs[0,0].bar(ind+(j-1)*.25,[audit["family_counts"][f][metric] for f in families],width=.25,label=label,color=color)
    axs[0,0].set(xticks=ind,xticklabels=families,ylabel="Passing cases",title="Declared behavioral tests (600 + 18 laws)")
    axs[0,0].legend(frameon=False)
    best=selection["phenomenological"][0]
    ar=np.load(EVIDENCE/"screen-v1"/"phenomenological-chain"/"trajectories.npz")
    ix=list(ar["ids"]).index(best)
    angles=np.arctan2(ar["state"][ix,...,3],ar["state"][ix,...,2])
    for i in (0,1,3,7,11):
        axs[0,1].plot(ar["time"],angles[:,i],label=f"packet {i}")
    axs[0,1].set(xlabel="Time (declared units)",ylabel="Velocity angle (rad)",title=f"Selected exposed chain: {best}")
    axs[0,1].legend(frameon=False,ncol=2)
    long=np.load(EVIDENCE/"refine-v1"/"conservative-long-ring"/"trajectories.npz")
    # Plot first selected conservative law, including any late failure.
    bestc=str(long["ids"][0])
    for i in range(long["state"].shape[2]):
        axs[1,0].plot(long["state"][0,:,i,0],long["state"][0,:,i,1],lw=.7,alpha=.7)
    axs[1,0].set(xlabel="x",ylabel="y",title=f"Extended prepared ring: {bestc}, t = 256")
    axs[1,0].set_aspect("equal",adjustable="datalim")
    for i,identifier in enumerate(long["ids"]):
        axs[1,1].plot(long["time"],long["radius"][i]/long["radius"][i,0],label=str(identifier))
    axs[1,1].axhline(.7,color="black",ls="--",lw=.8)
    axs[1,1].axhline(1.3,color="black",ls="--",lw=.8)
    axs[1,1].set(xlabel="Time",ylabel="Mean radius / initial",title="Extended conservative rings; dashed = bounds")
    axs[1,1].legend(frameon=False)
    fig.suptitle("GF-1: companion following is a mechanics screen, not a lensing solution",fontsize=14)
    fig.savefig(HERE/"results.png",dpi=170)
    plt.close(fig)
    checksums={}
    for p in sorted(EVIDENCE.rglob("*")):
        if p.is_file(): checksums[p.relative_to(HERE).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    write_json(HERE/"evidence-sha256.json",checksums)
    print(json.dumps({k:v for k,v in audit.items() if k in ("cases","family_counts","trajectory_count",
        "nonfinite_trajectories","all_evidence_checks_passed","max_prescribed_leader_speed_error")},indent=2))
    if not audit["all_evidence_checks_passed"]: raise SystemExit("Evidence audit failed; details preserved")

if __name__=="__main__": main()

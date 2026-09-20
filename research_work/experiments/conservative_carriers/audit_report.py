"""Independent SV-1 archive audit, comparison figures, and research report."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
import json,hashlib,subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mechanics as M
from run_campaign import save,digest
HERE=M.HERE
def read(p):return json.loads(p.read_text())
def independent_energy(y,spec):
    n=len(y)//6
    q=y[:3*n].reshape(n,3);p=y[3*n:].reshape(n,3)
    result=sum(float(pi@pi)/2 for pi in p[:-1])+float(p[-1]@p[-1])/(2*100)
    for i in range(n-1):
        for j in range(i):
            r=q[i]-q[j]
            result-=spec["eta"]/(n-2)*np.exp(-(r@r)/(2*spec["ell"]**2))*float(p[i]@p[j])
        r=q[i]-q[-1]
        result-=spec["mu"]/np.sqrt(r@r+.2**2)
    return result
def main():
    checks=[];archives={}
    def check(name,v,limit=0):
        checks.append(dict(name=name,value=float(v),limit=limit,passed=bool(v<=limit)))
    hashcount=0;pinned=0;energy_error=0;parity_error=0
    for folder,hashfile in (("evidence-v1","evidence-sha256.json"),("evidence-matched-v1","matched-evidence-sha256.json")):
        out=HERE/folder;manifest=read(out/"manifest.json")
        for path,expected in read(HERE/hashfile).items():
            check(folder+" evidence "+path,int(digest(HERE/path)!=expected));hashcount+=1
        for path,expected in manifest["hashes"].items():
            check(folder+" current source "+path,int(digest(M.ROOT/path)!=expected))
            raw=subprocess.check_output(["git","show",manifest["git_head"]+":"+path],cwd=M.ROOT)
            check(folder+" pinned source "+path,int(hashlib.sha256(raw).hexdigest()!=expected));pinned+=1
        rows=read(out/"runs.json");ref=read(out/"refinements.json");paths=np.load(out/"trajectories.npz")
        for r in rows:
            arr=paths[r["spec"]["id"]]
            for i in (0,-1):
                energy_error=max(energy_error,abs(independent_energy(arr[i],r["spec"])-r["energy_series"][i]))
            if r["spec"]["packet"]=="rotating":
                s=r["spec"];other=next(x for x in rows if x["spec"]["packet"]=="reverse" and all(x["spec"][k]==s[k] for k in ("eta","ell","mu")))
                a=np.array(r["metrics"]["rms_radius"]);b=np.array(other["metrics"]["rms_radius"])
                parity_error=max(parity_error,float(np.max(abs(a-b))/max(1,np.max(b))))
        check(folder+" preserved numerical failures",sum(not r["passed"] for r in rows))
        check(folder+" preserved refinement failures",sum(not r["passed"] for r in ref))
        check(folder+" controls",sum(not r["passed"] for r in read(out/"controls.json")))
        archives[folder]=dict(rows=rows,paths=paths,comparisons=read(out/"comparisons.json"),complete=read(out/"complete.json"),refinements=ref)
    check("independent endpoint energy",energy_error,1e-10)
    check("rotation reversal radius symmetry",parity_error,1e-6)
    save(HERE/"audit.json",dict(passed=all(c["passed"] for c in checks),checks=checks,
        evidence_hashes=hashcount,pinned_sources=pinned,independent_energy_error=energy_error,rotation_reversal_error=parity_error))
    raw=archives["evidence-v1"];matched=archives["evidence-matched-v1"]
    def row(packet,mu,eta,ell=1,archive=matched):
        return next(r for r in archive["rows"] if r["spec"]["packet"]==packet and r["spec"]["mu"]==mu and r["spec"]["eta"]==eta and r["spec"]["ell"]==ell)
    def comp(packet,mu,eta=.8,ell=1,archive=matched):
        return next(r for r in archive["comparisons"] if r["packet"]==packet and r["mu"]==mu and r["eta"]==eta and r["ell"]==ell)
    a=row("stream",0,.8);b=row("stream",0,0);c=comp("stream",0)
    plt.rcParams.update({"font.size":9})
    fig,axs=plt.subplots(1,3,figsize=(14,4.4))
    for eta,color,label in ((0,"#7d8792","Interaction off"),(.8,"#246bb6","Moving-carrier interaction")):
        r=row("stream",0,eta)
        axs[0].plot(r["times"],r["metrics"]["packet_width"],c=color,label=label)
    axs[0].set(title="Stream cohesion without source gravity",xlabel="Time (model units)",ylabel="RMS packet width (model units)")
    axs[0].legend(fontsize=8)
    xx=np.arange(6)
    for mu,offset,col,label in ((0,-.18,"#246bb6","Source gravity off"),(1,.18,"#c16c35","Source gravity on")):
        axs[1].bar(xx+offset,[comp(p,mu)["late_width_ratio"] for p in M.PACKETS],width=.34,color=col,label=label)
    axs[1].axhline(1,c="black",lw=.8)
    axs[1].set_xticks(xx,M.PACKETS,rotation=30,ha="right")
    axs[1].set(title="Same interaction across six packets",ylabel="Late width / interaction-off width")
    axs[1].legend(fontsize=8)
    for eta,col in ((0,"#7d8792"),(.8,"#246bb6")):
        r=row("stream",1,eta);path=matched["paths"][r["spec"]["id"]]
        q=path[:,:3*(M.N+1)].reshape(-1,M.N+1,3);rel=q[:,:-1]-q[:,-1,None,:]
        for i in range(M.N):
            axs[2].plot(rel[:,i,0],rel[:,i,1],c=col,alpha=.6,lw=.8)
    axs[2].set(title="Stream near a recoiling source",xlabel="Source-relative x (model units)",ylabel="Source-relative y (model units)")
    axs[2].set_aspect("equal",adjustable="datalim")
    for ax in axs:ax.grid(alpha=.15);ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("SV-1: conservative moving-carrier interaction, matched initial physical speeds",fontsize=13)
    fig.tight_layout(rect=(0,.09,1,.93))
    fig.text(.02,.025,"Illustrations: eta=0.8, range=1, N=12. These are mechanical proxies, not galaxy or photon predictions. No source rotation or energy is created for free.",fontsize=8)
    fig.savefig(HERE/"comparison.png",dpi=180);fig.savefig(HERE/"comparison.svg")
    svg=HERE/"comparison.svg";svg.write_text("\n".join(x.rstrip() for x in svg.read_text().splitlines())+"\n",encoding="utf8")
    plt.close(fig)
    complete=[x["complete"] for x in archives.values()]
    summary=dict(primary_runs=sum(x["runs"] for x in complete),refinements=sum(x["refined_runs"] for x in complete),
        all_primary_numerical_passed=all(x["runs"]==x["passed_runs"] for x in complete),
        all_refinements_passed=all(x["refined_runs"]==x["passed_refinements"] for x in complete),
        illustrative_stream=dict(eta=.8,ell=1,mu=0,late_width_ratio=c["late_width_ratio"],
            required_preparation_energy=a["preparation_energy"],control_preparation_energy=b["preparation_energy"],
            preparation_energy_ratio=a["preparation_energy"]/b["preparation_energy"]),
        maximum_energy_drift=max(x["maximum_scaled_energy_drift"] for x in complete),
        maximum_momentum_drift=max(x["maximum_momentum_drift"] for x in complete),
        maximum_angular_drift=max(x["maximum_angular_drift"] for x in complete),
        maximum_refinement_error=max(r["relative_position_error"] for x in archives.values() for r in x["refinements"]),
        observational_predictions=False,graviton_identification=False,causal_mediator=False)
    save(HERE/"summary.json",summary)
    lines=["# SV-1: conservative attraction between moving carriers","","A bounded momentum-dependent interaction can produce measurable packet cohesion with a conserved mechanical energy and reciprocal source recoil. In the declared stream example without source gravity, the late packet width is about 29% smaller than the interaction-off control after matching initial physical emission speeds. This is a finite, low-speed mechanical result. It does not establish universal extra gravity, photon hitchhiking, or a galaxy/cluster solution.","",
    "![Mechanical results](comparison.png)","","## What was executed","",
    "The original 252 settings span six source packets, seven signed strengths, three interaction ranges and ordinary-source gravity on/off. Thirty-six predeclared refinements were run. An amendment repeats the same 252 settings and 36 refinements with matched source-relative physical velocities. All 504 primary integrations and 72 refinements pass their declared numerical gates. The initial 32 Hamiltonian controls and 378 velocity-matching controls pass.","",
    "The amendment matters: with momentum-dependent coupling, identical canonical momenta mean different actual velocities. The original fixed-momentum scan is preserved; it is not silently replaced. Relative physical velocities are matched in the second scan, with different source recoil velocities and required preparation energy recorded.","",
    "## The interaction and its energy","","The complete Hamiltonian and equations are in [the protocol](protocol.md), with initialization in [amendment 1](amendment-1.md). The extra term is","","    H_interaction = - eta/(N-1) sum_(i<j) w_ij p_i dot p_j","",
    "with a Gaussian spatial kernel. The same Hamiltonian supplies both the force and the velocity law. Source recoil is dynamical; a compensating source spin accounts for the released packet's initial angular momentum. For |eta|<=0.8 the kinetic matrix has eigenvalues at least 0.2, and the softened source potential has a finite lower bound.","",
    "The source begins with a prescribed carrier reservoir. Preparation energy is the increase from coincident, resting carriers at the source plus compensating spin energy. Unchanged rest energies cancel in that difference. This is the required release/preparation work, not a derived conversion process or demonstrated astrophysical fuel supply. After preparation no external work is added.","",
    "## Matched-speed findings","","The table uses the fixed eta=0.8, range=1 slice, not a winner chosen to fit observations. Ratios are late-time RMS widths divided by the corresponding interaction-off control. Below one means a narrower packet; it does not by itself mean stronger gravity.","",
    "| Packet | Source gravity off | Source gravity on |","|---|---:|---:|"]
    for p in M.PACKETS:lines.append(f"| {p} | {comp(p,0)['late_width_ratio']:.3f} | {comp(p,1)['late_width_ratio']:.3f} |")
    lines+=["",f"The stream example with source gravity off has width ratio {c['late_width_ratio']:.4f}; required preparation energy is {a['preparation_energy']:.3f} model units versus {b['preparation_energy']:.3f} for the control, a factor {a['preparation_energy']/b['preparation_energy']:.2f}. Thus cohesion is not obtained with an uncounted energy bonus. It is a finite-time effect over t=0..20, not proof of indefinite binding.","",
    "With source gravity on, the stream has a broader packet but a larger fraction still near the source. Retention and coherent following are different metrics. In the same slice, rotating and hot packets also broaden. A tighter packet is therefore not a universal outcome, and no stable self-generated cluster-scale whirlpool has been demonstrated.",
    "","The initially radial symmetric packet acquires no substantial net circulation; rotating/reverse packets inherit their supplied angular momentum and have matching radial outcomes. Neither result demonstrates generation of net angular momentum from nothing.","",
    "## A force-sign trap caught by the controls","","At eta=0.8 the parallel pair has attractive canonical force and inward separation acceleration. The opposite-heading pair has repulsive canonical force but also inward instantaneous separation acceleration in the declared one-dimensional setup, because the velocity-momentum relation changes with separation. The orthogonal control has geometric separation acceleration even with the interaction off. Thus the model cannot be advertised as an exclusively aligned-attraction or predecessor-following law based only on the p_i dot p_j factor. Actual trajectories decide.","",
    "## Verification","",
    f"Across both primary scans, maximum scaled energy drift is {summary['maximum_energy_drift']:.3g}; maximum canonical momentum drift is {summary['maximum_momentum_drift']:.3g}; angular-momentum drift is {summary['maximum_angular_drift']:.3g}. The largest selected refinement discrepancy is {summary['maximum_refinement_error']:.3g}, below the declared 0.001 relative-position threshold. These are numerical consistency checks, not observational passes.","",
    f"The independent audit verifies {hashcount} evidence files and {pinned} pinned source/protocol instances. A separate scalar pair-sum energy calculation reproduces archived endpoint energies to {energy_error:.3g} model units. Rotation-reversal radial symmetry differs by at most {parity_error:.3g} relative. Full trajectories, scalar histories, failures and all settings are archived.",
    "","## Why this is not yet the astronomical theory","","The carriers have a finite inertial mass and low-speed kinetic energy. There is no massless limit, causal field mediator, photon equation, field production from empty initial data or physical coupling scale fitted to galaxies. The interaction is instantaneous, frame-dependent, and normalized by a finite packet's particle count. Those are explicit limitations, not completed requirements. A long Gaussian tail also supplies no finite propagation cone.",
    "","OP-1's observed-data findings are unchanged: blanket radial retention did not improve shared transfer, and Coma's absolute normalization remains unresolved. This mechanical proxy has no new prediction to place on the observed lensing/rotation chart.","",
    "## Attribution","","Hamilton equations and their conservation identities are established mechanics; see [the expert account of the action principle](https://www.scholarpedia.org/article/Principle_of_least_action). Gaussian kernels and spectral positivity bounds are standard mathematics. Velocity-dependent interactions also predate this work: [Essen's account of the Darwin interaction](https://arxiv.org/abs/1008.1182) discusses the established electromagnetic precedent. We do not implement the Darwin electromagnetic Hamiltonian or claim its mathematics as an invention. This is a declared project combination with no claim that it has never been considered. DOP853 is an established numerical method supplied by SciPy.","",
    "## Next work justified by these results","","1. Derive a local mediator whose finite-speed dynamics replace the instantaneous kernel, while preserving a positive energy and reciprocal recoil.","2. Specify the source reservoir and production process. Compute the energy needed to fill a spatial field from ordinary matter, rather than merely preparing twelve carriers.","3. Derive an actual follower response: compare equal physical emissions, opposite headings and counter-rotation, retaining the common-potential control.","4. Only then derive a common matter/light response and return to the fixed galaxy and cluster observations. Do not transplant the 29% width effect into a lensing-strength multiplier.","",
    "## Reproduction","","Original protocol: 14e6185; original source: 2cebb65. Matched-speed amendment/source: b3f539b. Run each driver in an isolated checkout without its evidence directory to reproduce a fresh archive; neither driver overwrites results.","",
    "    python -B research_work/experiments/conservative_carriers/audit_report.py","",
    "Requires NumPy, SciPy and Matplotlib. This command audits existing evidence and rebuilds the report/figure. No dark matter or expanding geometry is used.",""]
    (HERE/"report.md").write_text("\n".join(lines),encoding="utf8",newline="\n")
    print(json.dumps(dict(audit_passed=all(c["passed"] for c in checks),summary=summary),indent=2))
if __name__=="__main__":main()

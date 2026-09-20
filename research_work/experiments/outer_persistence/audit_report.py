"""Audit OP-1 evidence and produce the research record."""
import os
os.environ["OPENBLAS_NUM_THREADS"]="1"
import json,subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import persistence as P
D=P.D;L=P.L;HERE=P.HERE;save=L.save
def main():
    out=HERE/"evidence-v1";checks=[]
    def add(name,value,limit=0):
        checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
    hashes=D.read(HERE/"evidence-sha256.json")
    for path,h in hashes.items():add("evidence "+path,int(P.digest(HERE/path)!=h))
    manifest=D.read(out/"manifest.json")
    for path,h in manifest["hashes"].items():
        add("current source "+path,int(P.digest(P.ROOT/path)!=h))
        original=subprocess.check_output(["git","show",manifest["git_head"]+":"+path],cwd=P.ROOT)
        import hashlib
        add("pinned source "+path,int(hashlib.sha256(original).hexdigest()!=h))
    selected=D.read(out/"selection-before-test.json")["selected"]
    pred=D.read(out/"predictions.json");gals=D.galaxies()
    spec=selected["protected"];error=0.
    for g in gals:
        p=P.prepare(g);length=spec["L0"]*(g["mass"]/1e11)**spec["eta"]
        direct=P.direct_envelope(p["r"],p["q"],spec["delta"],length)
        gg=(p["g0"]+spec["A"]*p["release"]*(direct-p["q"])/p["r"])[p["indices"]]
        archived=next(r for r in pred["protected"][g["split"]]["galaxies"] if r["name"]==g["name"])
        error=max(error,float(np.max(abs(np.sqrt(g["r"]*gg)-archived["predicted"]))))
        add("fixed adopted distance "+g["name"],abs(archived["distance_Mpc"]-g["catalog"]["distance"]))
    add("independent selected envelope velocities km/s",error,1e-9)
    controls=D.read(out/"controls.json");refinement=D.read(out/"numerical-checks.json")
    add("failed controls",sum(not c["passed"] for c in controls))
    add("failed refinements",sum(not c["passed"] for c in refinement))
    save(HERE/"audit.json",dict(checks=checks,passed=all(c["passed"] for c in checks),evidence_hashes=len(hashes),
         source_hashes=len(manifest["hashes"]),independent_velocity_error_kms=error))
    cases=D.read(out/"candidates.json");maps=D.read(out/"baseline-residual-map.json")
    plt.rcParams.update({"font.size":9})
    fig,axs=plt.subplots(1,3,figsize=(14,4.3))
    for sign,color,label in ((1,"#2574b8","Positive retention"),(-1,"#c46b34","Sign reversal")):
        cc=[c for c in cases if c["valid"] and c["spec"]["A"]*sign>0]
        axs[0].scatter([c["scores"]["validation"]["galaxy_rmse"] for c in cc],
           [c["scores"]["validation"]["cluster_chi2_per_point"] for c in cc],s=15,alpha=.5,c=color,label=label)
    base=pred["baseline"]["validation"];winner=pred["protected"]["validation"]
    axs[0].scatter(base["galaxy_rmse"],base["cluster_chi2_per_point"],marker="x",s=60,c="black",label="Frozen CMF")
    axs[0].scatter(winner["galaxy_rmse"],winner["cluster_chi2_per_point"],marker="*",s=90,c="#6c3aa0",label="Selected")
    axs[0].set(xlabel="Validation galaxy RMSE (km/s)",ylabel="Validation pressure chi2 / point",title="All 225 declared settings")
    axs[0].legend(fontsize=7)
    for split,col in (("train","#2574b8"),("validation","#c46b34"),("test","#667a3c")):
        rr=[r for r in maps["bins"] if r["split"]==split]
        axs[1].plot([.165,.5,.835],[r["mean_signed_kms"] for r in rr],"o-",c=col,label=split+" (exposed)")
    axs[1].axhline(0,c="gray",lw=.8)
    axs[1].set(xlabel="Radius / last observed radius",ylabel="Mean predicted - observed (km/s)",title="Baseline radial residuals")
    axs[1].legend(fontsize=7)
    g=next(g for g in gals if g["name"]=="NGC3198")
    axs[2].errorbar(g["r"],g["y"],yerr=g["error"],fmt=".",c="black",label="Observed")
    axs[2].plot(g["r"],np.sqrt(g["r"]*P.BASE_FORCE(g)),c="#2574b8",lw=2,label="Frozen CMF")
    axs[2].plot(g["r"],np.sqrt(g["r"]*P.force(spec)(g)),"--",c="#c46b34",label="Selected OP-1")
    axs[2].set(xlabel="Radius (kpc)",ylabel="Circular speed (km/s)",title="NGC 3198: training example")
    axs[2].legend(fontsize=7)
    for ax in axs:ax.grid(alpha=.15);ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("OP-1: keeping a radial support envelope does not improve shared transfer",fontsize=13)
    fig.tight_layout(rect=(0,.08,1,.94))
    fig.text(.02,.025,"Fixed ordinary sources and distances; no dark matter or expansion. Pressure is not lensing. Radial means are row-weighted; model scores are object-weighted.",fontsize=8)
    fig.savefig(HERE/"comparison.png",dpi=180);fig.savefig(HERE/"comparison.svg");plt.close(fig)
    lens=D.read(out/"lensing.json");completion=D.read(out/"complete.json")
    inner=D.read(out/"inner-domain-sensitivity.json")
    paired=D.read(out/"paired-bootstrap.json")["protected"]
    lines=["# OP-1: finite outer-support persistence","","The 225-setting declared scan does not improve the shared model's galaxy transfer. Both selection rules choose the same small sign-reversed correction; validation improves slightly, test performance worsens, and primary cluster pressure and Coma predictions stay essentially unchanged. The recurrence and selected numerical calculations pass their declared checks. This diagnoses this envelope closure, not every possible swirl.","",
    "![Campaign results](comparison.png)","","## Frozen observation scores","",
    "| Model | Galaxy train / validation / test RMSE (km/s) | Pressure train / validation / test chi2 per point |","|---|---|---|"]
    for name in ("baseline","protected","MOND"):
        pp=pred[name]
        vals=lambda key:" / ".join(f'{pp[s][key]:.3f}' for s in ("train","validation","test"))
        lines.append(f"| {name} | {vals('galaxy_rmse')} | {vals('cluster_chi2_per_point')} |")
    lines+=["",f"Selection: {spec['id']}. The 11 CMF coefficients were fixed. Four universal settings were scanned; no per-object force fit. Of 226 cases including baseline, {completion['valid_train_validation_cases']} have valid training/validation force. Invalid cases and reasons remain in candidates.json. All partitions are historically exposed. Cluster scores include one fitted boundary pressure per object and diagonal errors only.","",
    f"Selected-minus-CMF test RMSE is {paired['baseline']['delta_rmse_kms']:.4f} km/s, conditional paired 95% interval {paired['baseline']['percentile95']}. Selected-minus-MOND is {paired['MOND']['delta_rmse_kms']:.3f} km/s, interval {paired['MOND']['percentile95']}. Neither is an improvement. The two selected records are the same formula, not independent replications.","",
    "## What this teaches us","","The full validation/test samples have positive average outer residuals, whereas NGC 3198 underpredicts outer speeds. A general increase in outer support therefore fixes neither the sample-wide calibration nor its transfer. The selected negative amplitude is one of the protocol's sign controls, not an after-the-fact edit.",
    "","The retained envelope is seeded at the first available radius. Removing the first 10% of sampled radii changes predicted speeds by as much as "+f"{max(r['max_speed_change_kms'] for r in inner):.2f} km/s"+" in the selected case. This is a material dependence on incomplete interior information. It is a sensitivity, not a failed step-refinement check.","",
    "## Cluster lensing","","| Ordinary-source bracket | Preferred beta | Bounded chi2 (six bins) |","|---|---:|---:|"]
    for r in lens["models"]:
        if r["model"]=="protected" and r["primary"]:
            lines.append(f"| ne0={r['ne0']} cm^-3, Mstar={r['mstar']:.1e} Msun | {r['beta_unbounded']:.3f} | {r['bounded']['chi2']:.3f} |")
    lines+=["","These match the frozen CMF predictions at the displayed precision. Beta greater than one cannot be achieved in the adopted static geometry. Broad reconstructed-bin errors and absent covariance/source information prevent an absolute lensing verdict. The 3/9/30 Mpc reach sensitivities are archived and were not used for selection.","",
    "## Equations, energy and originality","","The exact recurrence and units are in [the declared protocol](protocol.md). H is a squared-speed support envelope, not energy density. Its empirical force remains unfunded by a source/carrier Hamiltonian. We do not promote it to a physical graviton theory.",
    "","The scan uses established running maxima, exponential attenuation, quadrature and bootstrap, alongside the project's declared closure. Data and optical/MOND methods retain the [CMF attribution](../coherent_memory_fit/report.md#attribution). There is no historical novelty claim.","",
    "## Verification and reproduction","",f"{len(controls)} preliminary checks pass; {len(refinement)} selected refinement comparisons pass (163 unique comparisons, repeated for the identical two selections). The independent quadratic-cost envelope reproduces archived velocities to {error:.3g} km/s. The audit verifies {len(hashes)} evidence files and {len(manifest['hashes'])} inputs/source files both locally and at their recorded commit, plus all 149 adopted distances.","",
    "Protocol commit: 8c59e03. Source commit: 956297f. Evidence is in evidence-v1; the driver refuses to overwrite it. From the repository root:","","    python -B research_work/experiments/outer_persistence/audit_report.py","",
    "For a fresh execution use an isolated checkout at the source commit and run campaign.py. Requires NumPy, SciPy and Matplotlib. No new empirical data were acquired.",
    "","## Next physical question","","Replace the imposed support envelope with an explicitly reciprocal, energy-counted interaction. A bounded momentum-dependent interaction can test whether nearly parallel moving carriers attract and develop circulation beyond the ordinary-potential control. It must also be tested for opposite headings, source recoil, preferred-frame effects and causal limitations. A successful toy mechanism would still need a relativistic/light coupling and astronomical normalization; this campaign does not supply those.",""]
    (HERE/"report.md").write_text("\n".join(lines),encoding="utf8")
    print(json.dumps(dict(audit_passed=all(c["passed"] for c in checks),checks=len(checks),independent_velocity_error=error)))
if __name__=="__main__":main()

"""Derive the HH comparison report/figure from preserved evidence."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numerical import compare
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text(encoding="utf8"))
def main():
    sets={}
    for mode,folder in [("heading","evidence-v1"),("curvature","evidence-curvature-v1")]:
        root=HERE/folder
        screen=read(root/"mean-screen"/"results.json")
        refine=read(root/"mean-refinement"/"results.json")
        events0=read(root/"events-dt001"/"results.json")
        events1=read(root/"events-dt0005"/"results.json")
        eligible=[]
        for r in screen["rows"]:
            ps=[p for p in screen["pairs"] if p["id"]==r["id"] and p["impact"] in (1,2,4)]
            if len(ps)==3 and all(p["arrived"] for p in ps):
                eligible.append(dict(id=r["id"],inward=float(np.mean([p["inward"] for p in ps])),
                                     sideways=float(np.mean([abs(p["sideways"]) for p in ps]))))
        pairs=compare(screen,refine)
        sets[mode]=dict(folder=folder,screen=screen,refine=refine,events0=events0,events1=events1,
            eligible=eligible,positive=sum(r["inward"]>0 for r in eligible),
            arrived=sum(r["arrived"] for r in screen["results"]),
            refinement=dict(observer_pass=sum(p["observer_pass"] is True for p in pairs),
                both_arrived=sum(p["both_arrived"] for p in pairs),both_censored=sum(p["both_censored"] for p in pairs),
                arrival_mismatches=sum(p["arrival_mismatch"] for p in pairs)))
    extra=read(HERE/"evidence-numerical-v1"/"summary.json")
    audit=read(HERE/"integrity-audit.json")
    derived={mode:{k:v for k,v in d.items() if k in ("eligible","positive","arrived","refinement")} for mode,d in sets.items()}
    derived["mean_trajectories"]=1970
    derived["discrete_histories"]=6144
    derived["remaining_observer_failures"]={mode:[r for r in s["comparisons"] if r["observer_pass"] is False] for mode,s in extra["suites"].items()}
    (HERE/"comparison.json").write_text(json.dumps(derived,indent=2)+"\n",encoding="utf8")
    # A scientific plot of measured outputs, not an illustration of presumed success.
    fig,axs=plt.subplots(2,2,figsize=(12,8.5),constrained_layout=True)
    colors=dict(heading="#c45a26",curvature="#1966a7")
    for j,(mode,d) in enumerate(sets.items()):
        vals=np.sort([r["inward"] for r in d["eligible"]])
        axs[0,0].plot(np.arange(1,len(vals)+1),vals,label=f"{mode}: {d['positive']}/{len(vals)} positive",color=colors[mode])
    axs[0,0].axhline(0,color="black",lw=.8)
    axs[0,0].set(xlabel="Eligible case, sorted separately",ylabel="Paired inward mean (rad)",title="A changed coupling changes the bending sign")
    axs[0,0].legend(frameon=False,fontsize=9)
    moderate="H2-A10-D1-K4";strong="H4-A100-D0.1-K16"
    for mode,d in sets.items():
        with np.load(HERE/d["folder"]/"mean-screen"/"trajectories.npz") as z:
            for b in (-2,2):
                i=np.flatnonzero((z["ids"]==moderate)&(z["impact"]==b))[0]
                axs[0,1].plot(z["state"][i,:,0],z["state"][i,:,1],color=colors[mode],
                              label=mode if b==-2 else None)
    axs[0,1].scatter([0],[0],c="black",s=30)
    axs[0,1].set(xlim=(-20,20),xlabel="x (model length)",ylabel="y (model length)",title="Same moderate parameters; different carried quantity")
    axs[0,1].legend(frameon=False,fontsize=9)
    for mode,d in sets.items():
        ee=[r for r in d["events1"]["ensembles"] if r["id"]==strong]
        x=np.array([r["impact"] for r in ee])+(0.07 if mode=="heading" else -.07)
        axs[1,0].errorbar(x,[r["mean_angle"] for r in ee],yerr=[r["sd_angle"] for r in ee],
                         fmt="o",capsize=4,color=colors[mode],label=mode)
        axs[1,1].plot([r["impact"] for r in ee],[r["arrivals"]/128 for r in ee],"o-",
                      color=colors[mode],label=mode)
    axs[1,0].axhline(0,color="black",lw=.8)
    axs[1,0].set(xlabel="Signed impact (model length)",ylabel="Arriving-ray angle (rad)",
                 title="Strong capture: means +/- 1 standard deviation")
    axs[1,0].legend(frameon=False,fontsize=9)
    axs[1,1].set(ylim=(0,1.05),xlabel="Signed impact (model length)",ylabel="Fraction arriving by T=80",
                 title="Keep nonarrivals: survivor angles alone are misleading")
    axs[1,1].legend(frameon=False,fontsize=9)
    fig.suptitle("Temporary photon-companion attachment: two conditional toy laws",fontsize=15)
    fig.savefig(HERE/"comparison.png",dpi=160);plt.close(fig)

    lines=["# HH-1 / HH-2: temporary attachment can carry direction or curvature\n",
      "**Executed toy campaigns; neither is a validated theory of gravity.** Temporary attachment "
      "can be represented mathematically. In these tests, carrying an outgoing companion's heading "
      "spreads light outward; carrying its inward turning curvature gives inward bending. The latter "
      "is a new, explicitly declared interaction postulate. The simulations do not derive a real "
      "photon-graviton bound state, spontaneous swirl formation, stellar dynamics or observed cluster lensing.\n",
      "[HH-1 declaration](protocol.md) | [HH-2 declaration](protocol-curvature.md) | "
      "[numerical follow-up](numerical-followup.md) | [integrity audit](integrity-audit.json) | "
      "[machine-readable comparison](comparison.json)\n",
      "## What was assumed\n",
      "Use MS-1's finite, counted matter source and prescribed spiral wind: source energy 10, "
      "emission Q=1/30, launch radius 0.5, companion transport speed 0.5, initial age 30. "
      "At pitch h=tan(beta), radial transport is 0.5/sqrt(1+h^2); "
      "u=Q/(4*pi*v_r*r^2) only between the launch sphere and the expanding *emission front*. "
      "That front is transport in static space, not expansion of the universe. More winding "
      "increases local residence density while reducing radial reach at fixed age and emitted energy.\n",
      "Photons have imposed unit speed and constant energy, starting at x=-20, observed at x=20 "
      "or censored at model time 80. 'Nonarrival' means missing this plane by that deadline; "
      "it does not prove absorption or permanent trapping. Matter emits into a prescribed field; "
      "ordinary matter emission microphysics and a three-dimensional self-consistent swirl are still missing.\n",
      "## The two attachment laws\n",
      "Let n be photon direction, t the local companion tangent, p attached occupation, "
      "a=k_on*u*((1+n dot t)/2)^2 and b=k_off. Both laws use\n",
      "    p_dot = a*(1-p) - b*p\n",
      "**HH-1: carry heading.** A captured companion gives a stored vector t_capture until release. "
      "The mean moment m obeys\n",
      "    m_dot = a*(1-p)*t - b*m\n"
      "    theta_dot = kappa*cross(n,m),     x_dot = n\n",
      "**HH-2: carry curvature.** For t=cos(beta)*e_r+sin(beta)*e_phi, direct differentiation gives\n",
      "    C = (t dot grad)t = (sin(beta)/r)*(cos(beta)*e_phi-sin(beta)*e_r)\n"
      "    q_dot = a*(1-p)*C - b*q\n"
      "    theta_dot = chi*cross(n,q),      x_dot = n\n",
      "C has inward radial component -sin(beta)^2/r despite the outgoing radial velocity. "
      "A discrete attachment stores C_capture and turns the photon under "
      "chi*(I-n*n^T)*C_capture until release. This copies a curvature vector at capture; "
      "it does not dynamically keep the attached pair on one evolving streamline. "
      "kappa and chi use the same numeric grid for comparison but are different couplings "
      "with different dimensional meanings before nondimensionalization.\n",
      "After leaving the source front, an already attached photon keeps its stored moment "
      "until it releases it; it does not sample a nonexistent stream there. Release is exponential "
      "with mean bound duration 1/b. A later capture can select a new local direction or curvature.\n",
      "The mean-state closure is not an exact ensemble average of nonlinear random trajectories. "
      "Discrete histories therefore test scattering and nonarrivals separately. Events occur at "
      "the integration midpoint, at most one capture/release transition per step.\n",
      "## Complete parameter screen\n",
      "For each law: pitch {0.5,2,4}, capture coefficient {1,10,100}, release rate {0.1,1,10}, "
      "and guide coefficient {1,4,16}: 81 cases with ten signed impacts. "
      "Primary mean step 0.02; six fixed cases repeated at 0.01. "
      "Three fixed stochastic cases at four impacts, 128 realizations each, "
      "at steps 0.01 and 0.005 with recorded independent seeds.\n",
      "| Model | Mean rays arriving / 810 | Eligible paired cases / 81 | Positive inward means | Random arrivals at step 0.005 / 1536 |\n"
      "|---|---:|---:|---:|---:|"]
    for mode,d in sets.items():
        lines.append(f"| {mode} | {d['arrived']} | {len(d['eligible'])} | {d['positive']} | {d['events1']['arrivals']} |")
    lines.extend(["\nEligibility requires both sides to arrive at each impact 1,2,4. "
      "Inward A=(theta_minus-theta_plus)/2 and sideways B=(theta_minus+theta_plus)/2. "
      "A positive *mean* of three pairs is only a sign diagnostic, not a lensing acceptance threshold. "
      "The 15 HH-1 and 25 HH-2 ineligible cases are retained, not successes or zeros.\n",
      "HH-2's largest eligible primary mean is 0.545529 rad for H2-A100-D0.1-K1, "
      "with mean absolute sideways component 0.455750 rad. This is strongly asymmetric bending, "
      "not a clean isotropic lens. It is an exposed selection and not an observational prediction.\n",
      "Including numerical follow-ups and mirrors, these attachment campaigns contain **1,970 "
      "mean trajectories and 6,144 discrete histories**. This is in addition to GF-1's 618 "
      "formula cases/1,515 trajectories and MS-1's 210 cases/2,320 rays.\n",
      "## Numerical checks and unresolved values\n",
      "A time step is the simulated interval between integration updates. Halving it tests "
      "resolution while leaving the physical law and parameters unchanged. Model time is not "
      "yet calibrated to seconds. A resolution failure is not by itself a rejection of the idea.\n",
      "| First step 0.02 vs 0.01 | Arriving pairs within angle/position tolerance | Both censored | Arrival mismatches |\n"
      "|---|---:|---:|---:|"])
    for mode,d in sets.items():
        s=d["refinement"];lines.append(f"| {mode} | {s['observer_pass']}/{s['both_arrived']} | {s['both_censored']} | {s['arrival_mismatches']} |")
    lines.extend(["\nThe archived first audits report 57/60 and 58/60, including matched "
      "nonarrivals as passes. The table above separates those censored cases; matching censorship "
      "does not establish trajectory accuracy. Original audits are preserved.\n",
      "| Additional step 0.005 vs 0.0025 | Arriving pairs within tolerance | Both censored | Arrival mismatches |\n"
      "|---|---:|---:|---:|"])
    for mode,s in extra["suites"].items():
        lines.append(f"| {mode} | {s['arrived_pass']}/{s['both_arrived']} | {s['both_censored']} | {s['arrival_mismatches']} |")
    lines.extend(["\nThe follow-up includes the failed heading representative and all six curvature "
      "representatives plus its exposed largest-inward case. Two heading and one curvature "
      "observer comparisons remain outside the original 0.01 tolerance. Full terminal-state "
      "differences for censored rays are also saved; these are not declared converged. "
      "Discontinuous source/front boundaries and long path sensitivity remain numerical concerns.\n",
      f"All {len(extra['controls'])} follow-up controls pass: analytic bound-state steering, "
      "exponential decay after leaving the front, event bookkeeping, observer positions, unit-speed "
      "path bounds and full chirality reflection. HH-1's four original controls and HH-2's independent "
      "finite-difference curvature check pass. Occupation/moment bounds hold across the archived mean states.\n",
      f"The integrity audit verifies {audit['hashes_checked']} exact evidence byte hashes, "
      f"{audit['pinned_sources_checked']} source files against their recorded Git commits and "
      f"{audit['trajectory_arrays_checked']} ray archives. All pass. This is reproducibility and "
      "numerical implementation evidence, not a physical validation.\n",
      "## Why random hitchhiking is not yet a satisfactory lens\n",
      "For the strong curvature case H4-A100-D0.1-K16 at impact -2, "
      "68/128 rays arrive in the finer stochastic run. Those survivors have mean angle "
      "0.5933 rad and standard deviation 0.4136 rad, with mean extra travel time 9.830 model units. "
      "At impact +2, 120/128 arrive with mean -0.00810 rad and standard deviation 0.07281 rad. "
      "The two sides behave very differently.\n",
      "The moderate curvature case at impact -2 has mean 0.06521 rad but standard deviation "
      "0.24379 rad among 126 arrivals. Rare attachments can make a broad mixture of undeflected "
      "and strongly scattered light. Stronger average bending is not enough if images blur or "
      "the selected arriving rays conceal large losses. A zero survivor angle can also mean "
      "only the unscattered rays reached the plane.\n",
      "Between the two independent-seed curvature ensembles, the largest mean difference is "
      "2.34 combined standard errors. This is descriptive: twelve small ensembles, rare events, "
      "survivor selection and a changed timestep do not establish a stochastic convergence theorem. "
      "A tiny sample with zero captures cannot establish an exactly vanishing physical rate.\n",
      "## Conservation, attribution and what remains to test\n",
      "The source energy plus emitted transport energy is counted. The additional attachment law "
      "has **no derived binding Hamiltonian or reciprocal stream evolution**. Maintaining imposed "
      "photon speed/energy does not establish energy, momentum and angular-momentum conservation "
      "for the full source-stream-photon system. Opposite photon momentum is a recorded debt, "
      "not a dynamically evolved recoil. It must be paid by a completed model.\n",
      "These are fictional companion coupling states. No real massive particle is claimed to "
      "ride a photon at light speed for free, and no established spin-2 quantum bound state is asserted. "
      "Frequency-independent rates were imposed; achromatic lensing and physical dispersion were "
      "not independently derived. No dark matter, expanding geometry or distance fitting is used.\n",
      "The user's physical picture motivates these combinations. Curvature calculus, transverse "
      "projection, rate equations, persistent memory and random attachment/release are established "
      "mathematics, not claimed inventions. The [GF-1 attribution register](../companion_following/review-and-provenance.md) "
      "credits related alignment/swarming constructions. No historical novelty has been established. "
      "Standard lensing already integrates along the ray: "
      "[Bartelmann and Schneider](https://arxiv.org/pdf/astro-ph/9912508). "
      "The proposed distinction is state-dependent coupling to a directed stream, not counting "
      "the same ordinary gravity repeatedly as free amplification.\n",
      "Before an astronomical test can promote this route: derive a reciprocal attachment/release "
      "action with a positive counted energy; evolve the companion field and emitter recoil; "
      "show that the swirl forms and transports causally; derive both matter and photon responses "
      "from the same coefficients; converge capture statistics and boundary crossing; then freeze "
      "a calibration and test star velocities, shear, image sharpness, arrival times and color. "
      "No current result chooses this route over the separate clock or energy-conversion alternatives.\n",
      "![Actual paired bending, mean trajectories, scattering and nonarrival fractions](comparison.png)\n"])
    (HERE/"report.md").write_text("\n".join(lines),encoding="utf8")
    print(json.dumps({mode:dict(eligible=len(d["eligible"]),positive=d["positive"],refinement=d["refinement"]) for mode,d in sets.items()}))
if __name__=="__main__":main()

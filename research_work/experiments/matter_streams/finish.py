"""Preserve MS-1 first execution and refine its explicitly unresolved rays."""
import hashlib,json
from pathlib import Path
import numpy as np
from campaign import HERE,dump,trace

def main():
    root=HERE/"evidence-refinement";root.mkdir(exist_ok=False)
    import subprocess
    dump(root/"manifest.json",dict(git_head=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir() if p.suffix in (".py",".md")}))
    first=json.loads((HERE/"evidence-v1"/"screen"/"results.json").read_text())
    fine=json.loads((HERE/"evidence-v1"/"refinement"/"results.json").read_text())
    ids=[r["id"] for r in first["ranking"][:6]]
    positive=[r for r in first["ranking"] if next(c for c in first["rows"] if c["id"]==r["id"])["kappa"]>0]
    ids += [positive[0]["id"],positive[-1]["id"]]
    rows=[next(c for c in first["rows"] if c["id"]==i) for i in ids]
    finer,_=trace(rows,root/"dt0005",.005)
    finest,_=trace(rows,root/"dt00025",.0025)
    old={(r["id"],r["impact"]):r for r in finer["results"]}
    checks=[]
    for r in finest["results"]:
        q=old[(r["id"],r["impact"])]
        same=q["arrived"]==r["arrived"]
        da=abs(q["exit_angle"]-r["exit_angle"]) if same and r["arrived"] else None
        dy=abs(q["observer_y"]-r["observer_y"]) if same and r["arrived"] else None
        de=abs(q["exposure"]-r["exposure"])/max(1,abs(r["exposure"]))
        checks.append(dict(id=r["id"],impact=r["impact"],angle_difference=da,position_difference=dy,
            exposure_difference=de,passed=same and (da is None or da<.01) and (dy is None or dy<.01) and de<.01))
    summary=dict(comparisons=checks,passing=sum(r["passed"] for r in checks),count=len(checks),
        positive_coupling_cases=len(positive),positive_coupling_inward=sum(r["mean_inward"]>0 for r in positive),
        best_positive=positive[0],strongest_negative=first["ranking"][0],
        initial_unconverged=8,first_command_exit=1,
        first_command_error="Console JSON serialization after every evidence artifact and complete marker had been saved.")
    dump(root/"summary.json",summary);dump(root/"complete.json",dict(complete=True))
    lines=["# MS-1: cumulative steering in matter-fed streams\n",
      "The owner's along-the-river hypothesis was made explicit as directional turning accumulated along a ray. "
      "A counted source emits companions into a prescribed spiral wind with an evolving outer front. "
      "The stream geometry is assumed, not self-organized from the GF-1 packets.\n",
      "Ran 210 parameter cases and 2,100 primary rays, 60 first refinements, and 160 finer refinements. "
      "No dark matter, expansion, distance fit or observational lensing fit. All primary rays reached the observer plane.\n",
      f"**Following an outward stream did not give the required inward lensing in this law.** "
      f"All {len(positive)} positive-coupling cases have nonpositive mean paired inward deflection at impacts 1,2,4. "
      "The rays can turn and sweep sideways, but outward following is not automatically attraction toward the source.\n",
      f"The largest primary inward mean, {first['ranking'][0]['mean_inward']:.6g} rad, occurs for "
      f"{first['ranking'][0]['id']}: **negative coupling**, which turns light away from the outgoing heading. "
      "That is an exploratory counter-following law, not evidence that positive stream following solved lensing. "
      "These dimensionless angles are toy outputs, not predictions in arcseconds for real clusters.\n",
      "## What was varied\n",
      "Five spiral pitches, both handednesses, three alignment gates, and seven signed coupling values "
      "(-16,-4,-1,0,1,4,16). Each has ten opposite-side impact parameters. "
      "The direction law is theta_dot=kappa*u_g*F(cos(psi-theta))*sin(psi-theta). "
      "Straight aligned flow produces no bend; sustained deflection needs changing stream direction or heading mismatch. "
      "Inward paired deflection and common sideways motion are scored separately.\n",
      "At the same age every pitch contains exactly the same emitted energy Q*(age0+t). "
      "Tighter winding raises residence density but shortens radial reach. Source plus emitted energy remains 10. "
      "This is an energy transport accounting control, not a derived microscopic emission or torque mechanism.\n",
      "## Verification and preserved failures\n",
      "Seven preliminary controls pass; full-grid chirality reflection is exact, and zero coupling gives a straight "
      "ray to about 7e-13 numerical error. Only 52/60 initial refinement comparisons pass: the sharp source/front "
      "boundaries leave eight unresolved at dt0.01. Those failed scores remain unchanged.\n",
      f"Additional dt0.005 versus dt0.0025 comparisons pass {summary['passing']}/{summary['count']} "
      "under the original angle/position/exposure limits. The set includes all six selected counter-following laws "
      "and two positive-coupling representatives. Unconverged values must not be promoted.\n",
      "The first command returned exit1 solely while printing a NumPy Boolean as JSON, after all results, the plot "
      "and completion marker had been saved. The console serialization was corrected; the first simulation archive "
      "was not overwritten or rerun to hide it.\n",
      "## Limits and next extension\n",
      "Photon speed and energy are fixed by the steering ansatz. The opposite momentum impulse is recorded as a debt "
      "to the prescribed stream, not integrated as a reciprocal field response. Emission microphysics, neighbor-led "
      "wind formation, emitter torque, a common stellar-force law and actual cluster shear remain unresolved. "
      "This result limits this particular outgoing alignment law, not all coherent stream interactions.\n",
      "The requested attach-detach idea is a separate [HH-1 experiment](../photon_hitchhiking/protocol.md). "
      "It tests whether a carried direction memory differs from continually sampling the current local stream.\n",
      "Prior work is credited in [the protocol](protocol.md). Standard gravitational lensing already accumulates "
      "along a path: [Bartelmann and Schneider](https://arxiv.org/pdf/astro-ph/9912508). "
      "The proposed distinction here is an added direction-dependent coupling, not discovery of path integration.\n",
      "![Actual MS-1 rays](results.png)\n"]
    (HERE/"report.md").write_text("\n".join(lines),encoding="utf8")
    dump(HERE/"evidence-sha256.json",{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
        for folder in (HERE/"evidence-v1",root) for p in folder.rglob("*") if p.is_file()})
    print((root/"summary.json").read_text(encoding="utf8"))
if __name__=="__main__":main()

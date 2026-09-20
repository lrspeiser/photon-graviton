"""Render SR-1 evidence, explicitly separating qualification and physics."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent


def main():
    original=json.loads((ROOT/'spatial-v1/summary.json').read_text())
    current=json.loads((ROOT/'spatial-v2/summary.json').read_text())
    audit=json.loads((ROOT/'spatial-audit.json').read_text())
    controls=json.loads((ROOT/'spatial-v1/controls.json').read_text())
    if not audit['passed']:raise RuntimeError('Independent audit must pass before rendering')
    rows=current['runs'];byname={r['config']['name']:r for r in rows}
    lines=['# SR-1: reciprocal spatial response in coupled 3D evolution','',
           '20 September 2026. Synthetic dimensionless source experiment; no observational data or fitted cluster/galaxy parameters. No dark matter or cosmic expansion. The full twelve-item solution remains incomplete.','',
           '## What changed','',
           'The same constitutive spatial coefficient now enters particle energy, light speed, field kinetic/gradient energy, and all reciprocal field/source derivatives. s=0 recovers CC-2. s=1 adds equal weak temporal/spatial response; it is an established metric-optics choice, not an independently fitted lens multiplier or new historical invention. Lambda=0 in all evolutions; lambda=.5 is exercised only in derivative controls.','',
           'For stationary beta=0 and weak U=-C0/r, slow-particle acceleration is -grad U to leading order. The photon Hamiltonian is exp((1+s)U)|p|. Along a straight zeroth-order ray at impact b, integrating (1+s) C0 b/(b^2+z^2)^(3/2) over z gives deflection 2(1+s)C0/b. Therefore s=1 restores the usual factor 4C0/b at fixed orbital normalization. This analytic weak static limit does not imply that the dynamic source below yields an exact factor of two.','',
           'The frozen field principal block has a positive symmetrizer and characteristic speeds beta.n +/- C; photons have that same local cone. This is evidence of frozen principal hyperbolicity, not a proof of full nonlinear stability or stable galactic orbits. Finite source averaging gives an averaged cone and remains a nonlocal regulator.','',
           '## Numerical evidence','',
           f"- Local controls: {sum(r['passed'] for r in controls['local'])}/240; coupled-grid directional derivatives: {sum(r['passed'] for r in controls['grid'])}/24; four omitted-self-source negative controls distinguish the incorrect equations.",
           f"- CC-2 equivalence: energy error {controls['cc2_equivalence']['energy']:.6g}, maximum RHS error {controls['cc2_equivalence']['rhs']:.6g}.",
           f"- Original box: {sum(r['passed'] for r in original['runs'])}/8 trajectory qualifications pass. The failed conservative boundary clearance is preserved, not waived because energy drift is small.",
           f"- Enlarged boxes: {sum(r['passed'] for r in rows)}/8 trajectory qualifications pass; {sum(r['passed'] for r in current['comparisons'])}/10 time, space and domain comparisons pass.",
           f"- Independent archive audit: {audit['count']} checks pass; maximum raw-state energy reconstruction error {audit['max_independent_energy_error']:.6g}; maximum final Hamiltonian directional-derivative error {audit['max_independent_gradient_error']:.6g}.",
           f"- Maximum scaled H+Q drift {max(r['ledger_drift'] for r in rows):.6g}; averaged-cone error {max(r['cone_error'] for r in rows):.6g}; minimum geometric clearance {min(r['clearance'] for r in rows):.6g}.",
           f"- Maximum absolute total momentum drift {max(r['momentum_drift'] for r in rows):.6g}; angular momentum drift {max(r['angular_drift'] for r in rows):.6g}. Grid angular momentum is not an exact invariant.",'',
           'The clearance uses the largest sampled characteristic speed and source extent, and distance to the damping-layer entrance. It is a conservative continuum travel estimate plus a numerical edge-amplitude check, not a proof of a perfectly absorbing boundary. No long-time boundary claim is made.','',
           '| Source/control | s | Bend (mrad) | Field energy | Matter energy lost | Absorbed energy |','|---|---:|---:|---:|---:|---:|']
    for row in rows[:6]:
        lines.append(f"| {row['config']['name']} | {row['config']['s']} | {1000*row['bend']:.9g} | {row['final']['field']:.9g} | {row['initial']['matter']-row['final']['matter']:.9g} | {row['final']['absorbed']:.4g} |")
    lines+=['','Negative signed bends point toward the source for this geometry. Values are relative to the initial ray heading. The fast source has momentum/mass=.2 (speed about .196c); the slow source is set to 200 km/s. The photon retains its test energy 1e-4, so these slow-source runs are not yet a negligible-photon-backreaction astrophysical limit.','',
            '| Refinement/comparison | Probe position difference | Relative field-energy difference | Passed |','|---|---:|---:|---|']
    for row in current['comparisons']:
        lines.append(f"| {row['name']} | {row['probe_position_difference']:.6g} | {row['field_energy_relative_difference']:.6g} | {row['passed']} |")
    bend_resolution=abs(byname['s1-fast']['bend']/byname['s1-fast-space']['bend']-1)
    lines+=['',f'The photon bend itself differs by {100*bend_resolution:.6g}% relative to the finer-grid result. The declared position/field-energy gates passing does not establish percent-level lensing accuracy. More spatial refinements and source-radius tests are required before interpreting a precise bend prediction.']
    lines+=['','## Interpretation and remaining rejection tests','']
    for tag in ('scalar','fast','slow'):
        zero=byname['s0-'+tag]['bend'];one=byname['s1-'+tag]['bend']
        lines.append(f"- {tag}: changing s=0 to s=1 changes the bend magnitude by a factor {abs(one/zero):.6g} in this particular dynamical experiment.")
    for s in (0,1):
        scalar=byname[f's{s}-scalar']['bend'];vector=byname[f's{s}-fast']['bend']
        lines.append(f"- At s={s}, switching on the vector coupling in the fast-source pair changes bend magnitude by {100*(abs(vector/scalar)-1):.6g}%. This includes reciprocal source/field changes, not just a ray deflection in a frozen field.")
    lines+=['','These are controlled mechanism comparisons, not agreement with cluster observations. No eta=0 slow-source run was declared, so the slow-source values do not isolate the vector contribution. No source-sense reversal or rotated spatial-response run was included; CC-2 rotation results cannot be transferred automatically. The coarse/fine spatial comparison also changes domain length; the eight matched-resolution domain comparisons quantify that dependency.','',
            'The main unresolved structural problem survives: an isolated compact source in the linear exterior still gives Kepler/Yukawa declining rotation. Increasing the spatial optical response does not generate an extended, source-funded force. The direct-current suppression, amplitude-threshold sensitivity, physical source budget and long-lived swirl all remain open. The next meaningful development is to derive an extended-field support mechanism with a single universal law and test formation, stability, outer force and local limits; long runs still require a qualified massive/nonlinear 3D outgoing boundary. Do not proceed to claimed observational success by fitting an independent lensing multiplier.','',
            '## Reproduction and attribution','',
            'The first archives are immutable: run_spatial.py creates spatial-v1 and run_spatial_enlarged.py creates spatial-v2 only if absent. To reproduce without removing archived evidence, create a separate checkout at implementation commit 0ba3c92 (which predates both archives), then run Python with -B on run_spatial.py and run_spatial_enlarged.py in that order. The numerical source hashes are unchanged between that commit and the original runs. Use audit_spatial.py from the completed campaign to inspect the regenerated evidence. In the current checkout, run audit_spatial.py against checked-in evidence, then report_spatial.py for derived outputs. Manifests pin protocol and numerical source hashes and Git commits. The auditor reconstructs energy directly from raw arrays without calling solver energy/ingredients, and compares final-state Hamiltonian gradients against the solver RHS.','',
            '| Component | Attribution/status |','|---|---|',
            '| Hamiltonian matter/light evolution and optical metric | Established mathematics; [Gibbons et al.](https://arxiv.org/abs/0811.2877). |',
            '| Equal weak temporal/spatial response s=1 | Established metric response, explicitly adopted as a constitutive choice. |',
            '| Optional amplitude coupling | Related to [scalarization](https://arxiv.org/abs/gr-qc/9602056) and [vectorization](https://arxiv.org/abs/1706.01056); derivative-tested only here. |',
            '| This effective scalar/vector Hamiltonian and reciprocal discrete campaign | Specific candidate construction and measured numerical results; historical uniqueness is not claimed. |','',
            '![Synthetic spatial-response comparisons](spatial.png)','']
    (ROOT/'spatial-report.md').write_text('\n'.join(lines),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(11,5))
    x=np.arange(3)
    for s,color in ((0,'#617d98'),(1,'#007f73')):
        values=[1000*abs(byname[f's{s}-{tag}']['bend']) for tag in ('scalar','fast','slow')]
        axes[0].bar(x+(s-.5)*.34,values,width=.34,label=f's = {s}',color=color)
    axes[0].set_xticks(x,['Fast, scalar only','Fast, scalar + vector','Slow, scalar + vector']);axes[0].tick_params(axis='x',labelsize=8)
    axes[0].set_ylabel('Deflection magnitude (mrad)');axes[0].legend(ncol=2);axes[0].set_ylim(0,8.3);axes[0].set_title('One common response changes the bend')
    for s,color in ((0,'#617d98'),(1,'#007f73')):
        trace=json.loads((ROOT/f'spatial-v2/s{s}-fast.json').read_text())['trace']
        axes[1].plot([r['time'] for r in trace],[r['field'] for r in trace],label=f'Field energy, s={s}',color=color)
        axes[1].plot([r['time'] for r in trace],[trace[0]['matter']-r['matter'] for r in trace],ls='--',color=color,label=f'Matter loss, s={s}')
    axes[1].set_xlabel('Model time');axes[1].set_ylabel('Dimensionless energy');axes[1].legend(fontsize=8);axes[1].set_title('Field energy is supplied by matter')
    fig.suptitle('SR-1: synthetic 3D experiment — no observational comparison',fontsize=12)
    fig.text(.5,.015,f'Illustrative bends: the s=1 fast-source bend changes by {100*bend_resolution:.1f}% with spatial refinement.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.96));fig.savefig(ROOT/'spatial.png',dpi=170);plt.close(fig)
    print('Wrote spatial-report.md and spatial.png')


if __name__=='__main__':main()

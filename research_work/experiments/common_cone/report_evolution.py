"""Render the completed CC-2 archive without altering the scientific evidence."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;OUT=HERE/'evolution-v1'


def main():
    rows=json.loads((OUT/'runs.json').read_text());comparisons=json.loads((OUT/'comparisons.json').read_text())
    audit=json.loads((HERE/'evolution-audit.json').read_text());controls=json.loads((OUT/'controls.json').read_text())
    if len(rows)!=31 or not audit['audit_passed']:raise AssertionError('Require complete, independently audited campaign')
    byid={r['spec']['id']:r for r in rows}
    find=lambda source,eta,probe:byid[f'{source}-eta{eta:g}-{probe}']
    example=find('rotating',.08,'none')
    fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
    energy=np.array([[m[k] for k in ('field','matter','absorbed')] for m in example['metrics']])
    axes[0,0].plot(example['times'],energy[:,0]+energy[:,2],label='Field + accounted outlet',lw=3)
    axes[0,0].plot(example['times'],energy[0,1]-energy[:,1],'--',label='Matter energy lost',lw=2)
    axes[0,0].set(title='A. Matter funds the field',xlabel='Time',ylabel='Model energy');axes[0,0].legend(fontsize=8)
    y=np.load(OUT/(example['spec']['id']+'.npz'))['final_state'];n=32;dx=12/n
    f=y[:4*n**3].reshape((4,n,n,n));c=np.arange(n)*dx-6
    curl=(np.roll(f[2],-1,0)-np.roll(f[2],1,0)-np.roll(f[1],-1,1)+np.roll(f[1],1,1))/(2*dx)
    im=axes[0,1].pcolormesh(c,c,curl[:,:,n//2].T,cmap='RdBu_r',shading='auto')
    axes[0,1].set(title='B. curl A in the midplane at t=2.5',xlim=(-3,3),ylim=(-3,3),aspect='equal',xlabel='x',ylabel='y')
    fig.colorbar(im,ax=axes[0,1])
    colors={'rest':'#657687','rotating':'#128594','reverse':'#c16539'}
    for source in ('rest','rotating','reverse'):
        for eta in (0.,.08):
            r=find(source,eta,'photon')
            axes[0,2].plot(r['times'],[m['probe_angle']*1000 for m in r['metrics']],color=colors[source],ls='--' if eta==0 else '-',label=f'{source}, eta={eta:g}')
    axes[0,2].set(title='C. Derived photon direction',xlabel='Time',ylabel='Angle (mrad)');axes[0,2].legend(fontsize=7)
    for kind,color in (('time','#657687'),('space','#128594'),('domain','#c16539')):
        r=[x for x in comparisons if x['kind']==kind]
        axes[1,0].semilogy(range(3),[max(x['position'],1e-16) for x in r],'o-',label=kind,color=color)
    axes[1,0].set(title='D. Trajectory comparison',xticks=range(3),xticklabels=['rest','rotating','reverse'],ylabel='Maximum position difference');axes[1,0].legend(fontsize=8)
    for name,label in (('rotating-eta0.08-photon','32 cells'),('rotating-space','40 cells'),('rotated-n32','32 cells, rotated'),('rotated-n40','40 cells, rotated')):
        r=byid[name];l0=np.array(r['metrics'][0]['angular'])
        axes[1,1].plot(r['times'],[np.linalg.norm(np.array(m['angular'])-l0) for m in r['metrics']],label=label)
    axes[1,1].set(title='E. Continuum angular-momentum residual',xlabel='Time',ylabel='Absolute residual');axes[1,1].legend(fontsize=7)
    for name,label in (('radius-0.6','a=0.6'),('rotating-eta0.08-photon','a=0.9'),('radius-1.2','a=1.2')):
        r=byid[name];axes[1,2].plot(r['times'],[m['probe_angle']*1000 for m in r['metrics']],label=label)
    axes[1,2].set(title='F. Physical source-radius sensitivity',xlabel='Time',ylabel='Photon angle (mrad)');axes[1,2].legend(fontsize=8)
    fig.suptitle('CC-2: nonlinear 3D common-cone control campaign\nShort transient with finite source smoothing; no astronomical fit',fontsize=15)
    fig.savefig(HERE/'evolution.png',dpi=170);plt.close(fig)
    table=[]
    for source in ('rest','rotating','reverse'):
        scalar=find(source,0.,'photon')['metrics'][-1]['probe_angle']
        vector=find(source,.08,'photon')['metrics'][-1]['probe_angle']
        table.append(f'| {source} | {scalar*1000:.6f} | {vector*1000:.6f} | {(abs(vector)/abs(scalar)-1)*100:.3f}% |')
    comparison_table=[]
    for r in comparisons:
        label=(r.get('source','')+' '+r['kind']+(' N='+str(r['n']) if 'n' in r else '')).strip()
        energy=f"{r['field_energy_relative']*100:.6g}%" if 'field_energy_relative' in r else 'not a declared gate'
        comparison_table.append(f"| {label} | {r['position']:.6g} | {energy} | {'PASS' if r['passed'] else 'FAIL'} |")
    radius_table=[]
    for name in ('radius-0.6','rotating-eta0.08-photon','radius-1.2'):
        r=byid[name];last=r['metrics'][-1]
        radius_table.append(f"| {r['spec'].get('radius',.9)} | {last['field']:.8g} | {last['probe_angle']*1000:.6f} |")
    e=example['metrics'][-1];m0=example['metrics'][0]['matter']
    report=f'''# CC-2: nonlinear three-dimensional source and wave evolution

Completed 20 September 2026. [Protocol](evolution-protocol.md), [solver](evolution.py), [driver](run_evolution.py), [independent audit](evolution-audit.json), [chart](evolution.png).

All 31 declared evolutions completed. {audit['evolution_gate_passes']}/31 pass the energy-ledger and averaged-cone gates; {audit['comparison_passes']}/{audit['comparison_count']} refinement/domain/rotation comparisons pass. All {len(controls['checks'])} full discrete-Hamiltonian controls and {len(audit['checks'])} independent provenance/endpoint checks pass. This is an actual coupled 3D calculation, but a short dimensionless control experiment, not a persistent galactic swirl or a joint observational solution.

![CC-2 nonlinear 3D diagnostics](evolution.png)

## What the solver includes

Six moving ordinary-matter sources, a spherical finite-radius interpolation/deposition profile, four dynamic scalar/vector fields, optional massive or massless probe, all field self-interaction derivatives, source recoil, and a damping-energy ledger are evolved together. Fields initially contain zero energy. The Hamiltonian determines both the source term and the matter/light response. No independent optical multiplier, inserted halo, dark matter, expanding geometry or changed observed distance is present.

The primary cube has length12 and32 cells per side, step0.02, duration2.5 and source radius0.9. Six unit masses begin on a radius0.7 ring, resting or with tangential canonical momentum magnitude0.2 of either sign. The coupling is g=.08, kappa=.5, eta=0/.08 and omega=.2. Eighteen primary cases span source sense, vector coupling and probe type. Thirteen additional runs test time/grid/domain changes, a3D rotation and source-radius sensitivity.

The field energy averages forward/backward squared gradients and uses centered gradients in the drift term. Its canonical adjoint derivative, including coefficient self-sources, is implemented directly. The positive cellwise gradient/kinetic block follows from |beta|<alpha and the squared-gradient averaging bound. Particle energy samples alpha and beta using the same normalized spherical weights used for deposition; the derivative of weight normalization is retained in the force.

This regularization makes the particle cone depend on averaged alpha,beta, not the exact field cone at a point. A finite source profile remains nonlocal. The smallest recorded energy/cone residual does not remove that physical limitation. Continuum common-cone reasoning and local convexity from CC-1 remain distinct from proof of global nonlinear stability.

## Conservation and source funding

Maximum scaled drift in H+Q is {audit['maximum_energy_drift']:.6g}; the gate is1e-4. Q records energy removed by Pi_dot=-gamma Fdot, for which Qdot=integral gamma Fdot^2. Every integrator endpoint is saved. The maximum averaged-cone residual is {audit['maximum_cone_residual']:.6g}, consistent with roundoff. Independently reconstructing every final Hamiltonian agrees within {audit['maximum_independent_energy_error']:.6g}.

The omitted-self-source negative control gives an energy-ledger derivative error of {controls['negative_control']['omitted_self_source_rate']:.6g}, compared with {controls['negative_control']['correct_ledger_rate']:.6g} for the complete equations. This checks that energy accounting depends on the actual nonlinear terms, not merely the integrator tolerances.

For the rotating eta=.08 example without a probe, final field energy is {e['field']:.10g}, matter Hamiltonian energy loss is {m0-e['matter']:.10g}, and outlet energy is {e['absorbed']:.6g}. Funding includes field-dependent particle energy, not just mechanical kinetic energy. It does not establish a microscopic emission process, a galactic fuel supply or acceptable clock changes.

Maximum continuum momentum and angular-momentum diagnostic residuals over the whole campaign are {audit['maximum_momentum_residual']:.6g} and {audit['maximum_angular_residual']:.6g}. Field orbital terms and vector spin are included. The finite grid breaks exact continuous symmetry; these residuals are measured, not declared exact conservation identities.

## Refinement, domain and rotation

| Comparison | Maximum position difference | Final field-energy relative change | Declared result |
|---|---:|---:|---|
{chr(10).join(comparison_table)}

The time-only and spatial tests have position ceiling0.01 and field-energy ceiling5%. The enlarged-domain tests preserve cell spacing and use ceilings1e-4 and0.1%. Back-rotated trajectory comparisons use position ceiling0.01. Failure labels, if present, are retained rather than redefined. Two grid resolutions do not establish a continuum extrapolation.

The minimum conservative geometric clearance from the physical box boundary is {audit['minimum_boundary_clearance']:.6g}; the maximum measured boundary field amplitude is {audit['maximum_edge_amplitude']:.6g}. This uses the maximum measured characteristic speed, source support extent and run duration. Lattice tails need not have strict compact support. The larger-domain comparisons provide the direct domain-sensitivity test.

**Boundary qualification:** OB-1 separately showed that the original sponge fails after outgoing waves reach it. OB-2 provides a passing1D vacuum matched-wave prototype, not an integrated nonlinear3D boundary. CC-2 therefore supports the short tested interval only. Do not interpret this campaign as permission to use the same boundary for long-lived swirls. See [boundary report](boundary-report.md).

## Bending and physical smoothing

Final photon angles in milliradians from the initial path:

| Source | Scalar only | Scalar + vector | Change in bend magnitude |
|---|---:|---:|---:|
{chr(10).join(table)}

The sign depends on source rotation and the chosen off-axis light path. A transient vector contribution is not a universal extra scalar attraction. The massive probes are also evolved, but these paths are not circular stellar orbits. The fast source fixture is not a galactic-speed model: [CC-2S](slow-motion-report.md) quantifies the direct-current suppression at illustrative galactic speeds.

| Source radius in rotating photon case | Final field energy | Photon angle, mrad |
|---|---:|---:|
{chr(10).join(radius_table)}

Changing this radius changes the regulated physical source, not just mesh resolution. Sensitivity is retained rather than selecting the best-looking radius. Rotating a configuration and changing its smoothing radius answer different questions; a rotation pass does not make predictions independent of the regulator.

## What must change before observational claims

The current compact linear exterior fails the flat-rotation target and lacks the standard spatial-curvature optical contribution, as shown in [CC-2W](weak-field-report.md). Those findings are not repaired by the successful numerical evolution. [PF-0](phase-budget-report.md) tests an attributed reduced amplitude-feedback idea and finds substantial threshold sensitivity; it is not the full3D completion.

The next coupled candidate needs consistent spatial response, a derived mechanism for sustained outer support, source and angular-momentum budgets, and a validated long-duration boundary. Its parameters must then face held-out galaxies and clusters together, including local light/clock constraints. None of these broader requirements is marked achieved by CC-2. The [goal ledger](../../../research_plan/solution-goal-ledger.md) remains active (repository-root ledger path: research_plan/solution-goal-ledger.md).

## Attribution and reproduction

The Hamiltonian geometry, compact interpolation, finite differences and wave absorption methods are established mathematics; [CC-1](report.md) and [OB-2](boundary-report.md) record sources. These equations are a project candidate, not a claim of historical uniqueness, Einstein gravity or a generally covariant vector theory.

Protocol18dde4a and numerical sourceb69322f precede the original run. The archive pins all numerical sources and includes31 particle paths/momenta,31 final field states, every endpoint's energy/invariant/cone diagnostics, controls and declared comparisons. Full field movies are not stored. `python -B research_work/experiments/common_cone/audit_evolution.py` performs the independent audit; `report_evolution.py` renders this report/figure from the completed archive. The original driver refuses to overwrite evolution-v1. The historical repository-wide suite was not rerun or declared green.
'''
    # Link from experiments/common_cone climbs three levels to the repository root.
    (HERE/'evolution-report.md').write_text(report,encoding='utf8',newline='\n')
    print('Rendered evolution-report.md and evolution.png from audited complete evidence.')


if __name__=='__main__':main()

"""Assemble actual campaign results and verify every archived source/input digest."""
import ast,hashlib,json,subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import HERE,ROOT,read,save,hashes

def digest(data):return hashlib.sha256(data).hexdigest()
def main():
    directories=['small-v1','spatial1d-v1','spatial2d-v1','refinement-v1','explore-v1','mirror_audit-v1','readiness-v1']
    results={d:read(HERE/'evidence'/d/'results.json') for d in directories}
    sources=[];inputs=[];failures=[];manifests=[]
    for d in directories:
        folder=HERE/'evidence'/d;meta=read(folder/'manifest.json');manifests.append(meta)
        if meta.get('source_changed_during_run') or meta.get('excluded_modules_loaded') or meta.get('exception'):
            failures.append(dict(run=d,manifest=meta))
        for rel,expected in meta['source_sha256'].items():
            path=ROOT/rel
            gitpath=rel.replace('\\','/')
            committed=subprocess.check_output(['git','show',meta['commit']+':'+gitpath],cwd=ROOT)
            good=digest(path.read_bytes())==expected==digest(committed)
            sources.append(dict(run=d,path=rel,at_declared_commit=meta['commit'],passed=good))
            if not good:failures.append(dict(run=d,path=rel))
        for rel,expected in results[d].get('input_sha256',{}).items():
            good=digest((ROOT/rel).read_bytes())==expected
            inputs.append(dict(run=d,path=rel,passed=good))
            if not good:failures.append(dict(run=d,path=rel))
    imports=set()
    for path in HERE.glob('*.py'):
        tree=ast.parse(path.read_text(encoding='utf-8'),filename=str(path))
        compile(tree,str(path),'exec')
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):imports.update(a.name for a in node.names)
            elif isinstance(node,ast.ImportFrom) and node.module:imports.add(node.module)
    forbidden=sorted(imports&{'cl1','cl2','cl2_sources','astropy.cosmology'})
    if forbidden:failures.append(dict(forbidden=forbidden))
    small=results['small-v1'];one=results['spatial1d-v1'];two=results['spatial2d-v1']
    refined=results['refinement-v1'];wild=results['explore-v1'];mirror=results['mirror_audit-v1'];ready=results['readiness-v1']
    r={x['label']:x for x in refined['cases']};fine=r['primary'];lens=two['cases']['primary-256']
    homogeneous_count=len(small['rows'])+len(small['controls'])+2*len(wild['homogeneous'])
    spatial_count=len(one['cases'])+len(two['cases'])+len(refined['cases'])+len(wild['spatial'])+len(mirror['cases'])
    verification=dict(passed=not failures,source_hash_checks=len(sources),input_hash_checks=len(inputs),
                      declared_sources=sources,inputs=inputs,failures=failures,imports=sorted(imports),
                      run_directories=directories,homogeneous_cases=homogeneous_count,spatial_cases=spatial_count,
                      original_failed_gates_preserved={'spatial1d-v1':one['gates'],'spatial2d-v1':two['gates']},
                      evidence_sha256=hashes(sorted((HERE/'evidence').rglob('*.*'))))
    save(HERE/'publication-verification.json',verification)
    if failures:raise RuntimeError('Archive/source verification failed')
    # Figure: data only, no fabricated observational comparisons.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    raw=read(HERE/'evidence/refinement-v1/primary.json')
    history=np.array(raw['energy_history']);time=raw['time']
    axes[0,0].plot(time,history[:,0],label='Electromagnetic')
    axes[0,0].plot(time,history[:,1],label='Receiving wave')
    axes[0,0].plot(time,history.sum(axis=1),'k--',label='Total, including matter')
    axes[0,0].set(xlabel='Static coordinate time',ylabel='Counted energy',title='1D: energy transfers; the ledger closes')
    axes[0,0].legend(fontsize=8)
    labels=['Wavelength 1','Wavelength 2','Wavelength 4'];x=np.arange(3)
    axes[0,1].bar(x-.18,[r[k]['spectral_stretch']-1 for k in ['color-1','primary','color-4']],.36,label='Spectral stretch minus 1')
    axes[0,1].bar(x+.18,[r[k]['event_stretch']-1 for k in ['color-1','primary','color-4']],.36,label='Event stretch minus 1')
    axes[0,1].axhline(0,color='black',lw=.7)
    axes[0,1].set(xticks=x,xticklabels=labels,ylabel='Clock-corrected change',title='1D: timing and spectrum disagree')
    axes[0,1].legend(fontsize=8)
    fields=np.load(HERE/'evidence/spatial2d-v1/primary-final-field.npz')
    axis=fields['axis'];phi=fields['phi'];bound=np.max(abs(.2*phi))
    im=axes[1,0].pcolormesh(axis,axis,(.2*phi).T,cmap='RdBu_r',vmin=-bound,vmax=bound,shading='auto')
    for ray in lens['rays']:axes[1,0].plot(ray['x'],ray['y'],color='black',lw=.8,alpha=.7)
    axes[1,0].set(xlabel='x (fixture units)',ylabel='y (fixture units)',xlim=(-8,8),ylim=(-8,8),
                  title='2D: generated field and frozen rays')
    axes[1,0].set_aspect('equal');fig.colorbar(im,ax=axes[1,0],label='g phi at time 18',shrink=.8)
    names=['Baseline','Sign flip','Clock flip','Power 2','Power 3','Speed x2','Negative m2','Rest coupling']
    bends=[max(abs(ray['angle']) for ray in row['rays']) for row in wild['spatial']]
    axes[1,1].barh(names,bends,color=['#4878b0']*6+['#d48b36']*2)
    axes[1,1].set_xscale('symlog',linthresh=1e-4)
    axes[1,1].invert_yaxis()
    axes[1,1].set(xlabel='Maximum bend magnitude (radians)',title='Exploration: stronger bends are not joint passes')
    axes[1,1].text(.5,.44,'Orange: extra initial energy terms.\nSpeed x2 includes outward rays; powers 2/3 give zero.',
                   transform=axes[1,1].transAxes,ha='center',va='center',fontsize=8)
    fig.suptitle('CWC-1 | Fictional static-universe tests | No observational solution established',fontsize=14)
    fig.savefig(HERE/'overview.png',dpi=160);plt.close(fig)
    max_energy=max(x['relative_energy_error'] for x in small['rows']+small['controls'])
    max_time=max(x['timing_identity_error'] for x in small['rows']+small['controls'])
    report=f"""# CWC-1: coupled conversion, clocks and generated lensing
19 September 2026. Executed under the owner's static-universe, no-dark-matter
constraints. [Plan](protocol.md), [derivation](derivation.md),
[attribution](provenance.md), [follow-up attribution](provenance-follow-up.md),
[reproduction](README.md), [archive verification](publication-verification.json).

**Outcome:** The coupled Hamiltonian can transfer finite electromagnetic energy
into a receiving wave and produce attractive optical bending from a generated
profile. The tested family does not solve the joint timing, matter and cluster
problem. Sign/power/constant changes were actually run; none passes the whole
homogeneous scientific screen. No astronomical fit or historical novelty is claimed.

![Evidence overview](overview.png)

## Plan-to-execution record

| Stage | Executed evidence | Result |
|---|---|---|
| Attribution | Primary sources and prior project derivations; component claim register | Familiar ingredients credited; originality unestablished |
| Small systems | 27 parameter cases + 5 controls; Hamiltonian derivatives, reversal and missing-reciprocity negative control | Numerical pass; joint mechanism fail |
| Spatial conversion | 29 one-dimensional runs, including colors, intensities, direction, mass, wave speed and source-off controls | Original arrival-convergence gate failed; other numerical gates passed |
| Resolution follow-up | 6 runs at 4096, using original windows and thresholds | All four follow-up numerical gates pass; physical failures persist |
| Generated lensing | 6 two-dimensional runs, three grids, six signed ray impacts, refined rays, straight-field and weak-gradient controls | Inward bends; original mirror gate failed |
| Mirror diagnostic | Equal-action and reflected-source fixtures | Both symmetry diagnostics pass; original failed gate preserved |
| Unusual laws | 66 homogeneous laws + 66 source-off controls; 8 spatial screens and 8 Hamiltonian derivative checks | Numerical gates pass; 0/66 complete homogeneous mechanism passes |
| Observation readiness | 35 timing rows, 175 disks/3391 points, 6 lens systems, 6 Coma shear bins, 12 cluster density profiles | Input checks pass; all six prediction prerequisites remain unmet |

There are {homogeneous_count} homogeneous cases and {spatial_count} spatial runs,
plus the derivative, reversal, ray and input checks. The measurement inventory
is an executed readiness audit, **not** an observational model fit.
Historical jobs containing excluded geometry or halo comparisons were not run.

## What the tests establish
The starting receiving amplitude and momentum are zero in production runs.
Every reciprocal force follows the same Hamiltonian. Finite material internal
energy is included; material clocks are effective oscillators, not solved atoms.

Small-system energy error is at most {max_energy:.3g}, and the independent
arrival/clock identity error is at most {max_time:.3g}. The largest Hamiltonian
finite-difference discrepancy is {small['derivatives']['max_gradient_relative_error']:.3g}.
The deliberately omitted reciprocal source yields a
{small['derivatives']['missing_reciprocity_energy_defect']:.2%} energy defect and is rejected.

The 4096-point 1D fixture finishes with receiving energy
{fine['receiving_gain']:.6f} from initial electromagnetic energy 1.
Electromagnetic loss is {fine['em_energy_lost']:.6f}; counted material internal
energy loss is {fine['material_internal_energy_lost']:.6g}.
Maximum relative total-energy error is {fine['energy_relative_error']:.3g}.
This is classical field-energy transfer, not a proof of a microscopic photon
splitting process or a photon-number measurement.

The primary 2D snapshot bends all six rays inward. Maximum bend is
{lens['bend']['max_bend']:.7f} radians in dimensionless fixture units;
192-to-256 bend change is {two['spatial_bend_relative_change']:.3%},
ray-step refinement changes the answer by {lens['ray_refinement_relative']:.3g}
of the largest bend, and the weak-gradient comparison differs by
{lens['bend']['born_relative']:.3%}.
Total-energy error is {lens['energy_relative_error']:.3g}.
Only {lens['central_fraction_of_receiving']:.3%} of receiving energy remains
inside radius 3 at time 18. This is weak finite-duration central retention,
not a stable cluster reservoir. Bending convergence does not certify
convergence of every energy or field statistic.

These are frozen periodic 2D snapshots. There is no claim of evolved 3D
astrophysical light propagation, realistic cluster formation, actual source
geometry, physical mass normalization or universal gravity.

## What fails
The refined primary measured spectral stretch is {fine['spectral_stretch']:.6f}
(a blueshift under the defined diagnostic), while its event stretch is
{fine['event_stretch']:.6f}. Their relative mismatch is
{fine['timing_spectral_discrepancy']:.3%}, exceeding 1%.
The equal-energy carrier stretch spread is
{refined['color_stretch_fractional_spread']:.3%}, exceeding 1%.
The largest fixed-ruler, clock-normalized optical-speed change is
{fine['local_clock_speed_max_change']:.3%}, exceeding 0.1%.
Finite pulses are broad and distorted; these diagnostics are not precision
astronomical line measurements.

At a common position, the baseline force per inertial mass is proportional
to I/M. Doubling internal action at fixed inertial mass doubles the acceleration.
The b=0 case protects its clock by assumption but removes that material force.
The b=1 homogeneous closure keeps the fixed-ruler local speed constant but
cancels the measured homogeneous redshift. Changing the clock exponent is a
physical change requiring a material completion, not a coordinate relabeling.

For the generalized homogeneous F(phi) laws, the directly checked identity is

S_spectrum = S_events = exp[(1-b)(F_observed-F_emitted)].

This identity is for homogeneous propagation; the distorted spatial pulses
do not inherit it automatically. The baseline b=2 source-off run also
generates receiving energy from finite material internal energy. That is
counted and allowed as an energy source, but rules out calling it
photon-exclusive production. No broad photon-budget exclusion assumes a
universal cosmic age or a finite total illumination history.

## Numerical failures retained and investigated
The 1024-to-2048 arrival difference was {one['arrival_refinement_change']:.6f},
above 0.02. At 2048-to-4096 it is {refined['arrival_change']:.6f}, below the
unchanged 0.02 rule. Receiving-energy change is
{refined['receiving_change']:.6f} of initial EM energy.
The original failed directory remains [spatial1d-v1](evidence/spatial1d-v1/results.json).

The primary 2D mirror residual is {lens['bend']['mirror_relative']:.3%},
above 1%. Opposite material sources were deliberately assigned unequal actions.
Equalizing actions at the same total energy reduces the residual to
{mirror['cases'][0]['bend']['mirror_relative']:.3g}; reflecting the original
sources reproduces reflected rays to {mirror['reflected_covariance_error']:.3g}.
This diagnoses physical source asymmetry rather than discarding the failed
criterion. The original [spatial2d-v1](evidence/spatial2d-v1/results.json)
still has numerical_pass=false; the separate mirror audit passes.

## Deliberately unusual laws: all outcomes retained
| Change | Executed outcome | Interpretation |
|---|---|---|
| g to -g consistently | Identical ray observables; phi changes sign | Field relabeling for this symmetric potential, not a rescue |
| F=g phi^2, g phi^3, g phi^5 | Empty homogeneous state remains empty; 2D square/cube rays remain straight | F'(0)=0 prevents classical production without a new seed mechanism |
| b=-2,-1,0,0.5,1,2,4 | Some positive shifts; no joint homogeneous pass | Clock/force/source changes must all be counted |
| Receiving speed doubled | Smaller bending with mixed inward/outward rays | Changed propagation constant does not supply the missing completion |
| Negative squared mass plus bounded quartic | Stronger bend; 65.536 units of initial receiving potential energy in the 2D box | A counted initially unstable field reservoir, not an empty-energy receiving channel |
| Added material-rest coupling | Stronger bend and reduced composition contrast; material energy funds most generated waves | Additional matter coupling and source, not a photon-only derivation |

The material-rest 2D run starts with 4.003 material energy units and 1 EM unit.
It finishes with about 0.550 receiving energy, while EM loses only about 0.0116.
The larger optical signal cannot be attributed to photon losses alone.
Neither stronger-bending variant is promoted as a successful theory or a
precision result from its 128-point exploratory screen.

Negative kinetic-energy ghosts were rejected algebraically by the positive
energy requirement; no unbounded-energy model was used to manufacture a pass.
No constant, source term or sign was edited silently after a result.

## Attribution and next decision
Maxwell/Hamiltonian mechanics, dynamic-medium frequency conversion, scalar
couplings, scalar quartic potentials and symmetric splitting are established
mathematics. Photon/graviton mixing has its own established literature and is
not identical to this effective scalar receiving wave. The primary sources,
access limitations and permitted claim language are in the attribution files.
Code was written for this campaign; that does not confer originality on its equations.

Revisiting conversion and clocks together was useful: it gives an explicit,
conservative model with falsifiable cross-checks. This tested family should
not be promoted to a solution for galaxy clusters. A next model would need a
microscopic source interaction and a consistent matter/rod/clock completion
that predicts the same observed timing and material/lensing behavior without
an inserted unseen profile. The present results do not exclude all fictional
conversion or time alternatives.

## Checkpoints and integrity
Protocol/provenance checkpoint 7abb52b preceded execution; c11cbaa implemented
the Hamiltonian, f22aab2 declared 1D tests, f5d94ae published 1D evidence and
declared 2D tests, 1bfcf28 declared the unusual-law and resolution follow-ups,
and 13da431 published 2D/exploratory evidence and declared the final controls.
No force-push or rewrite of older PF5/RUT evidence was used.

Verification matches {len(sources)} source digests to the exact pre-run Git
commits and current files, and {len(inputs)} consumed-input digests to current
inputs. Every evidence file has a final digest. Failed scientific/numerical
gates are retained separately from successful process completion.
"""
    (HERE/'report.md').write_text(report,encoding='utf-8',newline='\n')
    print(dict(archive_pass=True,source_checks=len(sources),input_checks=len(inputs),
               homogeneous=homogeneous_count,spatial=spatial_count),flush=True)

if __name__=='__main__':main()

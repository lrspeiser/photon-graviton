"""Independent reconstruction of both RC-1 archives and a scoped report."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parent


def reconstructed(y,cfg):
    Q,P=y[:2];psi=y[2]+1j*y[3];other=y[4]+1j*y[5]
    population=np.abs(other)**2;norm=np.abs(psi)**2+population
    d=cfg['delta0']*np.tanh(Q);k=cfg['kappa0']*(1+cfg['A']*np.tanh(Q))
    environment=.5*(P**2+.04*Q**2)
    total=environment+(1+d)*np.abs(psi)**2+(1-d)*population+2*k*np.real(np.conj(psi)*other)
    return total,norm,population,environment


def main():
    checks=[]
    def check(name,value):checks.append(dict(name=name,passed=bool(value)))
    max_energy_reconstruction=0.;archives={}
    for version in ('reversible-conversion-v1','reversible-conversion-v2'):
        directory=ROOT/version;s=json.loads((directory/'summary.json').read_text());archives[version]=s
        manifest=json.loads((directory/'manifest.json').read_text());raw=np.load(directory/'trajectories.npz')
        for name,digest in manifest['sources'].items():
            check(version+' current '+name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
            committed=subprocess.check_output(['git','show',manifest['source_commit']+':research_work/experiments/common_cone/'+name],cwd=ROOT)
            check(version+' pinned '+name,hashlib.sha256(committed).hexdigest()==digest)
        check(version+' counts',len(s['runs'])==24 and len(s['controls'])==7 and len(s['negative_controls'])==12 and len(s['budgets'])==24)
        matches={}
        for row in s['controls']:
            cfg=row['config'];name=row['name'];y=raw[name];total,norm,pop,source=reconstructed(y,cfg)
            if name.startswith('no-conversion'):
                error=float(pop.max());matches[(cfg['delta0'],cfg['P0'])]=y
                check(version+' '+name,error<1e-12 and row['passed'])
            else:
                error=float(np.max(abs(pop-np.sin(cfg['kappa0']*raw['time'])**2)))
                check(version+' '+name,abs(error-row['population_error'])<1e-12 and error<1e-8 and row['passed'])
        for row in s['runs']:
            name=row['name'];cfg=row['config'];y=raw[name];total,norm,pop,source=reconstructed(y,cfg)
            max_energy_reconstruction=max(max_energy_reconstruction,abs(total[-1]-row['final_energy']))
            ed=float(np.max(abs(total-total[0]))/abs(total[0]));nd=float(np.max(abs(norm-1)))
            check(version+' '+name+' H/N',abs(ed-row['relative_energy_drift'])<1e-12 and abs(nd-row['norm_drift'])<1e-12 and abs(total[-1]-row['final_energy'])<1e-12)
            qdiff=float(np.max(abs(y[0]-matches[(cfg['delta0'],cfg['P0'])][0])))
            check(version+' '+name+' response',abs(qdiff-row['max_environment_difference_from_no_conversion'])<1e-12)
            global_bound=1-np.sqrt(cfg['delta0']**2+(cfg['kappa0']*(1+abs(cfg['A'])))**2)
            check(version+' '+name+' positivity',global_bound>0)
            if version.endswith('v2'):
                back=float(np.max(abs(raw[name+'-reversed']-y[:,0])))
                check(version+' '+name+' reversal',abs(back-row['reversal_error'])<1e-12 and back<1e-7)
            expected=ed<1e-8 and nd<1e-8 and row['reversal_error']<1e-7 and row['symmetric_no_force_identity']
            check(version+' '+name+' honest gate',expected==row['passed'])
        for row in s['negative_controls']:
            total,_,_,_=reconstructed(raw[row['name']],row['config'])
            error=float(np.max(abs(total-total[0]))/abs(total[0]))
            check(version+' '+row['name'],abs(error-row['relative_energy_drift'])<1e-12 and error>1e-5 and row['detected'])
        for i,row in enumerate(s['budgets']):
            required=(200000.**2*(6e20-6e18)/6.67430e-11)/(2e41*row['available_fraction'])
            check(version+' budget '+str(i),abs(required/row['required_effective_coefficient']-1)<1e-12 and row['equal_response_active_coefficient']==1+row['s'] and row['equal_response_budget_sufficient']==(1+row['s']>=required))
        expected=all(r['passed'] for r in s['runs']+s['controls']) and all(r['detected'] for r in s['negative_controls'])
        check(version+' aggregate gate',expected==s['passed'])
    first=archives['reversible-conversion-v1'];second=archives['reversible-conversion-v2']
    check('original failure retained',not first['passed'] and sum(not r['passed'] for r in first['runs'])==1)
    check('refinement succeeds',second['passed'])
    output=dict(passed=all(c['passed'] for c in checks),count=len(checks),max_energy_reconstruction_error=max_energy_reconstruction,checks=checks)
    (ROOT/'reversible-conversion-audit.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    if not output['passed']:raise RuntimeError('RC-1 independent audit failed')
    rows=second['runs'];symmetric=[r for r in rows if r['config']['delta0']==0];detuned=[r for r in rows if r['config']['delta0']!=0]
    report=f'''# RC-1: reversible conversion and environmental response

20 September 2026. A new effective interaction test in a nonexpanding setting, with no dark-matter component. Judged by its own Hamiltonian and conservation identities, not agreement with an older gravity formula. No observational fit or ordinary spin-2 identification is claimed.

## Direct answer to the proposed mechanism

Reversible conversion can be implemented without creating energy. In this local toy model, changing the relative response of the two channels also changes the environment's motion through an explicitly reciprocal force. Simply converting between identically responding channels does not give stronger gravity. Conversion, temporary binding, and subsequent travel on different paths are distinct steps; only local conversion and its environmental reaction are tested here.

Use the Hermitian matrix M(Q)=[[1+delta(Q),kappa(Q)],[kappa(Q),1-delta(Q)]], amplitudes psi=(a,b), and

H = P^2/2 + .2^2 Q^2/2 + psi^dagger M(Q) psi.

i psi_dot=M psi; Qdot=P; Pdot=-.2^2 Q-psi^dagger M'(Q) psi.

Here delta=.0 or .3 times tanh Q and kappa=kappa0(1+A tanh Q). The same interaction supplies conversion in both directions and the environmental force. A local state can influence the rate; no desired orbit or lens map is provided to it. The environment is one oscillator, not a completed spatial matter/field system. Its generalized force must not be advertised as a verified 3D momentum or angular-momentum budget.

## Executed evidence

- 100 finite-difference source-force/positive-matrix controls pass.
- First archive: 23/24 main trajectories meet every declared gate; run-19 fails backward-state accuracy (5.24614e-7 versus 1e-7). The failure is preserved.
- With unchanged equations and gates, stricter tolerances give 24/24 passes; all seven analytic/no-conversion controls and twelve deliberately missing-reaction controls pass their intended tests.
- Maximum relative energy drift {max(r['relative_energy_drift'] for r in rows):.6g}; norm drift {max(r['norm_drift'] for r in rows):.6g}; backward-state error {max(r['reversal_error'] for r in rows):.6g}.
- Highest sampled receiving-mode population {max(r['peak_receiving_population'] for r in rows):.12g}. Constant-rate tests recover sin^2(kappa t), including return conversion. This is normalized mode population, not extra energy created or a guaranteed microscopic conversion probability.
- Identical-channel cases change the environmental coordinate from their no-conversion controls by at most {max(r['max_environment_difference_from_no_conversion'] for r in symmetric):.6g}. Their interaction-force expectation is zero for the declared pure initial state.
- Different-channel cases change that coordinate by as much as {max(r['max_environment_difference_from_no_conversion'] for r in detuned):.6g} model units. These are response changes, not stable orbits or improved fits. Positive local mode energy does not prove full nonlinear stability.
- Omitting the environmental reaction produces relative energy errors from {min(r['relative_energy_drift'] for r in second['negative_controls']):.6g} to {max(r['relative_energy_drift'] for r in second['negative_controls']):.6g}. Counting only the two mode energies and ignoring signed interaction energy is also invalid.
- Independent raw-state reconstruction: {output['count']} audit checks pass; final-energy reconstruction difference <= {max_energy_reconstruction:.6g}.

## Enough gravitational response?

For a bookkeeping diagnostic with energy fraction f and derived active coefficients chi_gamma,chi_g:

M_active = [(1-f)chi_gamma+f chi_g] E_stored/c^2,

before any additional interaction source is included. In SR-1's weak branch, both coefficients are 1+s. The photon Hamiltonian derivative and the averaged wave-energy derivative therefore give no leading gain from conversion alone. This is a conditional result about SR-1, not a prohibition on a different invented interaction.

For the prior illustrative 200 km/s, R=6e20 m, M=2e41 kg budget, keeping only epsilon=.001 of source rest energy available requires an effective coefficient about 1780 in that same equivalent-energy interpretation, compared with SR-1's 2. The new theory could instead derive a different constitutive/interaction force; then its own force and energy relation must replace this budget. Choosing 1780 by hand does not explain observations. The budget's f is an energy fraction and must not be confused with the mode population when diagonal energies and interaction energy differ.

## Attachment and separate outgoing paths

The owner's extension is retained as the next spatial mechanism: a companion may couple to matter/radiation, travel in a bound or coherent state, and detach with a different outgoing direction. A complete implementation must derive binding energy, capture/release dynamics, the outgoing momenta and field/source reaction from one interaction. Deterministic or coherent release can in principle yield repeatable bending; random scatter is not assumed just because detachment occurs. Image-width and spectral predictions must decide which behavior the law actually produces.

Within any retained standard massless energy-momentum kinematics, an isolated massless particle cannot split into two non-collinear positive-energy massless particles without another participant: (p1+p2)^2=2 E1 E2(1-cos theta)/c^2. Matter, a finite-mass bound state or the surrounding field can supply the required balance. This is a stated kinematic assumption, not an old gravity law used as a veto. If the candidate changes it, specify and test the replacement conservation/propagation rules. No such bound state or spatial splitting dynamics is demonstrated by RC-1.

The next model must derive the field sourced by matter, spatial following and release, common signal propagation, an affordable long-term budget and both stellar and photon trajectories. Then test direct observations with shared constants. The current oscillator calculation supplies none of the missing galaxy/cluster profiles.

## Attribution and reproduction

Hermitian mode mixing and Hamiltonian reaction terms are established mathematical tools. Photon/graviton conversion has Gertsenshtein precedent; [Palessandro and Rothman](https://arxiv.org/abs/2301.02072) provide a derivation. Their magnetic-field mechanism and rates are not imported here. The local Q-dependent specialization is our declared effective candidate, not a claim of historical uniqueness. Existing CWC-1 radiation-to-scalar tests remain separate. [Weinberg's universal-coupling result](https://doi.org/10.1103/PhysRev.135.B1049) applies under stated framework assumptions; it is attribution/context, not a rejection gate for a different invented framework.

Use a separate checkout at b58ce4e, before the evidence archives are committed, and run reversible_conversion.py followed by refine_reversible_conversion.py. The first command intentionally exits unsuccessfully after preserving its complete failed first archive; proceed to the declared refinement. Run audit_reversible_conversion.py from the completed campaign against either reproduced or checked-in evidence. Source hashes and commits are preserved in both manifests.
'''
    (ROOT/'reversible-conversion-report.md').write_text(report,encoding='utf-8')
    print(json.dumps({k:v for k,v in output.items() if k!='checks'}))


if __name__=='__main__':main()

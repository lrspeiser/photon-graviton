"""Build the report and scientific plots from immutable evidence, not new fits."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
VERSIONS=['e1-v2','e2-v2','e3-v1','e4-v1','e5-v1','e3-colored-v1','audit-v2']

def load(name):
    return json.loads((HERE/'evidence'/name/'results.json').read_text())

def main():
    e1,e2,e3,e4,e5,col,audit=[load(n) for n in VERSIONS]
    assert all(x['numerical_pass'] for x in [e1,e2,e3,e4,e5,col,audit]), 'Report template requires passing final numerical gates'
    assert not any(x['mechanism_pass'] for x in [e1,e2,e3,e4,e5,col]), 'Update narrative for changed mechanism outcomes'
    fig,axes=plt.subplots(2,3,figsize=(15,9),constrained_layout=True)
    colors=['#216b85','#bd672d','#6950a1']
    ax=axes[0,0]
    for row,color in zip([x for x in e1['comparisons'] if x['softening']=='0.05'],colors):
        ax.plot(e1['radii_over_scale'],row['ratio'],'o-',label=f"height/Rd={row['height']}",color=color)
    ax.axhline(3,color='black',ls=':',label='Required contrast')
    ax.set(title='1. Angular activation fails the contrast',xlabel='Radius / characteristic scale',ylabel='Median cluster / disk activation')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[0,1]
    row=e2['primary'][1]['fits'][0]
    ax.errorbar(e2['radius_Mpc'][:-2],e2['observed'],yerr=e2['error'],fmt='o',color='black',label='Exposed Coma bins')
    for key,color,label in [('baryon_fit',colors[0],'Ordinary matter'),('fit',colors[1],'Maturing transport'),('old_fit',colors[2],'Old written response')]:
        ax.plot(e2['radius_Mpc'][:-2],row[key]['prediction'],color=color,label=f"{label}: chi2={row[key]['chi2']:.2f}")
    ax.set(xscale='log',title='2. Longer range does not fix the shear shape',xlabel='Adopted fixed radius [Mpc]',ylabel='Tangential shear')
    ax.legend(fontsize=7);ax.grid(alpha=.2)
    ax=axes[0,2]
    for data,shift,label,marker in [(e3,1,'White-time forcing','o'),(col,1.12,'Finite correlation','s')]:
        for avg in [.1,1,10]:
            rr=[r['radial_std']*100 for r in data['orbits'] if r['average_periods']==avg]
            ax.scatter(np.full(len(rr),avg*shift),rr,marker=marker,color=colors[0 if shift==1 else 1],label=label if avg==.1 else None)
    ax.axhline(5,color='black',ls=':')
    ax.set(xscale='log',yscale='log',title='3. Long averaging reduces orbit noise',xlabel='Averaging / orbital period',ylabel='Radial standard deviation [% of launch radius]')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1,0]
    for eta,color in zip([0,1,5],colors):
        row=next(r for r in e4['rows'] if r['eta']==eta and r['grid']==1024)
        snapshot=row['snapshots'][-1]
        ax.plot(snapshot['x'],snapshot['curvature'],color=color,label=f"erasure eta={eta}")
    ax.axvline(-2.2,color='black',ls=':');ax.axvline(2.2,color='black',ls=':',label='Passing lights')
    ax.set(xlim=(-3,3),title='4. Optical peaks shift; gas-central gate fails',xlabel='Static fixture coordinate',ylabel='Optical curvature [conditional units]')
    ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1,1]
    error=np.array(e5['selected']['fractional_bend_error'])*100
    labels=[r['name'].split('-')[0].split('+')[0] for r in e5['lenses']]
    ax.bar(labels,error,color=[colors[0]]*5+[colors[1]])
    ax.axhspan(-3,3,color='green',alpha=.12,label='Required +/-3%')
    ax.tick_params(axis='x',rotation=30)
    ax.set(title='5. Shared tensor fails lens transfer',ylabel='Predicted bend error [%]')
    ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
    ax=axes[1,2]
    ax.bar(labels,[l['motion']['chi2_per_bin'] for l in e5['lenses']],color=colors[2])
    ax.axhline(3,color='black',ls=':',label='Required <=3')
    ax.tick_params(axis='x',rotation=30)
    ax.set(yscale='log',title='5. Ordinary stellar motions also fail',ylabel='Motion chi2 / measured bin')
    ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
    fig.suptitle('PF5: static space, ordinary sources, no dark matter | numerical verification is not a successful theory',fontsize=14)
    fig.savefig(HERE/'overview.png',dpi=165)
    plt.close(fig)
    report=['# PF5 results: five photon-response hypotheses',
            '',
            'All five primary experiments and the finite-correlation extension were executed. Their final declared numerical gates pass; **none establishes a complete cluster-lensing solution**. Long-averaged fluctuation probes do show partial success. No dark-matter source or expanding-universe dynamics was used.',
            '',
            '![Six-panel results overview](overview.png)',
            '',
            '| Experiment | Executed coverage | Numerical result | Mechanism outcome |',
            '|---|---|---|---|',
            '| Angular activation | 175 disk light profiles, 12 cluster proxies, 3 heights, 3 softenings, nested source sampling | Refined pass; original sampling failure retained | Required threefold cluster/disk boost fails |',
            '| Maturing transport | 27 parameter settings, 6 primary sampling runs, ballistic/nonmaturing controls, 2 emissivity sensitivities | Pass | Longer range, worse primary Coma shape than ordinary matter |',
            '| Fluctuation amplitude | Independent covariance solves; 12 white-time and 12 finite-correlation orbit runs; 3 correlation-time covariance checks | Pass | Square-root scaling; only long averaging consistently quiet in these fixtures |',
            '| Shock erasure | 9 collision runs (3 grids x 3 strengths), source-off, uniform-flow and analytic-decay controls | Pass | Modest peak displacement; gas-central condition fails |',
            '| Directional optical strain | 4 reciprocal-mode energy controls, Hamiltonian rays, 5 shared lengths, 6 lenses and refined spatial solves | Pass | Lens bending and stellar motions fail |',
            '',
            '## Angular activation',
            '',
            f"The median nested-sampling score change fell from 0.05109 (failed 0.05 limit) to {e1['median_resolution_change']:.5f}. Angular isotropy/beam/plane, rotation, luminosity and subdivision controls pass. At one characteristic radius the cluster/disk median ratios are "+', '.join(f"{r['ratio'][2]:.3f}" for r in e1['comparisons'] if r['softening']=='0.05')+'. The required ratio is at least 3 across all declared radii and thicknesses. Full quantile overlap and softening sensitivity are retained in the JSON. These are disk-only light and thermal-emissivity proxies, not bolometric radiation reconstructions.',
            '',
            '## Maturing transport',
            '',
            f"The frozen primary gives Coma chi2={row_chi(e2,'fit'):.5f}, versus ordinary matter {row_chi(e2,'baryon_fit'):.5f} and the older pressure-written response {row_chi(e2,'old_fit'):.5f}. These are six exposed bins with one nonnegative geometry nuisance (five nominal degrees of freedom); they are weak-shear shape screens, not absolute lensing predictions. Both ordinary-mass brackets have the same shape here because their adopted gas and stellar profiles are proportional.",
            '',
            f"The primary contribution at 5 Mpc is {e2['primary'][1]['fits'][0]['extension_ratio']:.3f} times the old finite-footprint contribution. The maximum packet-doubling change is {max(e2['relative_resolution_changes'])*100:.2f}% on the normalized maximum response scale (limit 15%). All 27 settings are published without replacing the declared primary with a scan winner. A positive optical coefficient was calibrated at 1 Mpc from the old pressure response before scoring Coma; it is not an independently funded microscopic interaction. Packet propagation is causal; the Gaussian optical readout remains a conditional nonlocal rule.",
            '',
            '## Fluctuations and orbit stability',
            '',
            f"The white-time stationary covariance solves agree to {e3['checks']['relative_covariance_error']:.3g}; mean injection/damping agree to {e3['checks']['relative_energy_error']:.3g}. The amplitude exponent is {e3['checks']['scaling_exponent']:.12f}. The finite-correlation primary solves agree to {col['checks']['relative_covariance_error']:.3g}; the other declared correlation times are retained. The square-root law is conditional on linear covariance scaling and a square-root readout, not a demonstrated mass-speed relation for real galaxies.",
            '',
            '| Time forcing | Averaging/orbit | Quiet primary seeds /3 | Refined seed quiet |',
            '|---|---|---|---|']
    for label,data in [('White',e3),('Finite correlation 0.2',col)]:
        for avg in [.1,1,10]:
            count=sum(r['quiet'] for r in data['orbits'] if r['average_periods']==avg)
            fine=next(r for r in data['refined'] if r['average_periods']==avg)
            report.append(f"| {label} | {avg} | {count} | {fine['quiet']} |")
    report += ['',
            'Each orbit is measured for 20 periods after an initially empty field burns in for at least five averaging times. Quiet means radial standard deviation and median-radius drift each below 5%, with angular-momentum error below 1e-5. The saved aggregate mechanism flag requires all averaging-time variants to be quiet and is therefore false; this does not reject the successful long-averaging fixtures. Refinements compare statistical behavior, not matched Brownian trajectories. These are negligible-mass probes of a radial fixture. The field ledger in expectation and its explicit split-integrator residual are saved; signed stochastic exchanges do not establish a finite positive photon reservoir or matter backreaction.',
            '',
            '## Gas collision and erasure',
            '',
            f"The largest relative gas+field+fuel+boundary energy error is {max(r['relative_energy_error'] for r in e4['rows']):.3g}; no positivity floor adds energy. The source-off field stays zero and the uniform-flow entropy and constant-source decay controls pass.",
            '',
            '| Erasure strength | Positive optical peak | Positive gas peak | Distance from light at 2.2 | Noncompression entropy residual |',
            '|---|---|---|---|---|']
    for r in e4['rows']:
        if r['grid']==1024:
            report.append(f"| {r['eta']} | {r['optical_peaks'][1]:.4f} | {r['gas_peaks'][1]:.4f} | {r['star_peak_distance']:.4f} | {r['noncompression_residual_fraction']:.1%} |")
    report += ['',
            'Erasure moves optical peaks slightly closer to the passing lights, but gas peaks fail the predeclared central-region bound 0.75. About one third of the measured entropy residual lies outside compression, so this numerical residual cannot be treated as pure physical shock entropy. The simulation is 1D with prescribed light tracks, local field storage and a passive optical readout; no observed merger or full spatial energy-momentum closure is claimed.',
            '',
            '## Tensor optics and lens galaxies',
            '',
            f"The homogeneous reciprocal mode conserves photon+field+heat energy to {max(x['relative_energy_error'] for x in e5['homogeneous']):.3g}. Isotropic, source-off and zero-coupling controls vanish. The optical metric remains positive. The explicit maximum-step ray refinement agrees to {audit['ray_timestep_relative_change']:.3g}; agreement with the analytic weak Gaussian ray is {e5['analytic_ray_error']:.3g}. Spatial lens refinement changes the unit tensor bend by at most {max(x['relative_bend_change'] for x in e5['refined']):.3%}.",
            '',
            f"The five-lens training screen selects length {e5['selected']['length_kpc']} kpc and coupling {e5['selected']['coupling']:.6g}; maximum tensor eigenvalue magnitude is {e5['selected']['max_abs_tensor']:.3g}, below 1e-3. These settings are carried unchanged to exposed J1630. The largest primary bend error is {max(abs(x) for x in e5['selected']['fractional_bend_error']):.2%}, exceeding 3%.",
            '',
            '| Lens | Bend discrepancy | Stellar-motion chi2/bin |',
            '|---|---|---|']
    for l,err in zip(e5['lenses'],e5['selected']['fractional_bend_error']):
        report.append(f"| {l['name']} | {err:+.2%} | {l['motion']['chi2_per_bin']:.2f} |")
    report += ['',
            'The stellar mass is fitted only within the adopted Chabrier–Salpeter interval, with constant orbital anisotropy allowed in [-2,0.45]. Optical fitting cannot repair the severe ordinary-stellar motion failure. Light anisotropy uses stellar-mass-scaled luminosity proxies, not measured bolometric input. The reciprocal homogeneous test and spatial phenomenological response have not been derived as one closed theory.',
            '',
            '## Independent checks and reproducibility',
            '',
            f"An independent finite-packet survival test agrees with its analytic prediction within {audit['survival']['z']:.3f} sampling standard errors. Ballistic and constant-scattering displacement checks are within "+', '.join(f"{m['z']:.3f}" for m in audit['moments'])+f" standard errors. Point-source projected curvature error is {audit['gaussian_curvature_relative_error']:.3g}. The isolated CL-F1 analytic and independent 3D optics checks also pass. No historical joint suite containing excluded comparisons was run.",
            '',
            'All numerical thresholds were declared before their respective runs. The first E1 sampling failure and the E2 array-shape execution failure remain under evidence/e1-v1 and evidence/e2-v1. Refinement retained the original scientific thresholds. Manifests include source/input hashes and starting commits. An additional explicit maximum-step refinement avoids the identical-step limitation of the original tolerance comparison. Archive verification resolves all input/source hashes in Git history; ancillary report snapshots may be committed after a run begins. The latest runner inventories active repository modules and rejects the excluded CL1/CL2 source loaders. Historical data files contain unused model branches; only the documented ordinary-source measurements and static-registry values enter these calculations.',
            '',
            'Source and run instructions: [README](README.md). Original [protocol](protocol.md), [E1/E2 implementation](implementation-notes.md), [E3/E4 integrators](integrators-e3-e4.md), [E5 integrator](integrator-e5.md), [extra controls](extended-checks.md).',
            '',
            '## What remains unresolved',
            '',
            'A common photon interaction must still supply the response normalization, conserve energy and momentum with moving matter, predict the bolometric radiation history, and reproduce both cluster shear and galaxy/lens dynamics with shared parameters. Present evidence supports a few toy mechanisms and rejects the complete declared screens. Novel mathematical ingredients are proposals for this fictional universe; historical uniqueness has not been established.',
            '',
            'Evidence:', '']
    for name in VERSIONS:report.append(f"- [{name}](evidence/{name}/results.json) ([manifest](evidence/{name}/manifest.json))")
    (HERE/'report.md').write_text('\n'.join(report)+'\n',encoding='utf-8',newline='\n')
    print('Generated report.md and overview.png from immutable evidence')

def row_chi(data,key):
    return data['primary'][1]['fits'][0][key]['chi2']

if __name__=='__main__':main()

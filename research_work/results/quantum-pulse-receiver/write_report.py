from pathlib import Path
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'results.json').read_text())
fig,axes=plt.subplots(1,2,figsize=(11,4.8),sharex=True,sharey=True)
for ax,f in zip(axes,d['families']):
    for width,color in [(.02,'#d95f02'),(.1,'#7570b3'),(4.,'#1b9e77')]:
        curve=next(x for x in f['curves'] if x['receiver_width']==width)
        ax.plot(curve['time'],curve['profile'],label=f'Receiver energy width {width:g}',color=color)
    ax.set_xlim(-65,80);ax.set_ylim(0,.06);ax.set_title(f"Chosen frequency factor 1/{f['stretch']:g}")
    ax.set_xlabel('Fourier time coordinate (arbitrary units)');ax.grid(alpha=.2)
axes[0].set_ylabel('Normalized photon-number profile');axes[1].legend(fontsize=8)
fig.suptitle('Receiver coherence controls the pulse profile',fontsize=14)
fig.text(.5,.012,'Synthetic finite spectral channel; these are not causal detector-arrival calculations.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.04,1,.94]);fig.savefig(HERE/'profiles.png',dpi=170);plt.close(fig)
rows=[x for f in d['families'] for x in f['rows']]
lines=['# A quantum energy/momentum transfer can preserve a stretched pulse profile', '',
'A stipulated finite quantum transformation now connects a chosen photon-energy reduction to a stretched Fourier pulse profile while preserving total energy and spatial momentum. A sufficiently coherent initial receiver produces a profile close to ideal frequency dilation; removing receiver coherence while keeping its energy probabilities the same removes the temporal localization in this finite model.', '',
'This is a conditional spectral-channel result, not an astronomical propagation mechanism. The frequency map and special receiver spectrum are postulated. There is no spacetime-local interaction, causal signal-arrival calculation or escaping companion stream. Lost photon energy remains in the receiver. It therefore does not yet satisfy the proposed void transport rule or establish supernova timing.', '',
'## What is new relative to the earlier receiver check', '',
'The earlier finite receiver diagnostic tested two-frequency beats and total energy. This calculation uses whole finite-band pulse profiles and a positive-timelike receiver spectrum, checking momentum in every swapped basis state as well as energy. It computes how finite receiver coherence broadens the profile. These are additional conditional checks; they do not replace the missing spatial and repeated-use dynamics with a completed model.', '',
'## Formulas and provenance', '',
'**Known quantum and relativistic mathematics:** a permutation of orthogonal basis states is unitary; swaps within the same total-energy and total-momentum sector preserve those quantities. The receiving basis has E_R-P_R=U=100 in c=1 units, so its invariant mass satisfies m_R^2=U(2E_R-U)>0. This is a prescribed spectrum of internally different mass/momentum states, not the dispersion relation of an ordinary fixed-mass particle or a derived field Hamiltonian.', '',
'**Hypothetical finite-band interaction:** within designated disjoint photon bands, swap', '',
'    |E, E_R, P_R> <-> |E/S, E_R+(1-1/S)E, P_R+(1-1/S)E>.', '',
'The reverse swap is present; this is not a universally irreversible down-converter. Input/output labels are one-to-one and disjoint, with identity on other unpaired states. The tested factors are S=2 and S=1.5. The rule is selected, not predicted from gravity, void density, distance or an underlying time constant. It is not a claim of novel fundamental physics.', '',
'**Known Fourier scaling, conditionally retained by this receiver:** let input temporal amplitudes use the convention sum_E a_E exp(-iEt), with hbar=1. An input center t_b has spectral phase exp(i E t_b). Parameterize the initial receiver reference by exp(i E_R t_R). For an approximately translation-coherent receiving state, the outgoing center is', '',
'    t_profile,out = S t_b - (S-1)t_R.', '',
'Thus the same selected energy map yields pulse-center separations scaled by S when the initial receiver reference is held fixed across separate trials. Shifting both input and receiver references together by delta instead shifts the output by delta; the total system retains time-translation symmetry. The connection between energy conservation, receiving-state coherence and that symmetry is established quantum theory; see [Lostaglio et al., Physical Review X (2015)](https://arxiv.org/abs/1410.4572).', '',
'**Known Gaussian overlap/Fourier algebra applied to the postulate:** for a receiver energy wavefunction with probability standard deviation sigma_R, conditional receiver states overlap approximately as exp[-(deltaE-deltaE_prime)^2/(8 sigma_R^2)]. This suppresses photon coherence when the receiver can distinguish different energy transfers. With input photon energy width sigma_E,', '',
'    sigma_t,out^2 = [S/(2 sigma_E)]^2 + [(S-1)/(2 sigma_R)]^2.   [hbar=1]', '',
'The first term is ideal pulse stretching; the second is extra temporal blurring. Restore hbar in both numerators for dimensional units. The finite-grid calculation uses actual receiver overlaps, not this Gaussian approximation as its output. A broad incoherent energy distribution is not interchangeable with a coherent superposition.', '',
'## Executed results', '',
'The S=2 input band is 4–6 energy units and the output band 2–3; the S=1.5 bands are 4.5–6 and 3–4. Receiver mean energy initially is 100 units, with mean momentum approximately zero and positive invariant mass at every basis point. The initial coherent resource is included in the account; it is not created from an empty receiver.', '',
'| S | Receiver energy width | Output center for input center 4 | Output width | Fidelity to the ideal stretched pure pulse |',
'|---:|---:|---:|---:|---:|']
for f in d['families']:
    for r in f['rows']:
        if r['birth_phase_time']==4 and r['receiver_reference_time']==0:
            lines.append(f"| {f['stretch']:g} | {r['receiver_width']:g} | {r['output_profile_center']:.4f} | {r['output_profile_sd']:.4f} | {r['ideal_pulse_fidelity']:.6f} |")
lines += ['', 'For S=2, mean photon energy falls from 5 to 2.5 while receiver mean energy rises from 100 to 102.5 and its mean momentum rises from zero to 2.5. For S=1.5 the photon changes from 5.25 to 3.5 and the receiver gains 1.75 energy and momentum units. These are receiver gains, not additional energy that can also be counted as an escaping companion stream.', '',
'![Pulse profiles versus receiver coherence](profiles.png)', '',
f"Forty-eight coherent-state trials check three input centers, two receiver phase references, four receiver widths and two chosen dilation factors. Maximum output-center error is {max(abs(x['output_profile_center']-x['predicted_center']) for x in rows):.3g}; maximum fractional width error relative to the Gaussian expression is {max(abs(x['output_profile_sd']/x['gaussian_predicted_sd']-1) for x in rows):.3g}. Both direct joint-state partial-trace controls agree with the overlap construction. Density-matrix normalization/positivity and direct Fourier sums pass. Every basis-state swap preserves energy/momentum; the largest mean-energy residual is below 3e-14.", '',
'Dephasing the receiver produces a diagonal outgoing spectral density matrix and a uniform Fourier-time profile over its discrete recurrence period. This describes lost spectral phase information in the finite channel, not a simulation of photons physically appearing uniformly throughout cosmic time. Separate input-center trials initialize the receiver anew; they are not a continuous sequence through the same evolving receiver. The earlier repeated-use toy does not supply that missing spatial evolution here.', '',
'## Why this is not yet a causal propagation solution', '',
'Fourier time is the coordinate of an ideal photon-number profile at an unspecified port. The model has no physical input/output locations or switching dynamics. Energy/momentum conservation of asymptotic basis states does not establish local stress-energy flow during the interaction, Lorentz-covariant dynamics, or a causal arrival-time kernel.', '',
'For illustration, a naive broadband extension of the affine map with added delay D would give t_out=D+S*t_in-(S-1)*t_R. A causal retiming map must satisfy t_out>=t_in, requiring t_in>=t_R-D/(S-1) for S>1. An unqualified affine extension over all past times fails that condition. A finite-history implementation and its boundary behavior must be derived. This is not a claim that an advanced pulse peak alone proves acausality, nor a causality proof for the finite-band channel.', '',
'Supernova emission is an ensemble of photons with evolving spectra, not one coherent optical pulse lasting days. Its measured event stretching therefore still requires a spatially causal, repeatedly used receiving system and an explicit source/detector response. The near-ideal profiles here do not supply that astronomical prediction.', '',
'## Next physical requirement', '',
'Construct a local finite-history interaction whose phase evolution produces the frequency/time relation, whose receiving state evolves without hidden resets, and whose lost photon energy exits as companions with accounted momentum. The present receiver temporarily/ultimately retains that energy and has a specially prepared spectrum; it is not adopted as a permanent void deposit or as an identified graviton.', '',
'No redshift coefficient, galaxy parameter, observational likelihood or held-out score was fitted. Receiver preparation, locality, clock standards, outgoing companion transport, capture/support, lensing, stable stellar-population predictions and the deferred total photon-supply budget remain open. Run run.py and write_report.py in this directory to reproduce the finite-channel results and figure.']
(HERE/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
dest=Path('C:/Users/henry/Documents/Codex/2026-09-09/cr/outputs/quantum-pulse-receiver');dest.mkdir(parents=True,exist_ok=True)
for name in ['report.md','results.json','profiles.png']:shutil.copy2(HERE/name,dest/name)
print(dest)

# Matter, clocks and the autonomous receiver

The homogeneous universal-clock completion conserves its Hamiltonian energy but produces **zero measured redshift and zero event-duration stretch** between identical comoving local laboratories. Including ordinary matter also changes where the receiver's energy comes from: it cannot all be credited to photons. This closes a specific gap between earlier photon-only backreaction and prescribed universal-clock examples. It is a synthetic consistency test, not an astronomical fit or a rejection of all companion theories.

In everyday terms, changing the signal and the measuring clock together can leave the measured signal unchanged. An energy ledger written against an external reference clock does not by itself show that a detector sees red light.

## What was already established, and what is new in this calculation

The [earlier backreaction calculation](../companion-backreaction/report.md) already evolves a positive receiver field under a closed photon–field Hamiltonian. [Distributed modes](../distributed-modes/report.md) already explore finite three-dimensional source distributions. [Spherical propagation](../spherical-propagation/report.md) already shows localized free-field settling. None should be repeated or advertised as a missing first backreaction experiment.

The [universal-clock analysis](../universal-clock-coupling/report.md) separately states that universal coupling includes matter energy and can cancel a homogeneous redshift. Here those two requirements are combined with the autonomous receiver, its analytic solution, an explicit matter/photon energy split and finite two-event ray calculations. A partial-clock-response diagnostic tests whether weakening universality avoids the cancellation while keeping local light speed constant. The earlier cancellation argument was already known; this calculation does not claim its discovery.

## Conditional Hamiltonian and its source

**Proposed effective model using known Hamiltonian mathematics; originality unverified.** Use a homogeneous positive field n(t), fixed spatial rulers, c=1 and a finite periodic reference volume. Uniform matter is at rest and the background radiation is homogeneous and direction-balanced. Spatial gradients vanish, so matter at rest remains at rest. This is an exact homogeneous truncation of a stipulated model, not a realistic resolved galaxy/void distribution.

\[
q=1/n,\qquad d\tau=dt/n,\qquad
H=\frac{\Pi^2}{2K}+\frac{E_{\gamma0}+E_{m0}}{n},\qquad K>0.
\]

E_gamma0 is the sum of the radiation's constant local photon energies in this homogeneous model; E_m0 includes the matter rest energy. Atomic test clocks have negligible backreaction but share the same multiplier q. The clock/lapse kinematics use established geometric concepts described in [Gourgoulhon's 3+1 lecture notes](https://arxiv.org/abs/gr-qc/0703035). The choice of a physical field kinetic term and universal coupling here is additional hypothetical physics; it is not derived from Einstein's equations.

**Conditional deductions from that Hamiltonian, not new fundamental laws:**

\[
\dot n=\Pi/K,\qquad \dot\Pi=(E_{\gamma0}+E_{m0})/n^2.
\]

For n(0)=1 and Pi(0)=0, let U=E_gamma0+E_m0. Then

\[
\frac{K\dot n^2}{2}=U(1-1/n),\qquad
t(n)=\sqrt{\frac{K}{2U}}\left[\sqrt{n(n-1)}+\operatorname{asinh}\sqrt{n-1}\right].
\]

The radiation and matter reference-energy losses are respectively E_gamma0(1-1/n) and E_m0(1-1/n). Their sum equals the field kinetic-energy gain. The fraction supplied by radiation is therefore E_gamma0/U. Raising matter loading speeds up the evolution as U^(-1/2) at fixed target n and K; it does not leave the previous photon-only evolution intact.

At n=2 and E_gamma0=1:

| Assumed matter/radiation energy ratio | Photon reference-energy loss | Matter reference-energy loss | Field gain | Radiation share |
|---|---:|---:|---:|---:|
| 0 | 0.5 | 0 | 0.5 | 100% |
| 1 | 0.5 | 0.5 | 1 | 50% |
| 10 | 0.5 | 5 | 5.5 | 9.09% |
| 1,000 | 0.5 | 500 | 500.5 | 0.0999% |

These are deliberately chosen dimensionless inputs, **not measured cosmic energy ratios**. Local masses do not decrease in this prescription: their reference Hamiltonian changes with q. The result is about coupling and consistent accounting, not the deferred total starlight supply calculation. A theory cannot interpret the entire receiver gain as converted starlight while retaining this universal Hamiltonian.

## What the actual clocks measure

**Established ray and clock bookkeeping, conditional on the model:** photon momentum P is constant because the background is spatially homogeneous. Its reference energy H_gamma=P/n falls; its local energy H_gamma/q=P does not. For emission e and observation o,

\[
S_{\rm reference}=n_o/n_e,\qquad
1+z_{\rm measured}=S_{\rm reference}\,q_o/q_e=1.
\]

The same conclusion follows from an independent ray calculation. With dx/dt=1/n, a fixed source–detector distance D obeys

\[
D=\int_{t_e}^{t_o}\frac{dt}{n(t)}=\tau_o-\tau_e.
\]

For any two emissions, not just infinitesimally separated ones, Delta tau_o=Delta tau_e. Reference arrival spacing changes, but the receiving clock changes with it. There is no measured spectral or temporal stretch to compare favorably with redshifted astronomical sources. This statement needs no expansion premise or supernova template assumption.

The numerical test uses five matter/radiation ratios (0, 0.1, 1, 10, 1,000), two inertias (0.5, 2), four field checkpoints and nine ray/event comparisons per background. It solves 10 backgrounds, records 40 energy samples and 90 event pairs, then repeats with tighter tolerances. Distances and launch times scale with sqrt(K/U), deliberately comparing equivalent phases; these are not 90 independent astronomical observations.

- Maximum relative total-energy error: 7.54e-13.
- Maximum analytic/numerical time disagreement: 9.06e-13 in scaled time units.
- Reference redshifts: 0.0201–2.8962; maximum measured local redshift: 2.23e-16.
- Maximum finite-event proper-time stretch error: 3.50e-13.
- Numerical reference-arrival derivative agrees with n_o/n_e to 5.36e-6 relative; this is a one-sided finite-difference check.
- Tightening solver tolerances changes reference redshift by at most 2.62e-11.

## Can a weaker clock response repair this version?

**Proposed diagnostic parameterization, not a derived atomic model:** keep photon speed dx/dt=1/n and fixed spatial rulers, but assign clocks d tau=dt/n^beta. The universal case is beta=1. Standard ray bookkeeping gives

\[
1+z=(n_o/n_e)^{1-\beta},\qquad
c_{\rm local}=n^{\beta-1},\qquad
\frac{c_{{\rm local},o}}{c_{{\rm local},e}}=\frac{1}{1+z}.
\]

The finite event map follows the same proper-time rule; its infinitesimal stretch equals the displayed frequency stretch, although a finite event can have a varying stretch factor. A fixed physical ruler response is an explicit assumption here. Changing the ruler law would require a new matter/spatial completion and an audit against the user's nonexpansion and local-normality premises.

For the **earlier empirical model's prediction**, z=0.00766048 at approximately 100 million light-years, beta values 0, 0.25, 0.5, 0.75, 0.9 and 0.99 all require the same local-speed ratio 0.99239776, a 0.7602% change. At beta=1 no finite n ratio produces this nonzero shift. This number is a conditional consequence of our fitted distance relation, not a new measurement of remote light speed. Reducing beta cannot independently restore this redshift and retain the stipulated constant local c in this homogeneous, fixed-ruler family.

## Research decision and remaining scope

Do not promote the homogeneous universal receiver to an observed-redshift mechanism or use its total energy gain as photon-created deposits. The rising homogeneous field energy also does not implement the requested lossless traveling companions with no void accumulation. There is no capture or lensing solution here.

Keep the distinct routes explicit: a spatially inhomogeneous evolving field with derived matter and endpoint behavior, or a nonuniversal conversion interaction with derived phase/event transfer and a consistent atomic/ruler sector. Earlier spatial and microscopic candidates already exist; their documented failures and assumptions must guide revisions rather than resetting the program. Adding prescribed endpoint protection or a fitted per-source response would not close these requirements.

The immediate acceptance gate for a next candidate is one autonomous source–propagation–detector calculation that reports local clock standards, signal color, finite-event duration, local speed, all coupled energy sources, and companion transport together. The gravity response must then break the [known amplitude degeneracy](../theory-priority-review/report.md) before a favorable stellar fit can support photon origin. The bar-aware orbital likelihood and existing observational coverage remain required, not replaced by this gate.

No stellar/galaxy holdout scores were opened, no data were retuned, and the full research goal remains active. Reproduce with `python research_work/results/matter-clock-closure/run.py`; `results.json` retains every synthetic sample and the partial-response cases.

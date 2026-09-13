# Near and far sources as companion streams

This explores the user's river analogy through standard transport mathematics and the proposed photon-to-companion conversion. The exact-one-third retention law remains unchanged. New calculations concern source flow and geometry; no successful halo formation, deposited energy inventory, or cluster fit is claimed.

## What a stream means

Each source contributes an incoming energy intensity, labeled by its direction and emission time. Under straight paths, no intervening capture, isotropic emission and a stationary-source approximation:

    F_c,i(x) = L_i/[4*pi*d_i(x)^2] * [1-exp(-alpha*d_i(x))]
    q_dep(x) = eta_1/3(X) * kappa(x) * sum_i F_c,i(x)

The inverse-square and absorption formulas are known. The conversion alpha and companion capture interpretation are proposed. In the general time-dependent calculation replace the luminosity by the emission history appropriate to the retarded path, include intervening transport, and track sources inside the receiver separately. Do not apply the point-source singularity inside a resolved stellar system.

Capture uses the sum of incoming intensities, not the size of their net vector. Two equal streams arriving from opposite directions have zero net flux vector but both deliver energy to an absorber. The script checks this explicitly. Independent incoherent streams can cross without colliding or following the net-flux arrows; river junctions require an actual scattering/binding interaction, not just a visual analogy.

## Where nearby becomes far

There is no universal boundary in light-years. The controlling ratio is source distance D divided by receiver radius R. For our small alpha*D regime, companion intensity from a single source scales approximately as alpha*L/(4*pi*D), because conversion increases with distance while geometric dilution decreases it. This differs from the inverse-square photon flux alone.

Near-side versus far-side companion illumination is approximately (D+R)/(D-R):

| D/R | Near-side / far-side intensity | Consequence |
|---|---:|---|
| 2 | 3.00 | Strong gradient across the receiver |
| 5 | 1.50 | Source location can shape the input |
| 10 | 1.22 | Moderate gradient |
| 20 | 1.11 | Becoming nearly uniform |
| 100 | 1.02 | Almost uniform across the receiver |

The code computes the full exponential expression, not just this small-conversion limit. At conversion saturation, the intensity approaches inverse-square behavior and the ratio becomes the square of the expression above. Sources close to the receiver boundary require exact entry/chord geometry; this center-flux diagnostic does not replace it.

For a 100 kpc receiver, a source 1 Mpc away has D/R=10; a source 10 Mpc away has D/R=100. For a 1 Mpc cluster receiver, those same ratios occur at 10 and 100 Mpc. Distant sources can deliver directional beams that vary little across one galaxy. Directional is not the same as spatially concentrated into a galaxy-centered halo.

## Actual catalog calculation around Virgo/M87

Reuse the cached DustPedia snapshot, with 814 usable galaxy entries represented as sources, including 797 outside the illustrative 1 Mpc M87 receiver. These are galaxies' integrated stellar luminosities, not individually resolved stars. The luminosities are SED-model estimates with prior distance assumptions, not directly measured bolometric source histories. The receiver radius is the existing pilot convention, not an inferred halo edge.

External companion intensity at the center, THEMIS luminosity version:

| Source distance from M87 | Share within this catalog |
|---|---:|
| 1–2 Mpc | 30.8% |
| 2–5 Mpc | 22.4% |
| 5–10 Mpc | 14.3% |
| 10–20 Mpc | 17.2% |
| 20–50 Mpc | 15.0% |
| 50–100 Mpc | 0.36% |

The second luminosity model gives very similar fractions. The top ten sources supply about 24.5% of the catalog external intensity. Its angular distribution is uneven: the directional second-moment eigenvalues are approximately 0.127, 0.359 and 0.514, compared with 1/3 each for isotropic illumination. The net directional fraction is 0.195; this is not a capture efficiency or fraction of energy lost.

These sources extend only to about 72 Mpc separation. Therefore the result says nothing about how much the omitted billion-light-year population contributes. Nearby dominance in this incomplete catalog must not be extrapolated to the entire universe. The earlier exact-chord pilot's absorbed-power fractions by source-distance bin are included separately in the JSON. They are not identical to center-intensity shares and retain the old illustrative kappa*R=10 assumption.

No cluster X is invented here. At fixed receiver X, the one-third factor multiplies all incoming contributions equally and cancels from their fractional shares. Its absolute cluster application still requires the shared cluster proxy definition. A local varying X would change the deposited pattern and must be specified before fitting.

## Can distant shells outweigh nearby sources?

Yes, conditionally. For uniform luminosity per volume j, a shell has luminosity 4*pi*D^2*j*dD. Its geometric dilution cancels that area factor, giving

    dF_c = j*[1-exp(-alpha*D)] dD
    F_c(<Dmax)/j = Dmax - [1-exp(-alpha*Dmax)]/alpha.

These are known geometric/integral identities applied to our conversion postulate. For alpha*D small, the integrated result is approximately alpha*Dmax^2/2. Thus a growing collection of distant stars can indeed outweigh individual nearby sources. Once conversion saturates, each equally thick uniform shell contributes approximately equally.

This is a parameterized homogeneous illustration, not a measured luminosity field. An eternal, static, transparent uniform population extended infinitely far gives divergent intensity, so a finite physical source history, intervening capture, spatial decline, or another explicit transport limitation is needed. We do not impose a universe age or size; Dmax remains a calculation parameter and the actual source history/opacity must determine which energy arrives.

## Can streams be steered into the required halo?

Three physically distinct branches should be kept separate:

1. Straight streams with capture: no extra steering. Nearby sources establish gradients; the receiver's capture coefficient determines where energy stays. This is the reference calculation above.
2. Lightlike gravitational deflection: in weak gravity, the known ray equation has dn/ds approximately -2*grad_perpendicular(Phi)/c^2 under the usual equal-potential assumptions. It can focus passing beams but does not bind them. Our six lens-angle constraints correspond to total deflections of order 10^-5 radians in the retained geometry, not large river-like turns. This fact does not rule out focusing; it limits a literal picture of beams turning sharply into the galaxy. A full ray calculation has not been executed here.
3. Proposed companion-specific scattering or binding: streams can change direction or convert to retained states through a new interaction. Use a transport equation with an angular redistribution kernel and separate capture/release channels. Direction changes must account for momentum transferred to matter or fields. This branch is not established by naming a focusing coefficient or fitting a halo.

For all branches the target remains: predict q_dep(x,t), integrate its history with transport, obtain u_dep(x), then compare both motion and lensing from the same field. A sky pattern is an input to that calculation, not the halo density itself. In particular, sources illuminating a receiver uniformly do not generate an NFW radial profile unless the capture, steering or subsequent dynamics produces it.

## Practical next test and status

Use the resolved nearby-source field plus a separately normalized distant background; keep the one-third retention factor. Ray-trace their passage through one declared capture profile and compare the predicted angular and radial deposition with the finite halo targets. Fit any additional shared capture/steering parameter on designated training targets, and test on others. Do not tune distant source strengths to the halo that they are meant to predict.

Completed here: exact near/far contrast curves, catalog source-direction moments and distance contributions, comparison with prior exact-chord absorbed fractions, finite homogeneous-shell integrals, and a lightlike bending-scale calculation. Not completed: actual billion-year source histories, companion-specific steering, self-consistent binding, or a new rotation/lensing fit. Catalog hashes and normalization checks are recorded in companion-streams-results.json; the reproducible script is companion-streams.py. All six research goals remain open.

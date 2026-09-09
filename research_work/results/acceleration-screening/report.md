# A shared weak-acceleration switch changes gravity in the wrong direction

## Provenance and model

Modified-Poisson equations derived from actions are established in [Bekenstein and Milgrom 1984](https://adsabs.harvard.edu/pdf/1984ApJ...286....7B). Environmental scalar screening is also established; [Hinterbichler and Khoury 2010](https://arxiv.org/abs/1001.4525) is an example based on density. Neither is evidence that our proposed coupling is correct. The trial combination below is proposed here, with originality unverified. We do not adopt a cosmological expansion model or the cited theories' full phenomenology.

Proposed trial action density, with standard Newtonian/scalar ingredients and a stipulated interaction:

L = -|grad Phi|^2/(8 pi G) - rho Phi
    + K phi_dot^2/2 - K v^2 |grad phi|^2/2
    - K m0^2 phi^2 f(Y)/2,

Y=|grad Phi|^2/a0^2, f(Y)=Y/(1+Y).

K>0, m0^2>0 and a0>0. The environmental restoring term is nonnegative in the field energy. Its low-acceleration limit permits the flat-potential behavior examined earlier, while high acceleration restores the field toward zero. The switch is a proposed shape, not independently derived or calibrated. Acceleration is defined in a declared preferred frame; a covariant completion would need a physical observer field. A force-cancellation point also has Y=0, so this is not a reliable void identifier by itself.

## Gravity reaction cannot be omitted

Derived here by the established Euler-Lagrange method; originality unverified:

K(phi_tt-v^2 Laplacian phi)+K m0^2 f(Y) phi = Q_phi,

where Q_phi would have to come from the joined photon interaction, not a free fitted source. Varying Phi gives

div[mu(Y,phi) grad Phi]=4 pi G rho,

mu=1+B/(1+Y)^2, B=4 pi G K m0^2 phi^2/a0^2.

The two equations are consequences of one interaction. Keeping the scalar switch while omitting the extra term in mu would not be variation of this same action. This analysis does not yet include a complete photon action, a relativistic stress-energy source, or physical clocks/rulers, so it is not a finished energy-conserving cosmic model.

Established spherical flux integration, applied to the derived coefficient:

mu a = G M(<r)/r^2 = a_Newton.

Since mu>=1, this direct modification gives a<=a_Newton for fixed enclosed rho. It cannot itself supply the requested excess attraction. Additional gravitating deposited energy would be a separate computed source; it is not included or ruled out by this sign result. General nonspherical force directions are not determined by this simple algebraic relation.

## Static operator check

Established principal-symbol/ellipticity test applied to the candidate: transverse eigenvalue is mu and longitudinal eigenvalue is mu+2Y mu_Y. The latter reduces to

lambda_parallel=1+B(1-3Y)/(1+Y)^3.

At fixed B>0 its minimum for Y>=0 is 1-B/4 at Y=1, verified symbolically. Thus B<4 keeps this scalar gravitational operator strictly elliptic throughout the acceleration domain; B=4 is degenerate and B>4 gives a negative eigenvalue in part of it. This is loss of the usual elliptic boundary-value character, not by itself proof of a time-dependent ghost instability. A coupled dynamical stability analysis remains separate.

Increasing coupling indefinitely therefore does not simply increase screening safely. Reversing the interaction sign changes the potential and the field dynamics; it cannot be used as an unexplained fitting fix.

## Outcome

This first shared-action candidate fails to provide both the desired low-gravity activity and direct extra attraction. It is not a no-go theorem for modified gravity. A better physical selector (for example curvature measured in a specified local frame), another interaction or a derived deposited source may behave differently, but each requires its own conservation, source and stability analysis. No candidate is promoted and no astronomical parameters are fitted.

Run `python research_work/results/acceleration-screening/run.py` for the symbolic coefficient/eigenvalue checks. The protocol fixes the tested interaction before execution. Constant local light speed, endpoint clocks, line fidelity, lensing and an independently predictive redshift rate remain required and unresolved.

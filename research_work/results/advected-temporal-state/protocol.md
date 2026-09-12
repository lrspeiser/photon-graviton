# Directly advected temporal energy: local propagation audit

Replace separate stationary memory by a traveling state-energy density C with n=1+C/K, K>0. Keep radiation R=G+H and shared EM/GW ray speed v=1/n in units c0=1. Postulate receiving-state energy flux eta*v*C, first eta=1 and a sensitivity eta=0.5. This changes the prior memory closure explicitly; it is not a proven graviton law. The clock rule remains d tau=dt/n.

Equations: R_t+(vR)_x=-(n_t/n)R; C_t+(eta*v*C)_x=(n_t/n)R-Gamma*C; D_t=Gamma*C. Their sum conserves reference energy, and all temporal-state energy is carried by C. Substitution gives d*C_t+(eta/n^2)*C_x=-Gamma*C, with d=1-R/(K*n). Test whether this locally defined system has sensible propagation before a spatial simulation, rather than forcing numerical evolution through d=0.

For source-free perturbations in (R,C), derive M*y_t+B*y_x=0 with M=[[1,R/(K*n)],[0,d]], B=[[v,-R/(K*n^2)],[0,eta/n^2]]. Compare its numerical eigenvalues against v and eta/(n^2*d). Local clock speeds are these values times n; distinguish energy flux divided by density from the speed of a perturbation of the coupled system.

Before calculation use K=0.05,0.5,2; R=0.01,0.22,1; C=0.001,0.01,0.22,1; eta=1,0.5. Archive all 72 states. Add an exactly singular state K=0.5,C=0.01,R=0.51 and a coincident-speed state K=0.5,R=C=0.22,eta=1. For the latter verify the exact Fourier perturbation amplification for wavenumbers 1,10,100,1000 using an independent matrix exponential.

Also inspect the homogeneous capture sign: C_t=-Gamma*C/d. At d<0 the receiving sector grows despite the capture sink because radiation-to-state feedback overcompensates; at d=0 the algebraic evolution is singular unless the numerator satisfies a compatibility condition. Do not call negative characteristic speed inherently acausal or declare the model globally impossible; identify the domain and violated intended speed/capture assumptions. No observation fit, global stability proof or full momentum completion is claimed.

Follow-up declared after the local grid: test whether capture preserves the eta=1 forward/subluminal domain C>R. Start C=1, R=0.22, Gamma=1 at K=0.5 and 2. Derive the invariant R(K+C)=A, solve the C=R crossing algebraically, and independently compare direct ODE event time with quadrature to 1e-8. Stop direct integration at this crossing; do not integrate through a singularity. Track R+C+D to 1e-9. These two synthetic cases test domain preservation, not observation agreement.

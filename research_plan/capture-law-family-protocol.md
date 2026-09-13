# Small shared-law comparison: source-to-retention scaling

Before calculation, compare constant retention and five two-parameter power-law families: size, well depth, surface acceleration, mean-density proxy, and emitted flux per surface area. These are proposed response families using known algebra, not derived interactions. Do not call equivalent reparameterizations independent mechanisms.

First-stage target: the user-requested full-gravity Sun/Earth/Moon interpretation, extended to the other planets. Use the prior outgoing-power proxies and fixed D=R conversion geometry. All masses and source data have been exposed; this is exploratory transfer. Fit the normalization on Sun and one exponent on Earth, then freeze for Moon and the other eight tabulated bodies. Their existing gravitational masses are targets, not extra ordinary mass added on top.

Write mu=M/M_sun, r=R/R_sun, l=L/L_sun and p_c=P_c/P_c,sun. Candidate scaling is

    mu = p_c * [mu^a*r^b*l^d]^p.

Families (a,b,d): constant (0,0,0) with p=0; size (0,1,0); depth (1,-1,0); acceleration (1,-2,0); density (1,-3,0); flux (0,-2,1). Reference normalization is K_sun=M_sun*c^2/P_c,sun. K is eta*T for constant lossless storage or eta*tau for steady renewal, not separately identified efficiency or age.

For mass-dependent predictors, never insert the target test mass to make a prediction. Solve the implicit equation:

    mu_pred = exp{[ln p_c + p*(b ln r+d ln l)]/(1-a*p)}.

If 1-a*p vanishes, the relation is degenerate and cannot determine mass. This guards against using the measured gravity as a predictor of itself. Fit exponents algebraically on Sun/Earth only; no selection/refitting using Moon or other bodies. Report full prediction ratios and log-RMS. Also report whether stronger gravity/density increases or decreases the fitted retention. Negative exponents are allowed diagnostics, not silently declared deeper-well capture successes.

These scale laws describe effective retention only; capture probability must remain <=1 and a physical storage/history law must supply any effective K. Distinguish stable body-bound storage, independently self-bound states (no support law yet) and steady renewal. A static fit cannot distinguish the first and third. No frequency or graviton count is determined by fitting masses.

Galaxy analyses use ordinary matter plus additional gravity, unlike the full-mass target here. They cannot be silently combined into the same likelihood. Inspect the resulting first-stage candidates before claiming a universal law or applying one to galaxy residuals. The existing galaxy fit failures remain evidence, not a reason to use per-object rescue parameters.

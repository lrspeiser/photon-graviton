# CMF-2: ordinary-source profile geometry follow-up

Declared after CMF-1 and before these calculations. CMF-1's joint candidate has validation/test galaxy errors 23.99/21.72 km/s, versus the newly matched simple-MOND 26.14/16.09. Its galaxy-only candidate reaches 15.65 test-labelled km/s but cluster-pressure chi2/point is 195.64. The paired galaxy-only improvement interval includes zero. Primary Coma joint unbounded beta is 2.94/1.72, beyond static source geometry. These failures are retained, not overwritten.

Hypothesis: the slope of the ordinary-source acceleration profile helps distinguish source geometry and the response retained from smaller radii. This is an exposed, data-motivated extension. It is not an independent confirmation of CMF-1, a demonstrated vector swirl or an original mathematical method.

Add s=tanh[(d ln(max(g_b,1e-30 SI))/d ln r)/2], evaluated by numpy.gradient on the supplied increasing-radius source grid (first-order endpoints; zero derivative for a one-row profile). No observed velocity or pressure is differentiated. s is a force-profile shape proxy; interpreting it as actual mass-density slope is valid only for a spherical ordinary source and positive g_b.

G = RM basis plus s,x*s,y*s,z*s,s^2 (16 coefficients).
GH = RMH basis plus s,x*s,y*s,z*s,s^2,h*s (20 coefficients); ell=0.3,1,3.
Same positive bounded force law, coefficient bounds and high-acceleration release as CMF-1. No galaxy/cluster label switch. Four structures, three ridge values (0,0.001,0.01), four cluster weights (0,0.1,1,10): 48 fits, three seeded starts each, 144 attempts. All signs are permitted. RNG seed 200927.

Before fitting, independently check each new-family objective Jacobian. Preserve all original CMF-1 source files. For operational candidate selection, pool old and new optimizer-success fits, using CMF-1's fixed validation objective for the joint choice and validation galaxy RMSE for galaxy-only choice. Save selection before recomputing test/Coma. All prior test information is exposed and the follow-up is labelled accordingly. The best new-only candidate and pooled choice are both reported. Retain old MOND and IH-1 references; do not refit them differently. Repeat paired object bootstrap, primary Coma geometry and 3/30 Mpc reach sensitivities.

Repeat 600->1200 cluster grid and 4096->8192 source/memory plus optical refinement for the selected pooled joint law. Also test the selected galaxy-only memory law's inner-boundary sensitivity; the CMF-1 joint winner did not use memory, so its small boundary sensitivity was not evidence for boundary robustness of the galaxy-only law. A numerical refinement failure prevents a numerical pass, even if the score looks attractive. For galaxy-source slope discretization, add midpoint radii, interpolate log ordinary acceleration and reevaluate memory/derivatives at original radii; report resulting predicted-speed changes as input-grid sensitivity, not observational evidence.

Success criteria, fixed distances, no dark matter, static space, evidence preservation, existing attribution and missing field-energy/source/light-law derivations remain exactly those of protocol.md. No parameter-count advantage is claimed. No additional candidate families will be introduced into this campaign based on its test scores.

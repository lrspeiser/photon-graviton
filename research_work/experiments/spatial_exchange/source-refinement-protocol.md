# SE-R: convergence of generated scalar source

Declared 20 September 2026. SE-1's wave-source derivative changed 10.92% on
refinement despite passing channel-energy comparisons. Measure source convergence
directly for both positive-response SE-2 chi=200 emitter choices, theta=pi/4
and pi/2. No selection based on the better final source is allowed.

Keep all SE-2 equations and physical parameters fixed. Use the already pinned
emitter_model.py. Seven new T4 runs, L16,radius .9,chi200,mix .5,dt .01:
mixed emitter at n32,40,48,56; Y emitter at n48,56; mixed emitter n40 rotated
pi/3. Reuse the SE-2 mixed/Y n32 dt.02, Y n32 dt.01, Y n40 dt.01 and Y rotation
states when they finish; do not restart or overwrite those runs.

All new runs retain SE-2 energy/cone/edge/clearance gates. Record initial/final
states and .1-time metrics. After completion, compute wave and total coupling
source using the independent SE-T reconstruction and two derivative checks.
Archive all values, not just comparisons that pass.

For each angle require the n48-to-n56 relative wave-source difference <5%,
and smaller than the n40-to-n48 difference; total coupling source difference
<0.1%. Denominator is max(abs(finer source),1e-12). Compare n32 dt.02 to .01:
wave source difference <0.1%, total source <0.01%. At matching resolution,
rotation wave source difference <2%, total source <0.1%. These are new declared
accuracy gates; do not retroactively change SE-1/SE-2's original results.

This is still a pre-boundary source calculation. It cannot establish a correct
radial force, lasting swirl, orbit, lensing prediction or source lifetime.
All full-model observational requirements remain. Hamiltonian methods and grid
refinement are established mathematics; candidate constitutive laws are not
claimed to be established physics or previously unconsidered inventions.

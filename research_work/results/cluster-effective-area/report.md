# Effective collecting area from a smooth capture well

The external-capture model can use an effective area instead of a hard outer shell. For a spherical, fixed capture environment with straight companion trajectories, let kappa(r) be capture probability per unit path length. Established absorption geometry gives

    tau(b)=integral_-infinity^infinity kappa(sqrt(b^2+z^2)) dz
    sigma_eff=2 pi integral_0^infinity b [1-exp(-tau(b))] db
    deposited power=c u_companion sigma_eff
    R_eff=sqrt(sigma_eff/pi).

b is the trajectory's closest approach to the center. sigma_eff already includes capture efficiency, so it must not be multiplied by that efficiency a second time. For a fixed surrounding isotropic companion energy density, the result is exactly the user's collecting-area rule. The radius is now derived from a specified capture profile, not assigned as the physical edge of a gravity well. No internal stellar supply is required by this equation.

The geometry and exponential absorption are known mathematics; treating companions as capturable and specifying kappa from a gravity field are hypotheses. We have not yet derived kappa from well depth or a particle interaction. Gravitational focusing, moving wells and changes in the external bath are omitted here; they require trajectory/transport calculations rather than automatic reuse of the straight-ray area.

## Which profiles yield finite collecting area?

Suppose kappa falls as r^(-p) at large radius. For p>1 the chord optical depth scales as b^(1-p). Far from the center, absorption is weak and the area integrand scales as b^(2-p). Thus an isolated, unbounded, steady straight-ray model has finite area only for p>3. At p=3 the area grows logarithmically with the outer limit; for 1<p<3 it grows as a power. For p<=1 even the infinite chord optical depth diverges. These conclusions assume a nonzero asymptotic power-law coefficient.

This matters for possible gravity couplings. Outside a finite mass, a rate proportional to acceleration magnitude (~r^-2) would fail this isolated-area condition; a rate proportional to acceleration squared (~r^-4) could satisfy it. A rate proportional to tidal-field magnitude (~r^-3) is marginal; its square (~r^-6) converges. A potential-based power law has to fall faster than |Phi|^3 when |Phi|~1/r. These are asymptotic candidate constraints, not endorsements of a particular coupling, coordinate-invariant field action or local screening prescription.

Finite source ages, competing receivers and overlapping fields can alter the isolated-infinite-bath problem. Divergence here does not exclude every model with a shallow local rate, but it prevents assigning a finite isolated cross-section to that rate without additional physical transport structure.

## Smooth worked family

Use the optional positive profile

    kappa(r)=kappa0 [1+(r/a)^2]^(-p/2).

This smooth toy profile has no hard cutoff. a and kappa0 remain unspecified physical scales, not fitted observed cluster quantities. For p>1 the exact optical depth is

    tau(b)=T [1+(b/a)^2]^((1-p)/2)
    T=kappa0 a sqrt(pi) Gamma((p-1)/2)/Gamma(p/2).

For q=(p-1)/2>1, established integration gives

    sigma_eff/(pi a^2)=T^(1/q) gamma(1-1/q,T)-[1-exp(-T)],

where gamma is the lower incomplete gamma function. The weak-capture limit is sigma_eff/(pi a^2)~T/(q-1); the strong-capture leading term grows as T^(1/q). Stronger capture therefore increases the area, with diminishing response as central trajectories become fully absorbed. Doubling a gives four times the area only at fixed dimensionless strength kappa0 a and profile shape; at fixed kappa0 the dimensionless strength changes too.

| Tail p | kappa0 a | R_eff/a |
|---|---:|---:|
| 4 | 0.1 | 0.555 |
| 4 | 1 | 1.634 |
| 4 | 10 | 3.975 |
| 6 | 0.1 | 0.277 |
| 6 | 1 | 0.804 |
| 6 | 10 | 1.730 |

The p=2 and p=3 cases instead keep growing as impact limits increase from 10a to 100a to 1000a; they have no finite R_eff in this idealization. The test evaluates twelve profiles and independently integrates 36 chords. Finite-area analytic and numerical calculations agree within 1.12e-15 relative; chord differences are at most 4.40e-14. These numerical checks are not observational validation or a bound on modeling errors.

## Connection to the research goal

This makes the collecting-area proposal concrete while retaining the need for a common gravity-to-capture law across galaxies and clusters. In the preceding competing-receiver model, beta=sum(n_i sigma_i); increasing all collecting areas also changes the surviving companion bath. The fixed-bath area scaling is not permission to keep that bath unchanged when capture across the whole universe changes.

Next specify a candidate kappa tied to ordinary matter/field structure, calculate its spatial deposits and source competition, and test the same response against motion and lensing. No cluster radius, supply normalization or capture coefficient has been fitted to a desired mass. All six objectives remain open and no holdouts were opened.

Run `python research_work/results/cluster-effective-area/run.py`.

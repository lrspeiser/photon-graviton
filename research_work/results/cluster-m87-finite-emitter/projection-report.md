# Resolved projection of the M87 contribution

The observationally informed finite emitter now produces a resolved spherical deposit-power profile, rather than only a proof that its central value is bounded. Its central projected power is 2.401985e8 Lsun/Mpc^2. This is the M87 contribution alone, not the full cluster or an observed lensing surface density.

| Projected radius (kpc) | Surface power (Lsun/Mpc^2) | Power within sphere (Lsun) | Power within projected circular aperture (Lsun) |
|---|---:|---:|---:|
| 1 | 2.35813e8 | 32.66 | 746.85 |
| 10 | 1.87358e8 | 12,410 | 63,815 |
| 100 | 6.81058e7 | 1,364,494 | 2,962,034 |
| 300 | 2.45383e7 | 7,664,818 | 12,296,651 |
| 500 | 1.24384e7 | 14,992,717 | 20,878,254 |
| 1000 | 0 | 33,684,549 | 33,684,549 |

At the 1 Mpc boundary the line-of-sight chord length is zero, so surface power vanishes there; enclosed power remains the full total. A projected aperture includes deposits in front of and behind the corresponding sphere, explaining its larger enclosed power. These quantities are useful for a future common motion/lensing calculation, but neither is accumulated mass.

## Method and evidence

Inputs, photometry provenance and explicit spherical/bolometric assumptions are documented in [the finite-emitter report](report.md); the input result hash is stored in projection.json. Nothing is fitted to gravity or lensing. No additional galaxy, total luminosity, capture rate or time history is introduced.

The established convolution q(x)=integral j(y)K(|x-y|)d^3y is evaluated using the exact angular average of each spherical Gaussian. With r the deposition radius, s the emitter-to-deposit separation and sigma the Gaussian width, that average is proportional to

    exp[-(r-s)^2/(2 sigma^2)]
    * [1-exp(-2rs/sigma^2)] / [2rs/sigma^2].

The code handles the zero-argument limit and uses expm1 to avoid cancellation. Changing integration variable to x=(s-r)/sigma resolves narrow emitters even at large r. Integration covers x in [max(-12,-r/sigma),12]; omitted Gaussian tails are negligible at the tested precision. This retains the preceding approximation that light outside the receiver is negligible.

Positive radial densities are interpolated in log radius/log density on 300 and 600 nodes from 1e-8 to 1 Mpc. The innermost region uses the first sample; the analytic center is finite and the cutoff is far smaller than the smallest source component. Projection uses the known line integral

    S(b)=2 integral_0^sqrt(1-b^2) q(sqrt(b^2+z^2)) dz.

Surface-profile refinement differences are at most 2.49e-7 relative over the tabulated nonzero samples. The independently integrated three-dimensional power is 33,684,550.07 Lsun; the projected-image integral is 33,684,549.44 Lsun; the earlier emitter/ray calculation gives 33,684,549.22 Lsun. Relative discrepancies are 2.55e-8 and 6.61e-9. These are numerical consistency checks, not data uncertainties or validation of the physical mechanism.

Projected aperture power is also computed directly from spherical shells: each shell with radius r>b contributes fraction 1-sqrt(1-b^2/r^2), while shells inside b contribute fully. This produces both spherical and projected enclosed power from the same q. Nonnegative density, monotone tabulated surface power, aperture ordering, global conservation and the previous analytic upper bound are checked in project.py.

All these equations are established geometry/integration applied to hypothetical conversion and capture. They are not new gravity equations.

## Interpretation and remaining work

The artificial infinite M87 peak is resolved without changing the transport rates or adding energy. Other emitters and their directional contributions have not yet been incorporated into a resolved full-cluster map. A source history is needed to turn power into stored energy; the reservoir's stresses, support, movement and gravitational response are needed to turn that energy into predictions of stellar motion and bending. An arbitrary exposure time or lensing normalization has not been fitted. All six objectives remain open and final holdouts remain unopened.

Run `python research_work/results/cluster-m87-finite-emitter/project.py` after run.py.

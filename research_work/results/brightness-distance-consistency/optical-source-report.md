# Conditional energy-source bridge between optics and capture

The retained brightness geometry requires positive null focusing at all 960 evaluated supernova redshifts. Under the additional assumption that an isotropic massless companion bath supplies all of this focusing, its observer-limit energy density is **4.71405e-9 J/m^3**. At a constant bath of that density, the independently fitted galaxy capture exposure corresponds to **933.20 million years**. This is a conditional relation between two fitted model requirements, not evidence that the bath exists or that stellar energy supplies it. All six goals remain open.

## What is known and what is postulated

The Sachs optical equation relates beam focusing to spacetime curvature. With Einstein coupling, the curvature contraction is related to stress-energy. These are **known relations**, not new equations of this project; see [Perlick, Gravitational Lensing from a Spacetime Perspective, section 2.4.5](https://link.springer.com/article/10.12942/lrr-2004-9). We apply their focusing convention as D_A''=-R_opt D_A for a shear-free background, where positive R_opt focuses the beam. No expanding-universe solution is assumed.

Our **postulated, fitted** beam law remains

    f=z/(1+z), lambda=f/alpha,
    H(f)=alpha*D_A=-(1-f)*ln(1-f)*sqrt(1+f/(1+q*f)),
    alpha=0.0002488993286382367/Mpc, q=2.2601388650845373.

This analysis additionally assumes metric-geodesic photons, observer-normalized affine length lambda, Einstein coupling, zero background optical shear, and an isotropic perfect-fluid source. In particular, this is not a proof that a direct photon-conversion law has this affine/redshift relation.

Known stress-energy contraction gives, with rho a mass density and p a pressure,

    R_opt = -alpha^2 H''/H
          = (4*pi*G/c^2)*(rho+p/c^2)*(1+z)^2.

Thus the fitted optical geometry determines a required combination rho+p/c^2, not the density and pressure separately. The **additional companion hypothesis** p=u/3 and rho=u/c^2 yields

    u_required(z) = 3*c^4*R_opt(z)/[16*pi*G*(1+z)^2].

It assigns all background focusing to isotropic massless companions. Ordinary matter, other stresses, or optical shear would change that allocation. Collinear companions accompanying only one ray are not an isotropic bath and do not justify this formula by themselves.

## Results

| Redshift | Required all-companion energy density (J/m^3) |
|---|---:|
| 0, observer limit | 4.71405e-9 |
| 0.1 | 2.84165e-9 |
| 0.3 | 1.57481e-9 |
| 0.6 | 9.76690e-10 |
| 1 | 6.82651e-10 |
| 2 | 4.35682e-10 |
| 3, illustrative extrapolation | 3.45280e-10 |

These are inferred source requirements of the fitted formula, **not measured companion densities**. No extra fitting was done. The calculation uses zHEL, matching the fitted optical response; observer and representative values are tabulated separately from the actual sample. On the 960 sample rows, R_opt ranges from 3.28961e-7 to 5.67593e-7 /Mpc^2 and no value is negative. The analytic observer limit is 6.21393e-7 /Mpc^2, with rho+p/c^2=6.99344e-26 kg/m^3.

The companion-energy interpretation is a necessary source-sign pass only. It does not construct a solution to all Einstein equations, establish the required redshift/timing behavior, or prove stability and local conservation. It also does not remove the previous observed lens and brightness residuals.

## Connection to galaxy deposits

The retained intercepted galaxy fit requires

    F_capture = integral c*u_external(t) dt = 4.1618930494e16 J/m^2.

Under the additional constant-bath and shared-source identification,

    T_equivalent = F_capture/[c*u_required(0)]
                 = 9.331950518e8 years.

This uses the all-direction c*u convention of the existing capture cross-section calculation; it is not a one-sided planar flux and must not receive another factor of four. We have not assumed that this bath stays constant or that any galaxy has been capturing for that time. If the energy supply varies, its time integral must be calculated instead.

The redshift dependence in the table describes required source values along our past light cone. It cannot automatically be treated as a time history at one galaxy. A transport/spacetime solution must connect those locations and epochs. Streaming energy density is also distinct from permanent deposits in voids: this calculation does not introduce void capture, but the traveling bath's gravitational effect cannot be omitted under this branch's Einstein coupling.

## What this changes

There is now a quantitative conditional link between the optical fit's source requirement and the rotation fit's accumulated exposure. Neither an energy shortage nor sufficient stellar supply can be concluded from it. The next physical requirement is a source/transport history that produces the inferred bath and the capture exposure together, without double-counting energy. An arbitrary bath normalization is no longer enough once the all-companion focusing identification is made.

The independently known optical relation supplies the bridge; the source identity and fitted geometry remain hypotheses. The result is not evidence of a new interaction by itself.

Verification: symbolic observer expansion recovers the saved focusing limit exactly; SI/Mpc conversions are explicit; fitted input hashes are recorded. Reproduce with `python research_work/results/brightness-distance-consistency/optical-source.py`. Detailed results are in `optical-source-results.json`.

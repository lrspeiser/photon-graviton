# DRT-1: direction-resolved transport and a joint-energy test

**Exploratory, not adopted. No astronomical fits or regression grades changed.**

The source reviewed was main at `ee5a5f77d46371d05b5960fdf47d34cc3157ebb9`. The relevant static field functional is documented in `../code/eft_field_v24.py`. This experiment does not modify that code. The accompanying `drt1.py` was executed locally; its complete Git blob is `0d70b1795965cd5002c830db4f0c7e488eae2867`, verified against the uploaded source.

## Reproduce

Python 3.10+, NumPy and SciPy:

```sh
python research_work/results/hot-companion/drt1-transport/drt1.py --output-dir drt1-output
```

This generates a JSON report and nine crossing-stream snapshots at 160, 320 and 640 grid cells per dimension. Assertions test numerical identities and conservation, not acceptance as a gravity theory. Fixed seeds are included. No network, astronomical data or parameter fitting is required.

## Decisions

Retain direction-resolved transport for the hot, incoherent companion. It reproduces the existing hot scalar/vector kernels with an exact physical normalization and admits locally energy/momentum-conserving redirection. Do not silently replace the coherent cold sector with this intensity distribution.

Reject the tested minimal canonical kinetic completion of the static heat-repair energy: its freely variable counterstream occupations can lower the Hamiltonian without bound. This is a failure of the explicitly specified completion, not an instability run of the adopted QUMOND solver or a theorem against all possible actions.

## 1. Exact hot-moment normalization

Assume stationary transparent emission with power per volume `ell*k*rho` and common group speed `u`. The source luminosity gives

```
U = ell/(4*pi*u) integral k*rho/d^2 dV
F = ell/(4*pi) integral k*rho*n/d^2 dV
S_hot = (4*pi*G*u/ell) U
 g_hot = -(4*pi*G/ell) F
S_ex = (4*pi*G*u/ell) (U-|F|/u).
```

Here `n` points from source to observer. These are the existing hot kernels, not a new force prescription. Independent sums over 137 sources at 31 probes agreed to 3.27e-16 scalar and 4.69e-16 vector relative error. The identification is conditional on the specified luminosity and common-speed transparent transport; spectral transport requires frequency-dependent weighting.

## 2. Conservative redirection without ordinary-matter targets

Four directions support the reaction `(+x)+(-x) <-> (+y)+(-y)`. Both reaction directions are retained. Equal-energy packets stay on one wavenumber shell and conserve energy and vector momentum. The stipulated collision coefficient is not a derived microscopic cross section.

With `R=kappa*(f_xplus*f_xminus-f_yplus*f_yminus)`, the collision derivatives are `(-R,-R,+R,+R)`. An exact local update and exact lattice streaming evolve two colliding Gaussian packets.

| Collision coefficient | Transverse energy, n=160 | n=320 | n=640 |
|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 |
| 1 | 0.37583218 | 0.37598495 | 0.37602319 |
| 10 | 0.62802256 | 0.63164253 | 0.63263115 |

Maximum relative energy drift was 8.89e-16; momentum residual divided by initial energy was below 4.42e-16 in these normalized units. Twenty-four independent local ODE solutions agreed within 6.8e-13. A solitary beam does not turn through this particular reaction. Four angular bins and two spatial dimensions are a restricted test, not a general three-dimensional angular solution.

The central quadratic carrier has group speed u, phase speed u/2 and momentum per energy 2/u. The momentum ledgers retain that normalization.

## 3. A hot shell is bright inside without net flux or permanent trapping

For a transparent thin shell with radius R and luminosity P:

```
U(r) = P/(8*pi*u*R*r) log((R+r)/(R-r)), r<R
F(r) = 0.
```

The center limit is finite. At r/R=0, 0.5 and 0.9, the density relative to the center is 1, 1.098612 and 1.635799. A 512-node three-dimensional angular quadrature verifies the formula and flux cancellation. Switching the shell off leaves the center unchanged until R/u, after which its illumination disappears. At r=R/2 and t=R/u, 36.907% remains. Crossing intensity is not automatically a retained reservoir, and this density is not itself a force.

## 4. Gas: outside illumination and internal production differ

An exact reciprocal, elastic, two-direction slab has incoming intensities L and R, uniform production q, dimensionless thickness one and scattering optical depth tau. With u=1:

```
J'=q; U'=-tau*J
J0=[2*(L-R)-q*(1+tau/2)]/(2+tau)
J=J0+q*x
U=2*L-J0-tau*(J0*x+q*x*x/2).
```

For equal external illumination and q=0, U=2L everywhere at every tested optical depth: 0, 1, 3, 10 and 30. Scattering delays exit but also impedes entry. For internally produced energy above that bath, the mean density is `q*(1/2+tau/12)`, so opacity does enhance the internally supplied component.

**Gas still emits the model's cold baseline.** Equal-geometry diagnostic: M_gas/M_stars=20, stellar heat k=100, transparent stars. The stellar source supplies 50.5 normalized stored-energy units; the gas supplies 10 at tau=0, 26.667 at tau=10, 50.5 at tau=24.3 and 60 at tau=30. These are stipulated geometries/opacities, not measured cluster properties or lensing predictions.

The slab medium is fixed, with reaction momentum recorded. Separate 40,000 exact two-body collision checks at four target masses conserve combined energy/momentum while allowing finite recoil and frequency changes; they are not a full moving-gas transport simulation.

## 5. Dwarf and moving-source controls

Retaining a dwarf's tagged intensity under transparent superposition does not remove the external-field effect in the unchanged nonlinear force mapping. With a=1, internal g=1e-4, no exponential suppression and external/internal field=100, angular averaging gives 9.2409% of the isolated internal response. These are normalized controls, not a Draco fit. Directional transport alone does not repair this force response.

Source-velocity inheritance makes ballistic source-centered memory kinematically exact, but energy and momentum must transform too. For a quadratic packet with velocity V+u*n, mean laboratory energy is `(1+V^2/u^2)` times rest emission energy. At V/u=1.3 that is 2.69, with 62.8% in the bulk-motion term. This is a frame/energy-accounting requirement, not a computed drag rate or a complete source ledger. A newly changed monochromatic component travels only 34.658 kpc in 200 Myr relative to this transport frame.

## 6. Conditional join of the old spectrum and quadratic carrier

SGM-3 separately tested `w(z)=exp(-z/b)/sqrt(pi*b*z)`, b=pi/16, z=omega*L/u, and the carrier D=uL/2. Identifying this gate/source spectrum with the outgoing carrier is an ADDITIONAL hypothesis. Under it, group speed is `u*sqrt(2z)`.

The stationary density integrand `w/u_group` scales as 1/z near zero frequency. A source of finite age T supplies only frequencies `z >= zmin=(r/u/T)^2/2`. Therefore

```
U / [P/(4*pi*u*r^2)] = E1(zmin/b)/sqrt(2*pi*b).
```

At 30 kpc, ages 1, 5 and 13 Gyr give density ratios 1.86412, 4.69745 and 6.41564. These are energy-density ratios, not acceleration enhancements; the ages are illustrative source lifetimes, not an assumed expanding cosmology. Independent logarithmic-frequency integrals agree.

This same emitted spectrum has power-weighted mean speed 59.907 km/s, not 169 km/s. Only 1.3364e-15 of its power lies above 600 km/s. Slow-mode retention and age dependence are natural in this conditional join; a strong fast collision response requires a changed source spectrum or dispersion, not angular redirection alone.

## 7. The direct canonical energy completion fails

The fixed-source static gradient energy in the excess repair is

```
w_phi = g^3/(12*pi*G*a) - S_ex*g/(4*pi*G), g=|grad phi|.
```

Test an explicit minimal completion: add this energy to canonical free companion occupations with total free energy U. Using the normalization in section 1 and a=2*ell/u gives

```
H_local = g^3/(12*pi*G*a) + (1-2g/a)*U + (2g/a)*|F|/u.
```

On the monochromatic quadratic shell, `F/u=(u/2)*P`, where P is physical momentum density. No derivative at the cusp F=0 is assumed.

At balanced counterflow F=0, adding opposite packets preserves momentum but changes energy by `(1-2g/a)*delta_U`. It is negative for

```
g > a/2 = 3.148944695e-11 m/s^2.
```

This g is the companion potential gradient, NOT the ordinary field in the exponential screening factor. Holding g above threshold while allowing arbitrary occupation already gives no lower bound. Including and minimizing the cubic field energy does not cure it: with x=g/a and y=S_ex/a, normalized energy is `y+(2/3)*x^3-2*x*y`; its minimum is `y-(4/3)*y^(3/2)`. Seven independent numerical minimizations verify the expression.

Negative binding energy by itself is not a general exclusion. The problem here is unrestricted negative incremental energy without a stabilizing occupation constraint, reservoir, or backreaction. The static solver treats S_ex as supplied data; this audit does not claim that its numerical iterations show this instability.

A diagnostic positive saturation changes the force equation to `g^2=a*(Q+S*exp(-2g/a))`: at Q=S=a=1 it gives g=1.058487 instead of 1.414214. Saturation is not free: it changes the astronomical law and still does not resolve every symmetric-center issue. This alternative was not adopted or fitted.

## Bottom line

The transport picture now has exact normalization and executable conservative controls. It distinguishes counterflow from trapping, internal sources from external illumination, age from speed, and transport from the nonlinear force response. The most immediate attempted coupling to the existing static heat-repair energy fails its own boundedness requirement. Do not conceal that failure by evolving the transport without including its interaction energy.

External context, not evidence for the hypothetical companion: Weih et al., arXiv:2007.05718 (direction-resolved radiation transport); Savo et al., arXiv:1703.07114 and Pierrat et al., arXiv:1409.7229 (mean-path invariance under appropriate illumination). The calculations above are explicit toy-model results, not new observations.

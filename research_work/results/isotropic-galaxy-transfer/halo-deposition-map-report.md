# Fitted halo configurations and companion deposition targets

These are the exact numerical configurations of our restricted best fits, not uniquely measured halos. Some have poor stellar fits, extreme scale radii or parameter boundaries. They are not the published best-fit configurations from a full lens-image analysis.

## Halo definition

Known NFW density: rho_h(r)=rho_s/[x(1+x)^2], x=r/r_s. Known enclosed mass: M_h(<r)=A[ln(1+x)-x/(1+x)], A=4*pi*rho_s*r_s^3. Both r_s and rho_s specify the entire spherical fitting profile. r_s is its slope-transition radius, not the halo edge. The untruncated NFW total mass diverges logarithmically; A is a normalization, not total mass. A physical finite halo requires an outer truncation not measured here.

## Configurations under the standard benchmark geometry

| Galaxy | r_s kpc | rho_s Msun/kpc^3 | Halo mass inside 10 kpc, Msun | Stellar mass, Msun | Limits |
|---|---:|---:|---:|---:|---|
| J0037-0942 | 14.785 | 1.052e+07 | 4.8348e+10 | 7.8909e+11 | no bound hit |
| J1112+0826 | 63.901 | 6.2747e+06 | 2.0729e+11 | 4.2151e+11 | beta |
| J1204+0358 | 0.25163 | 1.0125e+11 | 5.5382e+10 | 2.8964e+11 | beta |
| J1402+6321 | 0.11176 | 1.135e+12 | 7.0002e+10 | 7.0651e+11 | beta |
| J1621+3931 | 924.95 | 1.0427e+05 | 5.9737e+10 | 6.1532e+11 | scale |
| J1630+4520 | 291.73 | 6.4713e+05 | 1.134e+11 | 5.6776e+11 | beta |

The table uses the standard geometry because it gave better NFW fits in the preceding comparison. This geometry is an alternative benchmark only. Our own geometry has separately fitted halo parameters in the JSON and should not be mixed with these numbers. J1402 still fits poorly; J1621 has a scale at its upper bound; J1204/J1402 have very compact fitted scale radii. No confidence intervals or realistic halo-population priors were obtained.

## Could our fixed deposit inventory occupy those regions?

Under our assumed ordinary gravitational coupling, rho_dep=u_dep/c^2. A deposited component matching the same density and relevant stresses would have the same gravitational effect regardless of its origin. That is a conditional equivalence, not evidence that capture produces that distribution.

| Galaxy | Original total deposit inventory, Msun | Target halo mass inside 5 Re, Msun | Inventory / target | Minimum inventory fraction moved to reproduce density inside 5 Re |
|---|---:|---:|---:|---:|
| J0037-0942 | 2.158e+12 | 3.2158e+11 | 6.711 | 33.096% |
| J1112+0826 | 6.3494e+11 | 1.6198e+12 | 0.392 | Insufficient fixed inventory |
| J1204+0358 | 1.8824e+11 | 7.0622e+10 | 2.665 | 34.537% |
| J1402+6321 | 1.4204e+12 | 9.9905e+10 | 14.218 | 45.344% |
| J1621+3931 | 1.4323e+12 | 1.2143e+12 | 1.180 | 36.520% |
| J1630+4520 | 9.5255e+11 | 1.5313e+12 | 0.622 | Insufficient fixed inventory |

This is a finite-aperture bookkeeping test, not a formation result. Five effective radii is a declared diagnostic region, not a measured halo boundary; portions extend beyond stellar constraints. The original profile is evaluated at the target physical radii without changing its existing inventory or scale. That comparison does not recalculate the one-third law under standard geometry. Both population cases and both halo geometries are retained, with masses at 0.5, 1, 2, 5 Re, 10 kpc, and the sphere at the Einstein radius. A sphere is not the projected cylinder measured by lensing.

Minimum moved mass is max(integrated positive density deficit, integrated positive density excess) inside 5 Re, assuming unrestricted rearrangement outside, a nonnegative reservoir, and enough total inventory. It is not an efficiency, migration probability or measured amount of movement. It applies only when the reported inventory suffices. Matching this finite radial interval does not reproduce the full line-of-sight lensing signal; mass beyond it and stellar mass must also be consistent.

## What the one-third law would need to produce

Our fixed law remains eta=X^(1/3)/(1+X^(1/3)), with rho0(r)=2*C_half*eta*J(r)/[1+(r/a)^2]^2. The Hill function and attenuation mathematics are known; their companion interpretation is proposed. The required redistribution is Delta_rho(r)=rho_h(r)-rho0(r). Positive Delta_rho means deposit more locally; negative means relocate some elsewhere. Cumulative target/original ratios in the JSON locate the mismatch without changing eta.

Equivalently a proposed effective capture-history model would need u_dep(r)=c^2*rho_h(r). If its local retained input power is q_ret(r), lossless accumulation requires integral q_ret(r,t)dt=c^2*rho_h(r), with transport and any work included. We have not measured that input or history. The fitted C_half inventory is not an independently established photon budget; passing its inventory test cannot establish sufficient real photon energy.

For a spherically symmetric redistribution flux j_E, energy continuity gives partial_t u+1/r^2 partial_r(r^2 j_E)=sources-sinks. Integrated over redistribution with no net local sources or sinks, 4*pi*r^2*integral j_E dt=-c^2*Delta_M(<r). This is known continuity mathematics: a positive required enclosed-mass change demands net inward energy transport. It does not supply a force law or a timescale for that transport.

Deposits also require a viable support/stress prescription. Placing positive energy at the fitted radii is not enough to ensure it stays there or produces the assumed lensing response. The model still lacks a demonstrated capture-and-migration law yielding these NFW-like targets. Merely naming the halo companion energy would relabel the inferred component without explaining it.

## Verification

The angular-capture density is integrated directly and compared with its analytic total capture-area expression, with relative tolerance 2e-4. Known spherical mass and E=Mc^2 conversions are used. SHA-256 hashes identify all three input results. No fit was rerun or silently selected on this inventory calculation. The prior geometry and stellar fit limitations continue to apply. All six goals remain open.

"""Render the executed potential construction and its scope from saved results."""
from pathlib import Path
import json
import math

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'results.json').read_text())
v=json.loads((HERE/'verification.json').read_text())
initial=json.loads((HERE/'initial-unsplit-results.json').read_text())
lines=[
'# A conservative three-dimensional extra potential for the Milky Way bar',
'',
'**The frozen empirical response now has an executed full-bar potential prototype.** Its radial, tangential and vertical accelerations are derivatives of one scalar potential. This supplies a missing forward-model component for the bulge comparison; it is not a fitted population of stellar orbits, a photon-derived gravity law or new holdout evidence.',
'',
f'The finest three-dimensional refinement changes the added force by at most **{100*r["force_refinement_maximum"]:.4f}%** on 60 fixed spatial probes (median {100*r["force_refinement_median"]:.6f}%). Separately refining the axisymmetric reference changes it by at most **{100*v["reference_force_change_max"]:.4f}%**; moving the outer source boundary from 80 to 100 kpc changes it by at most **{100*v["boundary_force_change_max"]:.4f}%**. These are finite numerical comparisons, not proofs of continuum accuracy or observational tolerances.',
'',
'## Equation and provenance',
'',
'Known QUMOND-style potential mathematics with the project\'s unchanged empirical response:',
'',
'`Q = A*(|grad(Phi_b)|/a_star)^(p-1)*grad(Phi_b)`,',
'',
'`Laplacian(Phi_extra) = div(Q)`, `g_extra = -grad(Phi_extra)`.',
'',
'The [QUMOND construction](https://arxiv.org/abs/0911.5464) is established literature. The parameters A=0.2422960666, p=0.4624587420 and a_star=7.2496089687e-10 m/s^2 are the earlier empirical galaxy fit. They are not newly fitted here, not unique first-principles laws, and not derived from the photon-loss rate, companion capture or stored energy. A shared fit can test this phenomenological gravity candidate while the physical photon bridge remains unresolved.',
'',
'The ordinary field contains the full published bar expansion through harmonic order 64, the existing stellar/gas disk expansion, nuclear components and softened central mass. No dark halo or synthetic companion ring/cap/shell is added. Its published component normalizations retain their model dependence. All components are evaluated in bar coordinates; a rotating-frame orbit model must supply the adopted pattern speed and observational orientation.',
'',
'## Three-dimensional Poisson construction',
'',
'We expand in orthonormal real cosine spherical harmonics Y_lm, retaining the even reflection symmetries of the adopted bar. All three components of Q enter. Define',
'',
'`q_lm(r) = integral Q_r*Y_lm dOmega`,',
'',
'`u_lm(r) = integral [Q_theta*partial_theta(Y_lm) + Q_phi*partial_phi(Y_lm)/sin(theta)] dOmega`.',
'',
'The projected source is `(1/r^2)*d_r(r^2*q_lm)-u_lm/r`. We integrate the spherical Poisson Green function by parts, using the existing stable radial power-kernel quadrature. This avoids estimating a noisy three-dimensional divergence on a Cartesian grid. Potential and force are evaluated from the same interpolated harmonic coefficients; the force is not assigned the algebraic vector -Q.',
'',
'These are applications of known harmonic analysis, Green functions and integration by parts. The implementation is a new project calculation, not a claimed new theorem. Analytically specified gradient sources with l=0,2,4 reproduce their known radial potential gradients with improving radial resolution. The check includes nonmonopole angular gradients and tests the sign and normalization independently of Milky Way data.',
'',
'## Recorded numerical repair',
'',
f'The initial unsplit angular construction changed forces by {100*initial["force_refinement_maximum"]:.3f}% between its two grids, failing the declared provisional 1% extra-force target. Its results remain in `initial-unsplit-results.json`; `run.py --unsplit` reproduces that calculation in `verification-unsplit-results.json`. This failure was not overwritten or reported as a passing result.',
'',
'The thin circular disk and the bar correction need different angular resolution. The revised calculation uses the exact linearity of the final Poisson solve:',
'',
'`Q_full = Q_axis + (Q_full-Q_axis)`.',
'',
'The first term receives an axisymmetric solution through l=128 with 512 angular quadrature nodes. The difference receives a full three-dimensional solution. Importantly, Q_axis is subtracted **before** projecting the nonlinear full-bar source: its nonzero mean difference remains in the three-dimensional correction. This is not silently replacing the true mean response by the response of the averaged baryons.',
'',
'The three correction grids use (radial nodes, polar nodes, azimuth nodes, harmonic order) of (192,64,96,24), (384,128,128,48), and (768,256,192,64). Azimuth sampling exceeds twice the largest retained m, avoiding a Nyquist-mode normalization ambiguity; numerical unit norms of the retained harmonics are checked. The reference extends to l=256 with 1024 angular nodes for its independent check. The disk representation in the input ordinary field remains the same throughout; it has not received a new convergence audit here.',
'',
'## Actual extra accelerations at representative positions',
'',
'Angles are relative to the bar x axis. Units for acceleration are (km/s)^2/kpc. These are predictions of the constructed extra potential only, not the total ordinary-plus-extra field, measured star velocities or circular-orbit speeds.',
'',
'| Radius (kpc) | Height (kpc) | Bar angle | Extra radial acceleration | Extra tangential acceleration | Extra vertical acceleration | Vector change from axisymmetric reference |',
'|---:|---:|---:|---:|---:|---:|---:|']
for row,check in zip(r['rows'],v['rows']):
    x,y,z=row['position_kpc'];R=math.hypot(x,y);phi=math.atan2(y,x)
    if abs(R-3)<1e-8 and z in [.1,.8]:
        ax,ay,az=row['fine_extra_acceleration'];ar=ax*math.cos(phi)+ay*math.sin(phi);ap=-ax*math.sin(phi)+ay*math.cos(phi)
        lines.append(f'| {R:g} | {z:g} | {math.degrees(phi):.0f} deg | {ar:.3f} | {ap:.3f} | {az:.3f} | {100*check["full_bar_extra_force_change_relative_axis"]:.3f}% |')
lines += ['',
'The bar-dependent correction is now a conservative gravitational force, unlike the earlier source-vector variation diagnostic. Its variation alone still does not predict a stellar velocity distribution. Stellar histories, the rotating bar, populations, distances and survey selection remain needed to compare the plane beneath the bulge with the upper and lower bulge.',
'',
'## Verification scope and remaining limits',
'',
f'- Direct finite differences of the potential agree with its force to {r["potential_derivative_error"]:.3g} maximum relative error on the probes. Reflection symmetry agrees to {v["reflection_error"]:.3g}.',
f'- All three checked maximum force changes are below the provisional 1% target: **{v["all_checked_force_changes_below_one_percent"]}**. A true value applies only to the listed spatial probes, numerical variations and fixed ordinary-matter representation.',
'- The 60 positions use radii 0.5,1,3,8,20 kpc, heights 0,0.1,0.8,1.5 kpc and bar angles 0,30,90 degrees. They are synthetic coordinates, not a held-out sample of measured stars. Domain and near-axis restrictions of the harmonic evaluator remain explicit.',
'- The main correction integration uses 0.001–80 kpc. The outer-boundary test uses 100 kpc and changes radial sampling as well as the boundary. It is a sensitivity check, not proof that unmodeled exterior material is negligible for every orbit. The inner boundary and ordinary disk expansion are not independently refined here.',
'- The uncut empirical force law has a divergent absolute monopole potential at infinity. We fix its additive gauge and integrate its mean radial force. We do not infer escape speeds, a finite total companion mass or a physical halo edge from this mathematical convention.',
'- The fixed bar is a test-particle source model. Its evolution and energetic backreaction are not solved. The construction is not a relativistic lensing completion or a local companion creation/capture action.',
'',
'## Next observational step',
'',
'The numerical field can now be taken into orbit-domain and Jacobi-conservation tests with the existing rotating-frame operator. A sufficiently broad shared orbit population must then be fitted to the measured position/velocity distributions with common selection and distance treatment, following the existing likelihood contract. The questionable cross-identifications remain a separate input issue; no orbit fit should turn them into evidence for stronger gravity.',
'',
'The previous frozen rotation, vertical, lensing and Cepheid predictions are preserved. No parameters or observed values were adjusted in this calculation, no old test roles were relabeled as fresh, and no holdouts were opened. The full scientific goal is not achieved by the solver checks. Photon conversion, observable time stretching, lossless transport, storage/response and the deferred total source budget still require a joint physical account.',
'',
'## Reproduction',
'',
'Run `run.py --unsplit` for the retained failed discretization, `run.py` for the split construction, `check.py` for the reference/boundary controls, and `report.py` for this report. Large coefficient caches stay in the ignored data cache; code, results, hashes and the preserved initial failure are tracked. `results.json` and `verification.json` contain the full 60-position predictions and check results.',
'']
(HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')

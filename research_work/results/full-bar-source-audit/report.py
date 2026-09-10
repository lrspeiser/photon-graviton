"""Report source-sign evidence separately from source accuracy and physical identity."""
from pathlib import Path
import json
import numpy as np
import run

HERE=Path(__file__).resolve().parent
d=json.loads((HERE/'results.json').read_text())
for path,expected in d['source_hashes'].items():
    assert run.digest(run.ROOT/path)==expected
rows=d['rows']
rho=np.array([r['rho_finer_Msun_kpc3'] for r in rows])
fine=np.array([r['rho_fine_Msun_kpc3'] for r in rows])
fd=np.array([r['rho_from_acceleration_steps'] for r in rows])
q=np.array([r['rho_from_requested_divQ_steps'] for r in rows])
checks={
    'fine_to_finer_source':abs(fine-rho)/abs(rho),
    'analytic_laplacian_vs_acceleration_divergence':abs(fd[:,-1]-rho)/abs(rho),
    'requested_divQ_vs_reconstructed_laplacian':abs(q[:,-1]-rho)/abs(q[:,-1]),
    'requested_divQ_step_change':abs(q[:,-1]-q[:,-2])/abs(q[:,-1]),
}
stats={}
for name,values in checks.items():
    index=int(np.argmax(values))
    stats[name]=dict(maximum_fraction=float(values.max()),median_fraction=float(np.median(values)),
                     worst_location={k:rows[index][k] for k in ['R_kpc','z_kpc','phi_rad']})
assessment=dict(samples=len(rows),negative_fine=int((fine<0).sum()),
                negative_finer=int((rho<0).sum()),negative_requested=int((q[:,-1]<0).sum()),
                checks=stats,global_positivity_proved=False,physical_deposit_law_derived=False,
                source_accuracy_1_percent_everywhere_on_grid=False,holdouts_opened=False,
                results_sha256=run.digest(HERE/'results.json'))
(HERE/'assessment.json').write_text(json.dumps(assessment,indent=2)+'\n',encoding='utf-8',newline='\n')
lines=[
    '# What source would produce our full-bar extra gravity?', '',
    '**All 240 sampled locations have positive equivalent source density**, both in the fine/finer reconstructed potentials and in a direct evaluation of the requested field-equation source. This is a limited compatibility result for a positive Newtonian source. It does not show that photons produce the source, establish positivity everywhere, or validate its detailed density distribution.', '',
    'This audit applies to the current conservative three-dimensional bar field used in the orbit integrations. The earlier source audit concerned different potential constructions; its positivity results could not automatically be transferred to this field.', '',
    '## Equations and their status', '',
    'The known Newtonian Poisson relation gives an equivalent source:', '',
    '`rho_extra = Laplacian(Phi_extra)/(4*pi*G) = -div(a_extra)/(4*pi*G)`.', '',
    'This is the density that would reproduce the extra potential if ordinary Newtonian sourcing applied. It is not automatically deposited energy divided by c squared. That identification would be an additional assumption about how companion energy gravitates; a spacetime-response interpretation needs its own field equation and stress/energy accounting.', '',
    'The implemented field uses the known [QUMOND construction (Milgrom 2010)](https://arxiv.org/abs/0911.5464) with the project\'s empirical power-law coefficients. Here Q denotes only the extra response; ordinary gravity is added separately:', '',
    '`Q = A (|grad(Phi_b)|/a_star)^(p-1) grad(Phi_b)`',
    '`Laplacian(Phi_extra) = div(Q)`.', '',
    'The direct div(Q) calculation tests the intended source, while the Laplacian of the cached potential tests the source actually represented by our numerical solver. These need not agree perfectly at finite resolution. Neither expression derives a photon-conversion rate, capture cross section or occupied storage capacity.', '',
    'For each real spherical harmonic of degree l, the known Laplacian identity used here is:', '',
    '`Laplacian[f_lm(log r) Y_lm] = [f_lm\'\' + f_lm\' - l(l+1) f_lm] Y_lm / r^2`.', '',
    'The primes refer to log-radius derivatives. This evaluates the second derivatives of the same spline potential used for the force. An independent finite difference of its acceleration checks the sign and normalization. These are established mathematical identities, not claimed novel formulas.', '',
    '## Grid and results', '',
    'The grid uses R = 0.5, 1, 2, 3, 5, 8, 12 and 20 kpc; z = 0, 0.1, 0.5, 1, 2 and 4 kpc; and five bar angles from 0 to pi/2. These are numerical probes, not observed stars. The retained even harmonics impose reflection symmetries; this grid cannot test real-galaxy asymmetries. No dark halo, new parameter adjustment or holdout observation is introduced.', '',
    f"The reconstructed equivalent density ranges from {rho.min()/1e9:.6f} to {rho.max()/1e9:.6f} solar masses per cubic parsec. The directly requested source ranges from {q[:,-1].min()/1e9:.6f} to {q[:,-1].max()/1e9:.6f}. These are conditional model outputs, not measurements of companions.", '',
    '| Numerical comparison | Largest fractional difference | Median fractional difference |',
    '|---|---:|---:|',
]
for name,values in stats.items():
    lines.append(f"| {name.replace('_',' ')} | {100*values['maximum_fraction']:.5f}% | {100*values['median_fraction']:.5f}% |")
lines += ['',
    'The largest fine/finer density change is 5.414% at R=1 kpc, z=4 kpc, on the bar axis. The largest intended-versus-reconstructed source discrepancy is 1.091% at R=12 kpc, z=0.5 kpc, on the bar axis. The latter comparison uses the same cached ordinary-matter field and is not an independent assessment of its mass accuracy. The axisymmetric reference, ordinary-matter model and outer boundary were not further refined in this audit.', '',
    'Finite differences use steps 0.004, 0.002 and 0.001 kpc. The final requested-source step change is at most 0.00414%; the cached-potential Laplacian and acceleration-divergence calculations agree to 0.00136%. Thus the larger source discrepancies are not removed merely by reducing this finite-difference step. More angular/radial field resolution and ordinary-matter uncertainty would be needed before claiming percent-level source accuracy throughout the grid.', '',
    '## What this changes—and what remains to be derived', '',
    'The sign check does not reveal a negative-source obstruction on this grid. It therefore leaves open the conditional possibility that a positive source produces the proposed extra field. However, a positive source is only a target: we still need a common transport/capture/support law that predicts this spatial distribution, rather than assigning the required density after fitting the gravity.', '',
    'The previously demonstrated amplitude degeneracy also remains: changing the photon-loss scale and compensating with a free gravity amplitude preserves the rotation predictions. Positive reconstructed density does not remove that freedom. A photon-origin claim requires an independently constrained response or capture law.', '',
    'Practical next steps are to keep the density-resolution limitation visible, specify the physical response that maps deposited energy to the scalar potential, and test whether that one law predicts radial, off-plane and lensing behavior. The total cosmic photon-supply budget remains deferred, not passed. No new observational agreement or global source admissibility is claimed here.', '',
    '## Reproduction', '',
    'Run `run.py`, then `report.py`. Results retain every sampled location, all derivative steps, both potential resolutions and the direct source. Input hashes are checked before and after calculation and again during reporting. A quadratic-potential control verifies divergence sign and normalization. No catalog, empirical coefficient, orbit cache or sample role is changed.', '',
]
(HERE/'report.md').write_text('\n'.join(lines),encoding='utf-8',newline='\n')
print(json.dumps(assessment,indent=2))

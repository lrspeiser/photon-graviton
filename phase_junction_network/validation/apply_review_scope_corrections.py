#!/usr/bin/env python3
"""Explicit, idempotent maintenance of claim labels; numerical models unchanged."""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]
MARK = 'Validation scope correction — 2026-09-20'
FROZEN = ROOT / 'microscopic/fermionic_multipair_vacuum_results.json'
EXPECTED = '9fc66061e31fe53c4b8a016c5bda9eda8a5f00e0d7efa044314a443a1b5d9400'
assert hashlib.sha256(FROZEN.read_bytes()).hexdigest() == EXPECTED

notes = {
    'README.md': 'The current finite electromagnetic evidence is a static 3D anchor, reduced-chain dynamics and a separate small full-cube control—not an established finite-spin 3D photon phase. The natural chain embedding fails an exact leakage test. Stage 3I reproduces, but its old cutoff comparison rebuilds Hamiltonians; use the new fixed-Hamiltonian test. Shared free-slab matter loops leave an uncancelled uniform frame-source response. Local nonlinear gravity is not closed.',
    'architecture_audit.md': 'This historical audit predates the integration tests. The new audit distinguishes reduced chains from a proven 3D phase, rebuilt cutoff models from fixed-Hamiltonian errors, source curvatures from physical speeds, and conditional nonlinear targets from microscopic closure.',
    'microscopic/README.md': 'The latest authoritative scope is the validation audit. Earlier Stage 6B chains are reduced models, not proven exact 3D sectors; the shared-speed root is a reduced-model calibration. Full nonlinear and interacting continuum claims remain open.',
    'microscopic/finite_em_dynamics.md': 'The chain Hamiltonian is exactly diagonalized as a reduced model, but no exact embedding in the finite-spin 3D parent has been established. The natural sheet embedding now fails on L=2,3,4. Two copies are constructed, not independently discovered 3D polarizations. The cubic tensor is a conditional Gaussian control. All numerical tables below are retained at reduced-model scope.',
    'microscopic/seagull_completion_stage3e.md': 'The speed crossing selects a parameter by matching a reduced-chain photon fit to a harmonic frame target. It is neither an independently predicted physical speed equality nor established equality of the finite regulators. The unit diagonal completion is a selected candidate rule; its symmetry protection remains unestablished.',
    'microscopic/fermionic_multipair_vacuum_stage3i.md': 'The full 6,336-state result is independently reproduced and retained. The historical 1.3965% comparison rebuilds the completed Hamiltonian after each cutoff; it is not a fixed-Hamiltonian truncation error. The new Q_N H_full Q_N comparison gives loop-correction errors of -52.7275% (one pair), -0.100742% (two), and -0.00004966% (three). These are percentages of the small correction, not the total gap. The 0.1243% number is a finite-momentum pole shift, not an infrared speed measurement.',
    'microscopic/nonlinear_volume_reduction_stage5a.md': 'The homogeneous identity remains valid under its stated target kinetic form. Appending determinant/trace conditions locally changes the old constraint classification: at nonzero momentum four conditions are second class and two remain first class, leaving two physical modes. Full nonlinear lattice closure and the microscopic origin of this kinetic geometry remain open.',
}
changed = []
for relative, note in notes.items():
    path = ROOT / relative
    text = path.read_text()
    if MARK in text:
        continue
    title, rest = text.split('\n', 1)
    link = 'README.md' if relative.startswith('validation/') else ('../validation/README.md' if relative.startswith('microscopic/') else 'validation/README.md')
    banner = f'\n> **{MARK}.** {note} See the [authoritative audit]({link}) and its executable results.\n'
    text = title + '\n' + banner + rest
    if relative == 'README.md':
        text = text.replace('contains a finite pure-gauge photon phase,', 'contains a finite static electromagnetic anchor with reduced-chain dynamics,')
        text = text.replace('A finite detuned region supports two transverse photon-like branches.', 'The reduced chains support linearly dispersing photon-like branches over the tested detuned region; two copies are constructed. A finite-spin 3D dynamical phase is not established.')
        text = text.replace('The accepted finite candidate root is', 'The selected reduced-chain/harmonic-frame matching root is')
        text = text.replace('The full result changes the one-pair photon correction by only', 'Across the historical rebuilt-cutoff models, the full result changes the one-pair photon correction by')
        text = text.replace('Reject the candidate if matter loops require a photon-only or frame-only counterterm.', 'Require microscopic source completion and controlled infrared kinetic matching; do not tune a photon-only or frame-only counterterm to manufacture agreement.')
    if relative == 'microscopic/finite_em_dynamics.md':
        text = text.replace('PASS at the pure-gauge dynamical stage.', 'PASS at reduced-chain dynamical scope, with a separate finite-cube control.')
        text = text.replace('## 3. Exact transverse propagation block', '## 3. Reduced-chain model (exact 3D embedding not established)')
        text = text.replace('For propagation along one principal axis, Gauss reduction leaves two independent transverse electric-flux fields. One polarization is the exact finite chain', 'The proposed principal-axis reduction is represented by the finite chain below. Its identification with a sector of the full 3D parent is unproven; it is now tested only as a reduced model')
        text = text.replace('whose exact principal-axis transverse sectors support', 'whose reduced principal-axis chain models support')
        text = text.replace('This completes the pure-gauge dynamical step requested after Stage 6A.', 'This completes the reduced-chain calculation, not the finite-spin 3D dynamical phase requirement.')
    if relative == 'microscopic/fermionic_multipair_vacuum_stage3i.md':
        text = text.replace('## 6. Convergence of the pair truncation', '## 6. Historical comparison across rebuilt cutoff Hamiltonians')
        text = text.replace('## 11. Finite common-cone diagnostic', '## 11. Finite-momentum pole-shift diagnostic (not an infrared cone)')
        text = text.replace('The shared candidate must be rejected if matter dressing preserves the electromagnetic Ward identity but independently shifts the photon and frame cones.', 'The shared candidate requires a controlled infrared comparison after both source and kinetic normalizations are derived. Unequal raw shifts at one large lattice momentum do not alone establish unequal limiting speeds.')
    if relative == 'microscopic/nonlinear_volume_reduction_stage5a.md':
        text = text.replace('retain the existing scalar and vector first-class constraints;', 'reclassify the existing scalar and vector constraints after adding local volume conditions; do not assume all four remain first class;')
        text = text.replace('At the required spin-2 value', 'At the selected target kinetic value')
    path.write_text(text)
    changed.append(relative)

path = ROOT / 'microscopic/finite_em_dynamics_impl.py'
text = path.read_text()
code_marker = '# Validation scope: chains are reduced models, not established exact 3D sectors.'
if code_marker not in text:
    replacements = {
        'an exact principal-axis transverse block of the same finite link operators;': 'a reduced principal-axis chain built from the same finite link amplitudes;',
        'Build one exact transverse flux block in the zero-total-flux sector.': 'Build one reduced chain in its zero-total-flux sector.',
        'two_polarization_first_gap_degeneracy': 'constructed_copy_count',
        'second_identical_transverse_copy_gives_twofold_gap': 'two_identical_copies_constructed_not_3d_degeneracy_measurement',
        'Stage 6B PASS: finite detuned pure-gauge transverse dynamics has a stable z=1 interval; dynamical matter and interacting QED remain open': 'Stage 6B reduced-chain PASS; finite-spin 3D phase, interacting matter and QED remain open',
        'zero total transverse flux; longitudinal link removed by Gauss law': 'zero total chain flux; not a proof of 3D Gauss reduction',
        'two identical independent transverse blocks for propagation along a principal axis': 'two identical reduced-chain copies; 3D embedding not established',
        'an exact finite principal-axis transverse block of the Stage-6A spin-link Hamiltonian': 'an exactly diagonalized reduced chain using Stage-6A finite-link amplitudes',
        'the cubic long-wavelength symbol has rank two and L^-2 anisotropy recovery': 'a separately assumed Gaussian cubic symbol has rank two and L^-2 anisotropy recovery',
    }
    for old, new in replacements.items():
        if old not in text:
            raise RuntimeError('Source anchor changed: '+old)
        text = text.replace(old, new)
    text = text.replace("'claim_boundary': {'established':", "'physical_3d_phase_established': False, 'claim_boundary': {'established':")
    text = text.replace("default=Path('phase_junction_network/microscopic/finite_em_dynamics_results.json')", "default=Path('phase_junction_network/microscopic/finite_em_dynamics_revalidated.json')")
    text += '\n' + code_marker + '\n'
    path.write_text(text)
    changed.append('microscopic/finite_em_dynamics_impl.py')
assert hashlib.sha256(FROZEN.read_bytes()).hexdigest() == EXPECTED
print('Scope corrections: '+', '.join(changed))

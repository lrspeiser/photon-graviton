#!/usr/bin/env python3
"""One-time, guarded source migration for the 2026-09-20 validation review.

Run only as explicit maintenance. Ordinary validation never edits source files.
The archived Stage-3I results are intentionally not changed.
"""
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'microscopic/check_fermionic_multipair_vacuum.py'
text = SOURCE.read_text()
marker = '# Validation baseline: reject unidentifiable logarithmic fits.'
if marker in text:
    print('Baseline migration already applied; no source changes.')
else:
    blob = hashlib.sha1(b'blob ' + str(len(text.encode())).encode() + b'\0' + text.encode()).hexdigest()
    if blob != 'c0c18687df6785a5ba010cc826e8e138bc83e4a5':
        raise RuntimeError('Unexpected source revision; reconcile rather than overwrite: ' + blob)
    start = text.index('def logarithmic_fit(xs, ys)')
    stop = text.index('\n\ndef build_report', start)
    replacement = '''def logarithmic_fit(xs, ys) -> dict[str, float]:
    # Validation baseline: reject unidentifiable logarithmic fits.
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape or x.size < 2:
        raise ValueError("A log fit requires at least two paired samples")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("Log-fit samples must be finite")
    if np.any(x <= 0) or np.any(y <= 0):
        raise ValueError("Log-fit samples must be positive")
    lx, ly = np.log(x), np.log(y)
    centered = lx - lx.mean()
    denominator = float(centered @ centered)
    if denominator <= 1.0e-24:
        raise ValueError("A log-fit exponent needs distinct abscissae")
    power = float(centered @ (ly - ly.mean()) / denominator)
    prefactor = float(np.exp(ly.mean() - power * lx.mean()))
    return {"power": power, "prefactor": prefactor}
'''
    text = text[:start] + replacement + text[stop:]
    old = '        if gap in (0.75, 1.0, 1.5):'
    if text.count(old) != 1:
        raise RuntimeError('Polarization scan anchor changed')
    text = text.replace(old, '        if quick or gap in (0.75, 1.0, 1.5):')
    # The full result remains an immutable historical record. New CLI runs have
    # a different default destination and state their cutoff interpretation.
    text = text.replace('"fermionic_multipair_vacuum_results.json"\n        ),',
                        '"fermionic_multipair_vacuum_revalidated.json"\n        ),')
    anchor = '    if not report["stage_pass"]:'
    note = '''    report["cutoff_semantics"] = (
        "Rebuilt completed-move Hamiltonians at each pair cutoff; not a "
        "fixed-Hamiltonian truncation error. See validation/check_fixed_cutoff.py."
    )
    report["not_tested"] = (
        ["full_four_pair_convergence", "coupling_scan", "full_sector_counts"]
        if quick else []
    )
    if quick:
        for key in ("multipair_series_converges", "coupling_scan_monotonic",
                    "neutral_fock_sector_counts"):
            report["checks"][key] = None
    report["stage_pass"] = all(
        value for value in report["checks"].values() if value is not None
    )
'''
    if text.count(anchor) != 1:
        raise RuntimeError('Report anchor changed')
    text = text.replace(anchor, note + '\n' + anchor)
    text = text.replace('                    if not passed\n', '                    if passed is False\n')
    SOURCE.write_text(text)
    print('Applied guarded fitting repair; legacy frozen output untouched.')

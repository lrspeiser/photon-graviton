"""Run current research diagnostics in a fresh directory, preserving saved evidence."""
from pathlib import Path
import argparse
import ast
import hashlib
import gzip
import json
import os
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'research_work/results'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot():
    manifest = json.loads((ROOT / 'SNAPSHOT_MANIFEST.json').read_text(encoding='utf-8'))
    mismatches, missing, present = [], [], {}
    for item in manifest['files']:
        p = ROOT / item['path']
        compressed = p.with_name(p.name + '.gz')
        if not p.is_file() and compressed.is_file():
            present[item['path']] = hashlib.sha256(gzip.decompress(compressed.read_bytes())).hexdigest()
            if present[item['path']] != item['sha256']:
                mismatches.append(item['path'])
        elif not p.is_file():
            missing.append(item['path'])
        else:
            present[item['path']] = sha(p)
            if present[item['path']] != item['sha256']:
                mismatches.append(item['path'])
    if mismatches:
        raise RuntimeError(f'Original snapshot differs: {mismatches}')
    expected_missing = {'unified_paper_v9/' + name for name in
                        ['build.py', 'design_tokens.json', 'math_format.py',
                         'verification.json', 'version8_original.docx']}
    if set(missing) != expected_missing:
        raise RuntimeError(f'Unexpected missing originals: {missing}')
    return present, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', action='store_true', help='Also rerun the three historical baseline scripts in an isolated copy.')
    parser.add_argument('--output-dir', type=Path, help='Fresh output directory; existing nonempty directories are refused.')
    args = parser.parse_args()
    out = (args.output_dir or ROOT / 'research_work/generated' /
           (time.strftime('%Y%m%d-%H%M%S') + '-' + uuid.uuid4().hex[:6])).resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        parser.error('Output directory must be new or empty.')
    original, missing = snapshot()
    saved = {p: sha(p) for p in RESULTS.rglob('*') if p.is_file()}
    for p in list(RESULTS.rglob('*.py')) + list((ROOT / 'research_work/energy').glob('*.py')):
        ast.parse(p.read_text(encoding='utf-8'), filename=str(p))
    out.mkdir(parents=True, exist_ok=True)
    jobs = [
        ('energy/compare_no_loss_catalog.py', []),
        ('energy/ledger_solver.py', ['--output', str(out / 'energy/ledger-checks.json')]),
        ('transport/check_transport.py', []),
        ('radial-capture/run_radial_capture.py', []),
        ('radial-capture/check_external_capture.py', []),
        ('microphysics/check_microphysics.py', []),
        ('frequency-transfer/check_frequency_transfer.py', []),
        ('scalar-wave/check_scalar_wave.py', []),
        ('gravity-response/check_gravity_response.py', []),
        ('data-audit/audit_observables.py', []),
        ('data-audit/enrich_observable_registry.py', []),
        ('conversion-first/run_conversion_first.py', []),
        ('interaction-rate/check_interaction_rate.py', []),
        ('matter-assisted/check_matter_assisted.py', []),
        ('collective-response/check_collective_response.py', []),
        ('oscillator-response/check_oscillator_response.py', []),
        ('soft-graviton/check_soft_graviton.py', []),
        ('capture-storage/check_capture_storage.py', []),
        ('deposit-support/check_deposit_support.py', []),
        ('capture-momentum/check_capture_momentum.py', []),
        ('isotropic-capture/check_isotropic_capture.py', []),
        ('broadband-capture/check_broadband_capture.py', []),
        ('source-receiver/check_source_receiver.py', []),
        ('matched-wave/check_matched_wave.py', []),
        ('thermal-conversion/check_thermal_conversion.py', []),
        ('background-replenishment/check_background_replenishment.py', []),
        ('thermalizer/check_thermalizer.py', []),
        ('slow-companion/check_slow_companion.py', []),
        ('dispersive-companion/check_dispersive_companion.py', []),
        ('electromagnetic-balance/check_electromagnetic_balance.py', []),
        ('companion-bath/check_companion_bath.py', []),
        ('cumulative-time/check.py', []),
        ('companion-backreaction/check.py', []),
        ('weak-signal-timing/check.py', []),
        ('environmental-screening/check.py', []),
        ('generated-wave-access/check.py', []),
        ('source-timescale/check.py', []),
        ('multiple-regions/check.py', []),
        ('radiation-train/check.py', []),
        ('sustained-illumination/check.py', []),
        ('spherical-propagation/check.py', []),
        ('companion-self-binding/check.py', []),
        ('bound-cloud-exchange/check.py', []),
        ('bound-pair-production/check.py', []),
        ('pair-production-balance/check.py', []),
        # Regenerate the latest canonical diagnostics into this run and compare with the archive.
        ('companion-extensions/mond-inventory.py', []),
        ('companion-extensions/mond-cross-scale.py', []),
        ('companion-extensions/coma-inverse.py', []),
        ('capture-to-orbit/checks.py', []),
        ('propagation-field/pf1.py', []),
        ('propagation-field/brightness.py', []),
        ('radiation-polarized-gravity/checks.py', []),
        ('collective-reservoir/checks.py', []),
        ('clock-completion/cc1.py', []),
        ('clock-gradient/timelike.py', []),
        ('supported-reservoir/checks.py', []),
        ('gravitational-focusing/checks.py', []),
        ('companion-supply/cc2a.py', []),
        ('companion-formation/checks.py', []),
        ('radiation-budget/rc1.py', []),
        ('companion-source/checks.py', []),
        ('donor-companion/checks.py', []),
        ('shared-field-bridge/checks.py', []),
        ('path-memory/checks.py', []),
        ('path-memory/pm2a.py', []),
        ('path-memory/pm3.py', []),
        ('path-memory/rut1.py', []),
        ('path-memory/rut3.py', []),
        ('path-memory/rut4_checks.py', []),
        ('path-memory/rut4r.py', []),
        ('path-memory/rut5.py', []),
        ('path-memory/rut6.py', []),
        ('path-memory/rut7_checks.py', []),
        # The owner's corrected annulus sampler (2f1d5ed, bde7a7c): a numerical-verification job in its own
        # right, outside the frozen results tree, run as a module from the repository root.
        ('module:research_work.annulus_sampling.checks', []),
        # CL-1, the written-track response at cluster scale and for light (path-memory/protocol-cl1.md): the
        # nine numerical gates with their negative controls and a regression anchor.
        ('path-memory/cl1_checks.py', []),
        # CL-2 stage 1, one written-field response tested jointly on galaxies, lenses and clusters
        # (path-memory/protocol-cl2.md and its two amendments): eight gates with controls and a regression anchor.
        ('path-memory/cl2_checks.py', []),
        # CL-F1: spherical cluster optics, analytic projection and image-position gates.
        ('path-memory/cluster_lensing_checks.py', []),
        # TF-1, the transverse coupling of the converted field (experiments/transverse_coupling/protocol.md): the fast
        # exact-property gates with their controls, anchored to the archived exact stage.
        ('../experiments/transverse_coupling/tf1_checks.py', []),
        # CL-2 stage 2 (path-memory/protocol-cl2-stage2.md and its amendments): the disk-force gates, the certificate of the
        # archived galaxy optimum recomputed from the inputs, and the SZ correlation matrices, anchored to cl2s2-results.json.
        ('path-memory/cl2s2_checks.py', []),
        # NL-1, the root spectrum (path-memory/protocol-nl1.md, amendment 1): the archived optimum's certificate recomputed
        # from the inputs and the local law's reproduction of the archived reference, anchored to nl1-results.json.
        ('path-memory/nl1_checks.py', []),
        # NK-1, shell footprints (path-memory/protocol-nk1.md): the kernel gates and per-galaxy anchors on three galaxies.
        ('path-memory/nk1_checks.py', []),
        # RW-1, the whirlpool beyond the rut (path-memory/protocol-rw1.md, amendment 1): the kernel gates and the archived
        # G-whole strengths of three galaxies recomputed from the inputs, anchored to rw1-results.json.
        ('path-memory/rw1_checks.py', []),
        # RUT-1 stage 9, the corrected response calculation (path-memory/protocol-rut9.md): the boundary terms
        # of the sharp truncation, every ring-limit annulus resolved on its own, and B13's own eigenfunction in
        # the time domain. Its status is its own and does not touch stage 8's.
        ('path-memory/rut9_checks.py', []),
        # RUT-1 stage 10, the seeded full-state experiment (path-memory/protocol-rut10.md): replays the exact
        # preparation of one seeded run, bodies and field, bit for bit, and re-derives every fit, gate and reading.
        ('path-memory/rut10_checks.py', []),
        # LAST, and expected to exit non-zero: by the owner's ruling on 3a80fec stage 8's numerical verification
        # status stays FAILED in its archive (the ring-limit gate, as declared and as amended). The job says so.
        ('path-memory/rut8_checks.py', []),
    ]
    if args.baseline:
        jobs.append(('baseline/run_baseline_001.py', []))
    # Bytecode caches under research_work/results are regenerated artifacts; a job refreshing one
    # would trip the saved-results integrity check, so jobs never write them.
    env = {**os.environ, 'PHOTON_GRAVITON_RESULTS': str(out), 'PYTHONUTF8': '1',
           'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'PYTHONDONTWRITEBYTECODE': '1'}
    records = []
    for name, extra in jobs:
        start = time.monotonic()
        target = (['-B', '-m', name[len('module:'):]] if name.startswith('module:')
                  else [str(RESULTS / name)])
        result = subprocess.run([sys.executable, '-X', 'utf8', *target, *extra],
                                cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=3000)
        (out / (name.replace('/', '__').replace(':', '__') + '.log')).write_text(
            result.stdout + result.stderr, encoding='utf-8')
        records.append({'script': name, 'exit_code': result.returncode, 'seconds': time.monotonic() - start})
        print(('PASS ' if result.returncode == 0 else 'FAIL ') + name, flush=True)
        if result.returncode:
            print(result.stderr, file=sys.stderr)
            break
    current, _ = snapshot()
    preserved = current == original and all(p.is_file() and sha(p) == digest for p, digest in saved.items())
    report = {'jobs': records, 'original_files_verified': len(original), 'missing_historical_files': missing,
              'originals_and_saved_results_preserved': preserved,
              'scope': 'Numerical diagnostics and integrity; not observational validation of the theory.'}
    (out / 'verification-summary.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(f'Results: {out}', flush=True)
    if not preserved or len(records) != len(jobs) or any(r['exit_code'] for r in records):
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

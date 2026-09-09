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
    ]
    if args.baseline:
        jobs.append(('baseline/run_baseline_001.py', []))
    env = {**os.environ, 'PHOTON_GRAVITON_RESULTS': str(out), 'PYTHONUTF8': '1',
           'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    records = []
    for name, extra in jobs:
        start = time.monotonic()
        result = subprocess.run([sys.executable, '-X', 'utf8', str(RESULTS / name), *extra],
                                cwd=ROOT, env=env, capture_output=True, text=True, encoding='utf-8', timeout=3000)
        (out / (name.replace('/', '__') + '.log')).write_text(result.stdout + result.stderr, encoding='utf-8')
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

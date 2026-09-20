"""SE-LR3 archive reconstruction independent of its evolution/ledger functions."""
import hashlib
import json
from pathlib import Path
import numpy as np
from audit_local_rotor import reconstruct, loop

ROOT = Path(__file__).resolve().parent


def main():
    out = ROOT/'local-transfer-v1'
    summary = json.loads((out/'summary.json').read_text())
    controls = json.loads((out/'controls.json').read_text())
    manifest = json.loads((out/'manifest.json').read_text())
    checks = [len(controls['cases']) == 48, controls['passed'] == all(r['passed'] for r in controls['cases'])]
    checks += [hashlib.sha256((ROOT/f).read_bytes()).hexdigest() == digest for f, digest in manifest['hashes'].items()]
    for row in controls['cases']:
        scale = max(1., abs(row['analytic']))
        errors = [abs(v-row['analytic'])/scale for v in row['numerical']]
        diff = abs(row['numerical'][0]-row['numerical'][1])/scale
        checks.append(row['passed'] == (max(errors+[diff]) < 2e-7))
    maximum = 0.
    derived = {}
    for row in summary['runs']:
        archive = np.load(out/(row['name']+'.npz')); trace = archive['trace']; h = 8/row['n']
        checks.append(row == json.loads((out/(row['name']+'.json')).read_text()))
        for name, t in (('initial', 0), ('final', -1)):
            state = archive[name]
            result = reconstruct(state, h, row['epsilon'])
            shifted = state.copy(); shifted[3] -= row['direct']*state[0]
            result[:3] = reconstruct(shifted, h, row['epsilon'])[:3]
            field = state.copy(); field[2:] = 0
            field_j = reconstruct(field, h, row['epsilon'])[3:6]
            error = float(max(np.max(abs(result-trace[t, 1:10])), np.max(abs(field_j-trace[t, 11:14]))))
            maximum = max(maximum, error); checks.append(error < 1e-12)
        checks.append(abs(field_j[2]/trace[0, 6]-row['A_angular_fraction']) < 1e-12)
        checks.append(abs(result[1]/trace[0, 1]-row['A_energy_fraction']) < 1e-12)
        edrift = np.max(abs(trace[:, 1]-trace[0, 1]))/trace[0, 1]
        jdrift = np.max(np.linalg.norm(trace[:, 4:7]-trace[0, 4:7], axis=1))/np.linalg.norm(trace[0, 4:7])
        edge = np.max(trace[:, 10])/trace[0, 1]
        checks.extend((abs(edrift-row['energy_drift']) < 1e-12, abs(jdrift-row['angular_drift']) < 1e-12, abs(edge-row['edge_fraction']) < 1e-12))
        for item in row['loops']:
            for key, count in (('coarse', 128), ('fine', 256)):
                checks.append(abs(loop(archive['final'][0], h, item['radius'], count)-item[key]) < 1e-12)
        passed = bool(np.isfinite(trace).all() and np.isfinite(archive['final']).all() and edrift < 1e-5 and jdrift < .01 and edge < 1e-5
                      and max(abs(r['coarse']-r['fine']) for r in row['loops']) < 1e-5)
        checks.append(passed == row['passed'])
        derived[row['name']] = (loop(archive['final'][0], h, 1.5, 256), field_j[2])
    replay = float(np.max(abs(np.load(out/'curl.npz')['final']-np.load(ROOT/'local-rotor-evolution-v1'/'rotating.npz')['final'])))
    checks.append(summary['comparisons'][0]['error'] == replay)
    checks.append(summary['comparisons'][0]['passed'] == (replay < 1e-12))
    for row, limit in zip(summary['comparisons'][1:], (.001, .05)):
        values = [abs(coarse-fine)/max(abs(fine), 1e-10) for coarse, fine in zip(derived['combined'], derived[row['name']])]
        checks.extend((abs(values[0]-row['circulation_relative_difference']) < 1e-12,
                       abs(values[1]-row['angular_relative_difference']) < 1e-12,
                       row['passed'] == all(v < limit for v in values)))
    checks.append(summary['passed'] == all(r['passed'] for r in summary['runs']+summary['comparisons']))
    result = dict(passed=bool(all(checks)), checks=len(checks), maximum_endpoint_error=maximum, campaign_passed=summary['passed'],
                  scope='Independent endpoint reconstruction and loop interpolation; intermediate extrema use archived traces.')
    (ROOT/'local-transfer-audit.json').write_text(json.dumps(result, indent=2)+'\n'); print(json.dumps(result))
    assert result['passed']


if __name__ == '__main__':
    main()

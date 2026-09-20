"""SE-2: apply the declared SE-1S derivative to completed emitter states."""
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
NAMES = ['no-emission', 'no-excitation', 'emission-only', 'exchange-0',
         'exchange-50', 'exchange-200', 'exchange-negative', 'sign-mirror']


def analyze(name):
    path = ROOT/'emitter-v1'/name
    record = json.loads(path.with_suffix('.json').read_text())['summary']
    cfg = record['config']; n = cfg['n']; length = cfg['length']; h = length/n
    chi = cfg.get('chi', 0.); k0 = cfg.get('mix', .5)
    state = np.load(path.with_suffix('.npz'))['final']
    size = 6*n**3
    f = state[:size].reshape(6,n,n,n); pi = state[size:2*size].reshape(f.shape)
    u = .08*f[0]; A = f[1:4]; waves = f[4:6]; momenta = pi[4:6]
    current = np.zeros((3,n,n,n)); gradient_energy = np.zeros((n,n,n))
    for j in range(3):
        forward = (np.roll(waves,-1,j+1)-waves)/h
        backward = (waves-np.roll(waves,1,j+1))/h
        gradient_energy += np.sum(forward**2+backward**2,axis=0)/4
        current[j] = np.sum(momenta*(forward+backward)/2,axis=0)
    def coefficients(background):
        a = np.exp(4*background)
        beta = np.exp(2*background)*.5*.08*A/np.sqrt(1+.08**2*np.sum(A*A,axis=0))
        k = k0*np.tanh(background/.001)
        m2 = .75**2*np.exp(2*chi*background)
        return a,beta,k,m2
    def energy(background):
        a,beta,k,m2 = coefficients(background)
        return float(np.sum(.5*a*np.sum(momenta**2,axis=0)+gradient_energy
                            -np.sum(beta*current,axis=0)+.5*m2*(waves[1]-k*waves[0])**2)*h**3)
    a,beta,k,m2 = coefficients(u)
    residual = waves[1]-k*waves[0]
    ku = k0/.001*(1-np.tanh(u/.001)**2)
    terms = {'kinetic': 2*a*np.sum(momenta**2,axis=0),
             'shift': -2*np.sum(beta*current,axis=0),
             'mixing': m2*(chi*residual**2-ku*waves[0]*residual)}
    density = sum(terms.values()); integrated = float(density.sum()*h**3)
    total_energy = energy(u)
    numeric = (energy(u+1e-7)-energy(u-1e-7))/(2e-7)
    archived_energy = sum(record['final'][key] for key in ['radiation_propagation','companion_propagation','mixing_potential'])
    coordinates = np.arange(n)*h-length/2
    radius = np.sqrt(sum(x*x for x in np.meshgrid(coordinates,coordinates,coordinates,indexing='ij')))
    shells = []
    for lower,upper in zip([0,1,2,3,4],[1,2,3,4,np.inf]):
        selection = (radius >= lower)&(radius < upper)
        shells.append(dict(lower=lower,upper=None if np.isinf(upper) else upper,
                           source=float(density[selection].sum()*h**3),
                           absolute_source=float(abs(density[selection]).sum()*h**3)))
    null = name in ['no-emission','no-excitation']
    passed = abs(numeric-integrated) <= 1e-9+1e-5*abs(integrated) and abs(total_energy-archived_energy) < 1e-12
    if null: passed = passed and abs(total_energy)<1e-12 and abs(integrated)<1e-12
    return dict(name=name,energy=total_energy,source=integrated,
                response_per_energy=integrated/total_energy if total_energy>1e-20 else None,
                positive_source=float(np.maximum(density,0).sum()*h**3),
                negative_source=float(np.minimum(density,0).sum()*h**3),
                contributions={key:float(value.sum()*h**3) for key,value in terms.items()},
                numerical_source=numeric,derivative_error=abs(numeric-integrated),
                energy_error=abs(total_energy-archived_energy),shells=shells,
                minimum_U=float(u.min()),maximum_U=float(u.max()),
                internal_energy_change=record['final']['internal_rest']-record['initial']['internal_rest'],
                state_sha256=hashlib.sha256(path.with_suffix('.npz').read_bytes()).hexdigest(),
                passed=bool(passed))


def main():
    manifest = json.loads((ROOT/'emitter-v1/manifest.json').read_text())
    names = [c['name'] for c in manifest['configurations']]
    rows = [analyze(name) for name in names if (ROOT/'emitter-v1'/f'{name}.json').exists()]
    baseline = json.loads((ROOT/'source-response-v1/results.json').read_text())
    base = {r['name']:r for r in baseline['rows']}
    for row in rows:
        row['source_relative_to_emission_only'] = row['source']/base['emission-only']['source']
        row['source_relative_to_X_chi200'] = row['source']/base['exchange-200']['source']
    summary = dict(passed=all(row['passed'] for row in rows), complete=len(rows)==len(names),
                   completed_runs=len(rows), declared_runs=len(names),
                   source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   baseline_sha256=hashlib.sha256((ROOT/'source-response-v1/results.json').read_bytes()).hexdigest())
    # Derived checkpoint may be regenerated; immutable evolution states are never changed.
    (ROOT/'emitter-source-response.json').write_text(json.dumps(dict(summary=summary,rows=rows),indent=2)+'\n')
    lines = ['# SE-2 generated-source checkpoint', '',
             f"Completed states: {len(rows)}/{len(names)}. Diagnostic checks pass: {summary['passed']}.", '',
             '| Case | Wave energy | Scalar source D | D/energy | D / emission-only D |',
             '|---|---:|---:|---:|---:|']
    for row in rows:
        response='undefined' if row['response_per_energy'] is None else f"{row['response_per_energy']:.8g}"
        lines.append(f"| {row['name']} | {row['energy']:.8g} | {row['source']:.8g} | {response} | {row['source_relative_to_emission_only']:.8g} |")
    lines += ['', 'These are instantaneous wave-sector source derivatives at time 4, not total',
              'gravity, ray deflection or observational fits. Same finite initial internal',
              'excitation does not imply equal emitted energy. X/Y propagation plus mixing',
              'energy is counted together; kinetic and interaction contributions and shell',
              'integrals are preserved in the JSON. Candidate/source differences can change',
              'the subsequently evolved background as well as mode composition.', '',
              'No statement of convergence is made before the full spatial campaign and its',
              'independent audit finish. A passing partial diagnostic is not campaign success.',
              'Derived by emitter_source_response.py using the preregistered SE-1S derivative',
              'and SE-2 protocol; established Hamiltonian differentiation is credited.']
    (ROOT/'emitter-source-report.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary))
    for row in rows:
        print(row['name'],row['source'],row['response_per_energy'])
    assert summary['passed']


if __name__ == '__main__':
    main()

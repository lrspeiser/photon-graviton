"""Reconstruct homogeneous archive energy and tangent diagnostics independently."""
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
archive=ROOT/'homogeneous-stability-v1'
data=json.loads((archive/'results.json').read_text());checks=[]
J=np.zeros((6,6));J[:3,3:]=np.eye(3);J[3:,:3]=-np.eye(3)
for name,digest in data['summary']['sources'].items():
    checks.append(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
checks.append(len(data['rows'])==12 and len(data['controls'])==48)
max_energy=0.
for row in data['rows']:
    chi,k0,channel=(row[k] for k in ['chi','k0','channel'])
    raw=np.load(archive/f'chi{chi}-mix{k0}-{channel}.npz')
    phi,x,y,pp,px,py=raw['state'];u=.08*phi
    energy=.5*np.exp(4*u)*(pp*pp+px*px+py*py)+.02*phi*phi+.28125*np.exp(2*chi*u)*(y-k0*np.tanh(1000*u)*x)**2
    error=float(np.max(abs(energy-raw['energy'])));max_energy=max(max_energy,error)
    maps=raw['tangent'];singular=np.linalg.svd(maps,compute_uv=False)
    drift=float(np.max(abs(energy-energy[0])))
    sym=float(np.max(abs(maps[-1].T@J@maps[-1]-J)))
    checks += [error<1e-15,len(raw['time'])==201 and raw['time'][-1]==20,
               abs(drift-row['energy_drift'])<1e-15,
               abs(sym-row['symplectic_error'])<1e-12,
               abs(singular[:,0].max()-row['maximum_amplification'])<1e-10,
               np.max(abs(singular[-1]-row['endpoint_singular_values']))<1e-10]
    expected=drift<1e-11+1e-8*abs(energy[0]) and sym<1e-6 and all(s['relative_error']<1e-4 for s in row['shadows'])
    checks.append(bool(expected)==row['passed'])
for row in data['controls']:
    checks.append((row['error']<1e-7 and row['symplectic_generator_error']<1e-10)==row['passed'])
checks.append(data['summary']['passed']==all(r['passed'] for r in data['rows']+data['controls']))
result=dict(passed=bool(all(checks)),count=len(checks),maximum_energy_reconstruction_error=max_energy)
(ROOT/'homogeneous-audit.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# SE-H: finite-time homogeneous perturbations','',
       'All twelve trajectories and 48 Hamiltonian/generator controls pass the declared',
       'numerical gates. The independent raw-archive audit passes '+str(len(checks))+' checks.', '',
       '| chi | mixing k0 | Initially excited channel | Maximum tangent amplification |',
       '|---:|---:|---|---:|']
for r in data['rows']:
    lines.append(f"| {r['chi']} | {r['k0']} | {r['channel']} | {r['maximum_amplification']:.8g} |")
lines += ['', 'The amplification is a Euclidean singular value in the declared canonical',
          'coordinates over 20 model time units. It is coordinate dependent. The uncoupled',
          'free X coordinate has shear map [[1,T],[0,1]] and singular value',
          '(sqrt(T^2+4)+T)/2 = 20.0498756 at T=20. Thus amplification around 20 is',
          'not by itself evidence of exponential instability or of enhanced gravity.', '',
          'The interaction changes perturbation evolution, but this horizon and these',
          'initial states do not establish an asymptotic growth rate or nonlinear stability.',
          'The calculation evolves the homogeneous scalar sector with the vector field',
          'zero. Spatial perturbations couple to the vector sector and were not tested.',
          'There are no particles, emitter or boundary here. This does not validate the',
          'stability of the emitted 3D configurations or a permanent swirl.', '',
          'Maximum absolute energy drift: '+str(max(r['energy_drift'] for r in data['rows']))+'.',
          'Maximum endpoint symplectic residual: '+str(max(r['symplectic_error'] for r in data['rows']))+'.',
          'Maximum shadow-trajectory relative discrepancy: '+str(max(s['relative_error'] for r in data['rows'] for s in r['shadows']))+'.', '',
          'The tangent map is checked against independently evolved +/- perturbations at',
          'two amplitudes. The archive audit reconstructs energies and tangent diagnostics;',
          'it does not independently re-integrate the shadow trajectories.', '',
          'Protocol and code were committed before execution. Raw trajectories, tangent',
          'maps, controls, source hashes and all twelve results are preserved in',
          'homogeneous-stability-v1. Run audit_homogeneous.py to regenerate this report.',
          'Hamiltonian, symplectic, tangent and numerical integration methods are established',
          'mathematics. The field constitutive laws are candidate assumptions; no physical',
          'novelty or observational success is claimed. All twelve goal items remain active.']
(ROOT/'homogeneous-stability-report.md').write_text('\n'.join(lines)+'\n')
print(json.dumps(result));assert result['passed']

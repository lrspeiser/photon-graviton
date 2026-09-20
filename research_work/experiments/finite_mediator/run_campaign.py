"""Execute and preserve FM-1's 96 declared trajectories and controls."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import hashlib
import json
import platform
import subprocess
import time
from datetime import datetime, timezone
import numpy as np
import field_model as M


def plain(x):
    if isinstance(x, dict): return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)): return [plain(v) for v in x]
    if isinstance(x, np.generic): return x.item()
    return x


def save(path, obj):
    path.write_text(json.dumps(plain(obj), indent=2, allow_nan=False)+'\n', encoding='utf8', newline='\n')


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def integrate(spec, n=48, dt=.02):
    model = M.Model(spec, n)
    state = model.initial()
    e0 = model.energies(state)
    d0 = model.diagnostics(state)
    times, energies, paths, momenta, diagnostic = [], [], [], [], []
    boundary = 100.
    start = time.monotonic()
    for step in range(round(4/dt)+1):
        t = step*dt
        f, pi, q, p = model.unpack(state)
        boundary = min(boundary, float(8-np.max(abs(q))-model.a-(4-t)))
        # Every integrator endpoint is checked and retained, not just plotted snapshots.
        times.append(t)
        energies.append(model.energies(state))
        paths.append(q.copy())
        momenta.append(p.copy())
        diagnostic.append(model.diagnostics(state))
        if step < round(4/dt): state = M.rk4(model.rhs, state, dt)
    energies = np.array(energies)
    drift = np.max(abs(energies[:, 0]-e0[0]))/max(1, abs(e0[0]))
    md = max(np.linalg.norm(d['momentum']-d0['momentum']) for d in diagnostic)
    ad = max(abs(d['angular']-d0['angular']) for d in diagnostic)
    gates = dict(energy=bool(drift<=1e-4), boundary=bool(boundary>0), finite=bool(np.all(np.isfinite(state))))
    row = dict(spec=spec, n=n, dt=dt, seconds=time.monotonic()-start, times=times, energies=energies,
               diagnostics=diagnostic, energy_drift=drift, momentum_residual=md, angular_residual=ad,
               boundary_clearance=boundary, max_speed=max(d['max_speed'] for d in diagnostic),
               field_gain=energies[-1, 1], matter_loss=e0[2]-energies[-1, 2],
               final_position=paths[-1], final_momentum=momenta[-1], gates=gates, passed=all(gates.values()))
    arrays = dict(q=np.array(paths), p=np.array(momenta), final_state=state)
    return row, arrays


def main():
    out = M.HERE/'evidence-v1'
    out.mkdir(exist_ok=False)
    sources = sorted(list(M.HERE.glob('*.py'))+list(M.HERE.glob('*.md')))
    save(out/'manifest.json', dict(git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=M.ROOT, text=True).strip(),
         utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(), numpy=np.__version__,
         hashes={p.relative_to(M.ROOT).as_posix(): digest(p) for p in sources}))
    controls = M.controls()
    save(out/'controls.json', controls)
    if not all(r['passed'] for r in controls): raise AssertionError('Preserved Hamiltonian control failure')
    waves, wave_paths = M.wave_controls()
    save(out/'wave-controls.json', waves)
    np.savez_compressed(out/'wave-paths.npz', **wave_paths)
    print(f'{len(controls)} Hamiltonian controls pass; wave results: {waves}', flush=True)
    rows, archive = [], {}
    for source in ('rest', 'rotating', 'reverse'):
        for g in (0., .08):
            for eta in (-.08, 0., .08):
                for probe in ('none', 'photon', 'massive'):
                    sid = f'{source}-g{g:g}-eta{eta:g}-{probe}'
                    spec = dict(id=sid, source=source, g=g, eta=eta, probe=probe)
                    row, arrays = integrate(spec)
                    rows.append(row)
                    for k, v in arrays.items(): archive[sid+'__'+k] = v
                    save(out/'primary.json', rows)
                    print(f'primary {len(rows)}/54 {sid}: energy={row["energy_drift"]:.3g}', flush=True)
    np.savez_compressed(out/'primary-trajectories.npz', **archive)
    refined, comparisons, refined_archive = [], [], {}
    for coarse in rows:
        spec = coarse['spec']
        if spec['g'] != .08 or spec['eta'] not in (0., .08): continue
        for n, dt in ((48, .01), (64, .01)):
            fine, arrays = integrate(spec, n, dt)
            refined.append(fine)
            key = spec['id']+f'__n{n}'
            for k, v in arrays.items(): refined_archive[key+'__'+k] = v
            discrepancy = np.linalg.norm(archive[spec['id']+'__q']-arrays['q'][::2], axis=2).max()
            relative_energy = abs(coarse['field_gain']-fine['field_gain'])/max(fine['field_gain'], 1e-10)
            comparisons.append(dict(id=spec['id'], n=n, dt=dt, position_discrepancy=discrepancy,
                                    field_energy_relative_change=relative_energy,
                                    passed=bool(discrepancy<=.01 and relative_energy<=.03)))
            save(out/'refined.json', refined)
            save(out/'refinement-comparisons.json', comparisons)
            print(f'refinement {len(refined)}/36 {key}: position={discrepancy:.3g}, field={relative_energy:.3g}', flush=True)
    np.savez_compressed(out/'refined-trajectories.npz', **refined_archive)
    doubled, colors, color_archive = [], [], {}
    for old in rows:
        spec = old['spec']
        if spec['g'] != .08 or spec['eta'] not in (0., .08) or spec['probe'] != 'photon': continue
        row, arrays = integrate(dict(spec, photon_p=2e-4))
        doubled.append(row)
        diff = row['diagnostics'][-1]['probe_angle']-old['diagnostics'][-1]['probe_angle']
        diff = abs(np.arctan2(np.sin(diff), np.cos(diff)))
        colors.append(dict(id=spec['id'], angle_difference=diff, passed=bool(diff<=.001)))
        for k, v in arrays.items(): color_archive[spec['id']+'__'+k] = v
    save(out/'doubled-photons.json', doubled)
    save(out/'color-comparisons.json', colors)
    np.savez_compressed(out/'doubled-trajectories.npz', **color_archive)
    # Eta reversal: exact symmetry of full configuration, including off-axis probe.
    parity = []
    for row in rows:
        s = row['spec']
        if s['eta'] != .08: continue
        other = next(r for r in rows if all(r['spec'][k] == s[k] for k in ('source', 'g', 'probe')) and r['spec']['eta'] == -.08)
        err = np.max(abs(archive[s['id']+'__q']-archive[other['spec']['id']+'__q']))
        parity.append(dict(kind='eta_sign', id=s['id'], position_error=err))
    # Ring reflection changes labels theta -> -theta; probe-free cases only.
    order = [0, 5, 4, 3, 2, 1]
    for row in rows:
        s = row['spec']
        if s['source'] != 'rotating' or s['probe'] != 'none': continue
        other = next(r for r in rows if r['spec']['source']=='reverse' and r['spec']['probe']=='none' and r['spec']['g']==s['g'] and r['spec']['eta']==s['eta'])
        mirrored = archive[other['spec']['id']+'__q'][:, order].copy()
        mirrored[:, :, 1] *= -1
        parity.append(dict(kind='rotation_reflection', id=s['id'], position_error=np.max(abs(archive[s['id']+'__q']-mirrored))))
    save(out/'symmetries.json', parity)
    save(out/'evidence-hashes.json', {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print('FM-1 complete; all failures retained.', flush=True)


if __name__ == '__main__': main()

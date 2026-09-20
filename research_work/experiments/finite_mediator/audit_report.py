"""Independent endpoint energy, archive checks, and FM-1 figures."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import hashlib
import json
import subprocess
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import field_model as M

HERE = Path(__file__).resolve().parent
OUT = HERE/'evidence-v1'


def read(name):
    return json.loads((OUT/name).read_text())


def sha(data):
    return hashlib.sha256(data).hexdigest()


def independent_energy(state, spec, n):
    dx = 16/n
    size = 3*n*n
    count = 6 if spec['probe']=='none' else 7
    f = state[:size].reshape(3, n, n)
    pi = state[size:2*size].reshape(f.shape)
    q = state[2*size:2*size+2*count].reshape(count, 2)
    p = state[2*size+2*count:].reshape(count, 2)
    field = .5*dx**2*np.sum(pi**2+.2**2*f**2)
    for a in (1, 2): field += .5*np.sum((np.roll(f, -1, axis=a)-f)**2)
    coordinate = np.arange(n)*dx-8
    matter = 0.
    for i in range(count):
        wx = np.maximum(1-(((coordinate-q[i, 0]+8)%16-8)/.75)**2, 0)**3
        wy = np.maximum(1-(((coordinate-q[i, 1]+8)%16-8)/.75)**2, 0)**3
        weights = np.outer(wx, wy)
        weights /= weights.sum()
        sampled = np.sum(f*weights, axis=(1, 2))
        mass = 1. if i<6 else (0. if spec['probe']=='photon' else 1e-4)
        e = np.sqrt(mass**2+np.dot(p[i], p[i]))
        matter += e*np.exp(spec['g']*sampled[0]+spec['eta']*np.dot(sampled[1:], p[i])/e)
    return np.array([field+matter, field, matter])


def main():
    checks = []
    def check(name, error, tolerance):
        checks.append(dict(name=name, error=float(error), tolerance=tolerance, passed=bool(error<=tolerance)))
    for name, expected in read('evidence-hashes.json').items():
        check('archive_hash_'+name, int(sha((OUT/name).read_bytes()) != expected), 0)
    manifest = read('manifest.json')
    for name, expected in manifest['hashes'].items():
        check('working_source_'+name, int(sha((M.ROOT/name).read_bytes()) != expected), 0)
        committed = subprocess.check_output(['git', 'show', manifest['git_head']+':'+name], cwd=M.ROOT)
        check('pinned_source_'+name, int(sha(committed) != expected), 0)
    rows, refined, doubled = read('primary.json'), read('refined.json'), read('doubled-photons.json')
    primary = np.load(OUT/'primary-trajectories.npz')
    refinement = np.load(OUT/'refined-trajectories.npz')
    double_paths = np.load(OUT/'doubled-trajectories.npz')
    for name, count, expected in [('primary', len(rows), 54), ('refinement', len(refined), 36), ('photon', len(doubled), 6)]:
        check(name+'_count', abs(count-expected), 0)
    endpoint_errors, edges = [], []
    for group, archive in ((rows, primary), (refined, refinement), (doubled, double_paths)):
        for row in group:
            key = row['spec']['id']+(f'__n{row["n"]}' if group is refined else '')
            state = archive[key+'__final_state']
            independent = independent_energy(state, row['spec'], row['n'])
            err = np.max(abs(independent-np.array(row['energies'][-1])))
            endpoint_errors.append(float(err))
            check('independent_energy_'+key, err, 1e-12)
            check('trajectory_endpoint_'+key, np.max(abs(archive[key+'__q'][-1]-np.array(row['final_position']))), 1e-13)
            check('funding_account_'+key, abs(row['field_gain']-row['matter_loss']-(row['energies'][-1][0]-row['energies'][0][0])), 1e-12)
            f = state[:3*row['n']**2].reshape(3, row['n'], row['n'])
            edge = max(np.max(abs(f[:, 0, :])), np.max(abs(f[:, -1, :])), np.max(abs(f[:, :, 0])), np.max(abs(f[:, :, -1])))
            edges.append(float(edge))
    # Independent energy gradient at an evolved, source-funded state.
    example = next(r for r in rows if r['spec']['id']=='rotating-g0.08-eta0.08-photon')
    state = primary[example['spec']['id']+'__final_state'].copy()
    model = M.Model(example['spec'])
    r = model.rhs(state)
    expected = np.concatenate((-r[model.nf:2*model.nf]*model.dx**2, r[:model.nf]*model.dx**2,
                               -r[2*model.nf+2*model.k:], r[2*model.nf:2*model.nf+2*model.k]))
    rng = np.random.default_rng(8801)
    coordinates = list(rng.choice(2*model.nf, 24, replace=False))+list(range(2*model.nf, len(state)))
    for j in coordinates:
        step = 1e-8 if j>=len(state)-2 else 1e-6
        plus, minus = state.copy(), state.copy()
        plus[j] += step; minus[j] -= step
        fd = (independent_energy(plus, example['spec'], 48)[0]-independent_energy(minus, example['spec'], 48)[0])/(2*step)
        check(f'evolved_gradient_{j}', abs(fd-expected[j])/max(1, abs(expected[j])), 1e-5)
    all_rows = rows+refined+doubled
    comparisons = read('refinement-comparisons.json')
    colors = read('color-comparisons.json')
    symmetries = read('symmetries.json')
    selected = []
    for source in ('rest', 'rotating', 'reverse'):
        for eta in (0., .08):
            for probe in ('none', 'photon', 'massive'):
                row = next(r for r in rows if r['spec']['source']==source and r['spec']['g']==.08 and r['spec']['eta']==eta and r['spec']['probe']==probe)
                selected.append(dict(id=row['spec']['id'], field_energy=row['field_gain'], matter_loss=row['matter_loss'],
                                     angle=row['diagnostics'][-1]['probe_angle'], curl=row['diagnostics'][-1]['curl_max'],
                                     circulation=row['diagnostics'][-1]['circulation_inner'],
                                     final_probe_position=row['final_position'][-1] if probe!='none' else None))
    summary = dict(checks=checks, audit_passed=all(c['passed'] for c in checks),
        primary_count=len(rows), primary_numerical_passes=sum(r['passed'] for r in rows),
        refined_count=len(refined), refinement_numerical_passes=sum(r['passed'] for r in refined),
        refinement_comparison_passes=sum(c['passed'] for c in comparisons),
        color_passes=sum(c['passed'] for c in colors),
        maximum_energy_drift=max(r['energy_drift'] for r in all_rows),
        maximum_momentum_residual=max(r['momentum_residual'] for r in all_rows),
        maximum_angular_residual=max(r['angular_residual'] for r in all_rows),
        maximum_speed=max(r['max_speed'] for r in all_rows),
        minimum_boundary_clearance=min(r['boundary_clearance'] for r in all_rows),
        maximum_final_boundary_field=max(edges), independent_energy_max_error=max(endpoint_errors),
        maximum_position_refinement=max(c['position_discrepancy'] for c in comparisons),
        maximum_field_energy_refinement=max(c['field_energy_relative_change'] for c in comparisons),
        maximum_photon_energy_angle_difference=max(c['angle_difference'] for c in colors),
        maximum_symmetry_error=max(c['position_error'] for c in symmetries), selected=selected)
    (HERE/'audit.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n', encoding='utf8', newline='\n')
    if not summary['audit_passed']: raise AssertionError('Independent audit failure preserved')
    plot(rows, primary, comparisons)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('checks', 'selected')}, indent=2))
    print(json.dumps(selected, indent=2))


def plot(rows, paths, comparisons):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
    find = lambda source, eta, probe: next(r for r in rows if r['spec']['source']==source and r['spec']['g']==.08 and r['spec']['eta']==eta and r['spec']['probe']==probe)
    example = find('rotating', .08, 'none')
    energy = np.array(example['energies'])
    axes[0, 0].plot(example['times'], energy[:, 1], label='Field energy gained', lw=3)
    axes[0, 0].plot(example['times'], energy[0, 2]-energy[:, 2], '--', label='Matter energy lost', lw=2)
    axes[0, 0].set(title='A. Matter funds the initially empty field', xlabel='Time', ylabel='Model energy')
    axes[0, 0].legend()
    model = M.Model(example['spec'])
    f = model.unpack(paths[example['spec']['id']+'__final_state'])[0]
    curl = (np.roll(f[2], -1, 0)-np.roll(f[2], 1, 0)-np.roll(f[1], -1, 1)+np.roll(f[1], 1, 1))/(2*model.dx)
    limit = np.max(abs(curl))
    mesh = axes[0, 1].pcolormesh(model.x, model.y, curl, cmap='RdBu_r', vmin=-limit, vmax=limit, shading='auto', rasterized=True)
    axes[0, 1].quiver(model.x[::3, ::3], model.y[::3, ::3], f[1, ::3, ::3], f[2, ::3, ::3], color='#333333', alpha=.6)
    axes[0, 1].set(xlim=(-4,4), ylim=(-4,4), aspect='equal', title='B. Directional field and curl at t=4', xlabel='x', ylabel='y')
    fig.colorbar(mesh, ax=axes[0, 1], label='curl A')
    colors = {'rest':'#516575', 'rotating':'#147c8d', 'reverse':'#bb553d'}
    for index, probe in enumerate(('photon', 'massive')):
        ax = axes[0, 2] if index==0 else axes[1, 0]
        for source in ('rest', 'rotating', 'reverse'):
            for eta in (0., .08):
                row = find(source, eta, probe)
                angle = np.array([d['probe_angle'] for d in row['diagnostics']])*1e3
                ax.plot(row['times'], angle, color=colors[source], ls='--' if eta==0 else '-', label=f'{source}, eta={eta:g}')
        ax.set(title=('C. Photon direction' if index==0 else 'D. Massive-probe direction'), xlabel='Time', ylabel='Angle from initial path (mrad)')
        ax.axhline(0, color='gray', lw=.6)
        if index==0: ax.legend(fontsize=8)
    ax=axes[1, 1]
    waves=read('wave-controls.json')
    ax.loglog([w['n'] for w in waves], [w['relative_max_error'] for w in waves], 'o-', color='#147c8d')
    ax.axhline(.03, ls='--', color='#bb553d', label='Declared 3% limit')
    ax.set(title='E. Wave accuracy: coarse grid fails', xlabel='1D grid cells', ylabel='Relative maximum error')
    ax.legend()
    ax=axes[1, 2]
    for n, label, color in ((48, 'Time refinement only', '#516575'), (64, 'Space + time refinement', '#147c8d')):
        vals=[c['field_energy_relative_change'] for c in comparisons if c['n']==n]
        ax.semilogy(np.arange(1,19), np.maximum(vals,1e-16), 'o-', label=label, color=color)
    ax.axhline(.03, ls='--', color='#bb553d', label='Declared 3% limit')
    ax.set(title='F. Field-energy refinement', xlabel='Declared case index', ylabel='Relative change')
    ax.legend(fontsize=8)
    fig.suptitle('FM-1: matter-generated waves and derived light bending\nDimensionless 2D consistency experiment; no galaxy or cluster observations fitted', fontsize=16)
    fig.savefig(HERE/'results.png', dpi=170)
    fig.savefig(HERE/'results.svg')
    p=HERE/'results.svg'
    p.write_text('\n'.join(line.rstrip() for line in p.read_text(encoding='utf8').splitlines())+'\n', encoding='utf8', newline='\n')
    plt.close(fig)


if __name__=='__main__': main()

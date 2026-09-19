"""TF-1 stage 2: the longitudinal scan on CWC-1's 1D fixture (CWC-1's simulate reused unchanged)."""
import sys
import numpy as np
from common import CWC, save, read

sys.path.insert(0, str(CWC))
from spatial import simulate   # noqa: E402

ARCHIVE = CWC/'evidence/spatial1d-v1'
SCAN = [.005, .01, .02, .05, .1, .2]          # 0.05 and 0.2 at wavelength 2 exist in the archive; 0.5 is read from it
WAVELENGTHS = [1., 2., 4.]
N = 1024
KEEP = ['config', 'initial_sectors', 'final_sectors', 'sector_order', 'energy_relative_error', 'momentum_error_over_initial_em',
        'min_optical_speed', 'max_optical_speed', 'local_clock_speed_max_change', 'em_energy_lost',
        'material_internal_energy_lost', 'receiving_gain', 'max_final_phi', 'pulse_metrics']


def summary(result):
    return {k: result[k] for k in KEEP}


def archived(label):
    return summary(read(ARCHIVE/(label + '.json')))


def stretches(row, vacuum):
    pm, vm = row['pulse_metrics'], vacuum['pulse_metrics']
    row['spectral_stretch'] = vm['frequency_mean']/pm['frequency_mean']
    row['event_stretch'] = pm['clock_separation']/vm['clock_separation']
    row['timing_spectral_discrepancy'] = abs(row['event_stretch']/row['spectral_stretch'] - 1)
    return row


def power_law(g, y):
    lg, ly = np.log(np.asarray(g, float)), np.log(np.asarray(y, float))
    slope, intercept = np.polyfit(lg, ly, 1)
    resid = ly - (slope*lg + intercept)
    return dict(exponent=float(slope), intercept=float(intercept), log_scatter=float(np.sqrt(np.mean(resid**2))),
                max_log_residual=float(np.max(np.abs(resid))))


def run(out):
    cases = {}

    def execute(label, **kw):
        result, _ = simulate(n=N, **kw)
        row = summary(result)
        save(out/(label + '.json'), row)
        cases[label] = row
        print('1D', label, 'energy', row['energy_relative_error'], 'receiving', row['receiving_gain'], flush=True)
        return row

    vac = {}
    for wl in WAVELENGTHS:
        vac[wl] = execute(f'vacuum-{wl}', g=0, wavelength=wl)
    repro = execute('primary-0.2-2.0', g=.2, b=2, wavelength=2.)
    for g in SCAN:
        for wl in WAVELENGTHS:
            if g == .2 and wl == 2.:
                continue
            execute(f'coupling-{g}-{wl}', g=g, b=2, wavelength=wl)
    cases['coupling-0.2-2.0'] = repro
    cases['coupling-0.5-2.0'] = archived('coupling-0.5-clock-2')

    # reproduction against CWC-1's archive at the same grid
    arch = dict(vacuum=archived('vacuum-1024'), vacuum_1=archived('vacuum-color-1.0'), vacuum_4=archived('vacuum-color-4.0'),
                primary=archived('primary-1024'), coupling_0_05=archived('coupling-0.05-clock-2'))

    def rel(a, b):
        return abs(a - b)/max(abs(b), 1e-30)
    repro_rows = dict(
        receiving_gain=rel(repro['receiving_gain'], arch['primary']['receiving_gain']),
        frequency_mean=rel(repro['pulse_metrics']['frequency_mean'], arch['primary']['pulse_metrics']['frequency_mean']),
        clock_separation=rel(repro['pulse_metrics']['clock_separation'], arch['primary']['pulse_metrics']['clock_separation']),
        vacuum_frequency=rel(vac[2.]['pulse_metrics']['frequency_mean'], arch['vacuum']['pulse_metrics']['frequency_mean']),
        vacuum_1_frequency=rel(vac[1.]['pulse_metrics']['frequency_mean'], arch['vacuum_1']['pulse_metrics']['frequency_mean']),
        vacuum_4_frequency=rel(vac[4.]['pulse_metrics']['frequency_mean'], arch['vacuum_4']['pulse_metrics']['frequency_mean']),
        coupling_0_05_gain=rel(cases['coupling-0.05-2.0']['receiving_gain'], arch['coupling_0_05']['receiving_gain']))
    reproduction_worst = max(repro_rows.values())

    # per-coupling metrics
    table = []
    for g in SCAN + [.5]:
        rows = {}
        for wl in WAVELENGTHS:
            key = f'coupling-{g}-{wl}'
            if key in cases and cases[key]['pulse_metrics'] is not None:
                rows[wl] = stretches(dict(cases[key]), vac[wl])
        if 2. not in rows:
            continue
        main = rows[2.]
        if g == .5:                          # archived case: its vacuum is the archived 1024 vacuum
            main = stretches(dict(cases['coupling-0.5-2.0']), arch['vacuum'])
        spread = None
        if len(rows) == 3:
            s = [rows[w]['spectral_stretch'] for w in WAVELENGTHS]
            spread = float(np.ptp(s)/np.mean(s))
        table.append(dict(g=g, receiving_gain=main['receiving_gain'], em_energy_lost=main['em_energy_lost'],
                          spectral_stretch=main['spectral_stretch'], event_stretch=main['event_stretch'],
                          timing_discrepancy=main['timing_spectral_discrepancy'], colour_spread=spread,
                          fixed_ruler_speed_change=main['local_clock_speed_max_change'],
                          energy_relative_error=main['energy_relative_error'],
                          screens=dict(timing=main['timing_spectral_discrepancy'] < .01,
                                       colour=(spread < .01) if spread is not None else None,
                                       fixed_ruler=main['local_clock_speed_max_change'] < .001,
                                       conversion=main['receiving_gain'] > .001)))
    fit_rows = [t for t in table if t['g'] <= .1]
    gs = [t['g'] for t in fit_rows]
    fits = dict(receiving_gain=power_law(gs, [t['receiving_gain'] for t in fit_rows]),
                timing_discrepancy=power_law(gs, [t['timing_discrepancy'] for t in fit_rows]),
                colour_spread=power_law([t['g'] for t in fit_rows if t['colour_spread']], [t['colour_spread'] for t in fit_rows if t['colour_spread']]),
                fixed_ruler_speed_change=power_law(gs, [t['fixed_ruler_speed_change'] for t in fit_rows]))
    passing = [t['g'] for t in table if t['screens']['timing'] and t['screens']['colour'] and t['screens']['fixed_ruler']]
    largest_pass = max(passing) if passing else None
    gain_at_pass = next((t['receiving_gain'] for t in table if t['g'] == largest_pass), None)
    gates = dict(energy=all(t['energy_relative_error'] < 1e-3 for t in table),
                 reproduction=reproduction_worst < 1e-6,
                 power_laws=all(f['log_scatter'] < .1 for f in fits.values()))
    return dict(stage='scan', gates=gates, numerical_pass=all(gates.values()), reproduction=repro_rows,
                table=table, power_laws=fits, largest_coupling_passing_all_screens=largest_pass,
                receiving_gain_at_that_coupling=gain_at_pass,
                cosmological_rate_per_Mpc=70/299792.458,
                note='CWC-1 fixture at n = 1024, b = 2; the transverse channel is silent in one dimension, so every cost here is the longitudinal channel. No fixture-to-galaxy mapping.')

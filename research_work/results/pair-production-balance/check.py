"""Same-operator pair creation, inverse scattering and radiative loss balance."""
from pathlib import Path
import hashlib
import json
import os
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.sparse import bmat, csc_matrix, diags

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = Path(os.environ.get('PHOTON_GRAVITON_RESULTS', ROOT/'research_work/generated'))/HERE.name


def rates(energy, gap, polarizability):
    coefficient = polarizability**2/(6*np.pi)
    down = coefficient*energy*max(energy-gap, 0.)**3
    up = coefficient*energy*(energy+gap)**3
    annihilation = polarizability**2*gap**7/(1680*np.pi**3)
    return down, up, annihilation


def phase_space_checks(cfg):
    single_angle = 2*np.pi*quad(lambda u: 1-u*u, -1, 1, epsabs=1e-13)[0]
    pair_angle = 4*np.pi*2*np.pi*quad(lambda u: 1+u*u, -1, 1, epsabs=1e-13)[0]
    values = []
    energy = cfg['photon_energy']; lam = cfg['transition_polarizability']
    for ratio in cfg['gap_ratios']:
        gap = ratio*energy
        sigma_num = 2*np.pi*lam**2/(4*(2*np.pi)**3)*energy*(energy-gap)**3*single_angle
        integral = quad(lambda w: w**3*(gap-w)**3, 0, gap, epsabs=1e-60, epsrel=1e-12)[0]
        decay_num = .5*2*np.pi*lam**2/(4*(2*np.pi)**6)*pair_angle*integral
        sigma, _, decay = rates(energy, gap, lam)
        error = max(abs(sigma_num/sigma-1), abs(decay_num/decay-1))
        assert error < cfg['checks']['relative_phase_space_error']
        values.append({'gap': gap, 'scattering_cross_section': sigma,
                       'single_pair_decay_rate': decay, 'relative_phase_space_error': float(error)})
    return values


def evolve(initial, birth, inverse, annihilation, times, solver):
    count = len(initial)
    raw_birth = birth.copy(); birth = birth.copy(); birth[-1] = 0
    death = inverse+annihilation
    generator = diags([birth[:-1], -birth-death, death[1:]], [-1, 0, 1], format='csc')
    omitted = np.zeros(count); omitted[-1] = raw_birth[-1]
    counter_rows = csc_matrix(np.vstack([birth, inverse, annihilation, omitted]))
    augmented = bmat([[generator, csc_matrix((count, 4))],
                      [counter_rows, csc_matrix((4, 4))]], format='csc')
    state = np.r_[initial, np.zeros(4)]
    solution = solve_ivp(lambda t,y: augmented@y, (times[0], times[-1]), state,
                         t_eval=times, method='BDF', jac=augmented,
                         rtol=solver['rtol'], atol=solver['atol'])
    assert solution.success, solution.message
    p = solution.y[:count]
    n = np.arange(count)
    return p, n@p, solution.y[count:]


def population_case(ratio, gamma, solver, cfg):
    energy = cfg['photon_energy']; gap = ratio*energy
    sigma_down, sigma_up, decay = rates(energy, gap, cfg['transition_polarizability'])
    flux = decay/(gamma*sigma_down)
    creation = flux*sigma_down
    inverse_ratio = sigma_up/sigma_down
    r = 1/(inverse_ratio+gamma)
    stationary_mean = r/(1-r)
    n = np.arange(solver['maximum_pairs']+1, dtype=float)
    initial = np.zeros(len(n)); initial[0] = 1
    birth = (n+1)**2; inverse = inverse_ratio*n*n; annihilation = gamma*n*n
    times = np.array(cfg['illumination_times_in_creation_units'])
    p, mean, counters = evolve(initial, birth, inverse, annihilation, times, solver)
    norm_error = float(max(abs(p.sum(axis=0)-1)))
    assert norm_error < cfg['checks']['probability_error']
    assert p.min() > -cfg['checks']['probability_error']
    ledger = gap*(mean-counters[0]+counters[1]+counters[2])
    ledger_error = float(max(abs(ledger)))
    assert ledger_error < cfg['checks']['absolute_exchange_ledger_error']
    stationary_error = abs(float(mean[-1])-stationary_mean)
    assert stationary_error < cfg['checks']['final_stationary_mean_error'], stationary_error
    assert counters[3,-1] < cfg['checks']['integrated_omitted_births'], counters[3,-1]
    # Effective source area fixes a finite incident budget. It is large enough
    # for this optically thin single-store calculation over the numerical range.
    area = cfg['effective_source_area']
    assert max(sigma_down*(n+1)**2+sigma_up*n*n) < area
    input_energy = flux*area*energy*(times[-1]/creation)
    outgoing = input_energy-gap*counters[0,-1]+gap*counters[1,-1]+gap*counters[2,-1]
    assert outgoing >= 0
    # Isotropic independent incidents have zero mean momentum transfer. Bound
    # omitted leading finite-M kinetic heating using the maximum event impulse.
    recoil_bound = ((2*energy+gap)**2*(counters[0,-1]+counters[1,-1])
                    +gap**2*(counters[2,-1]+mean[-1]))/(2*cfg['nominal_seed_mass'])
    assert recoil_bound < cfg['checks']['recoil_energy_bound'], recoil_bound
    dark_times = np.array(cfg['dark_times_in_single_pair_decay_units'])
    dark_p, dark_mean, dark_counters = evolve(p[:,-1], np.zeros(len(n)), np.zeros(len(n)), n*n, dark_times, solver)
    dark_norm = float(max(abs(dark_p.sum(axis=0)-1)))
    dark_ledger_error = float(max(abs(gap*(dark_mean+dark_counters[2]-mean[-1]))))
    assert dark_norm < cfg['checks']['probability_error']
    assert dark_p.min() > -cfg['checks']['probability_error']
    assert dark_ledger_error < cfg['checks']['absolute_exchange_ledger_error']
    assert np.all(dark_mean <= mean[-1]*np.exp(-dark_times)+1e-8)
    stationary_p = (1-r)*r**n
    stationary_distance = float(np.sum(abs(p[:,-1]-stationary_p)))
    return {'gap_ratio':ratio, 'annihilation_to_creation_ratio':gamma, 'solver':solver,
            'photon_flux':float(flux), 'empty_mode_creation_rate':float(creation),
            'single_pair_dark_decay_rate':float(decay), 'inverse_to_creation_ratio':float(inverse_ratio),
            'stationary_pair_mean':float(stationary_mean), 'stationary_stored_energy':float(gap*stationary_mean),
            'stationary_mean_error':float(stationary_error), 'stationary_probability_l1_error':stationary_distance,
            'maximum_probability_normalization_error':max(norm_error,dark_norm),
            'maximum_exchange_ledger_error':max(ledger_error,dark_ledger_error),
            'omitted_birth_count_bound':float(counters[3,-1]),
            'expected_finite_mass_recoil_energy_bound':float(recoil_bound),
            'finite_incident_energy':float(input_energy), 'total_outgoing_energy':float(outgoing),
            'seed_mass_accounted_separately':cfg['nominal_seed_mass'],
            'illumination_times':times.tolist(), 'mean_pairs':mean.tolist(),
            'creation_counts':counters[0].tolist(), 'inverse_counts':counters[1].tolist(),
            'radiative_pair_loss_counts':counters[2].tolist(),
            'dark_times':dark_times.tolist(), 'dark_mean_pairs':dark_mean.tolist(),
            'dark_radiative_pair_loss_counts':dark_counters[2].tolist()}


def main():
    cfg = json.loads((HERE/'protocol.json').read_text(encoding='utf-8'))
    phase = phase_space_checks(cfg)
    cases = []
    for ratio in cfg['gap_ratios']:
        for gamma in cfg['annihilation_to_empty_creation_ratios']:
            runs = [population_case(ratio,gamma,solver,cfg) for solver in cfg['runs']]
            a,b = runs
            mean_change = max(max(abs(np.array(a[key])-b[key])) for key in ['mean_pairs','dark_mean_pairs'])
            count_change = max(max(abs(np.array(a[key])-b[key])/(1+np.array(b[key])))
                               for key in ['creation_counts','inverse_counts','radiative_pair_loss_counts'])
            assert mean_change < cfg['checks']['refinement_mean_error'], mean_change
            assert count_change < cfg['checks']['refinement_relative_count_error'], count_change
            cases.append({'runs':runs,'maximum_mean_refinement_change':float(mean_change),
                          'maximum_relative_count_refinement_change':float(count_change)})
            print('gap/E',ratio,'C/A',gamma,'mean',b['stationary_pair_mean'],
                  'stored/E',b['stationary_stored_energy']/cfg['photon_energy'],'refinement',mean_change)
    single_solver = cfg['runs'][-1]
    p,m,c = evolve(np.array([0.,1.]), np.zeros(2), np.zeros(2), np.array([0.,1.]),
                   np.array(cfg['dark_times_in_single_pair_decay_units']), single_solver)
    single_error = float(max(abs(m-np.exp(-np.array(cfg['dark_times_in_single_pair_decay_units'])))))
    assert single_error < cfg['checks']['single_pair_decay_error']
    soft = [{'gap_ratio':q,'limiting_stored_energy_over_photon_energy':
             float(q/np.expm1(3*(np.log1p(q)-np.log1p(-q))))} for q in [.1,.01,.0001,1e-8]]
    result = {'scope':cfg['scope'],'checks_pass':True,'phase_space_checks':phase,'population_cases':cases,
              'single_pair_exponential_decay_error':single_error,'negligible_annihilation_soft_gap_limit':soft,
              'astronomical_storage_lifetime_determined':False,'coherent_propagation_time_field_derived':False,
              'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'check.py',HERE/'protocol.json']}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'pair-production-balance-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__ == '__main__':
    main()

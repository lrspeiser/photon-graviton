"""PF-1 driver: T1 driven simulation, T2 coupled energy exchange, T3 analytic observables.

    python pf1.py [--output-dir DIR] [--canonical]

Regenerates pf1-results.json into a fresh directory and compares it with the archived copy;
--canonical overwrites the archive. Exits nonzero if a tolerance declared in protocol.md fails
or the regenerated numbers differ from the archive. T4 is brightness.py.
"""
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np
from scipy.integrate import quad, solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent/'companion-extensions'))
import evidence_io  # noqa: E402
import wave as W  # noqa: E402

C_KMS = 299792.458
ALPHA = 0.0002488993          # per Mpc: archived illustrative conversion constant (paper equation 1)


def t1(ndot=.005, nu=1., sigma=3., emissions=(-100., -76.), points=32):
    """Driven pulses in physical time; linear index with fixed separation D giving 1+z = 2.
    Pulses sit 8 sigma apart and are measured in adjacent +-4 sigma windows."""
    out, ok = {}, True
    D = np.log(2)/ndot
    const = dict(n=lambda t: 1. + 0*np.asarray(t, float), tau=lambda t: np.asarray(t, float),
                 t_of_tau=lambda s: np.asarray(s, float))
    for label, field in [('linear_index', W.linear(ndot)), ('constant_index', const)]:
        n, tau = field['n'], field['tau']
        arrivals = [float(W.arrival(field, te, D)) for te in emissions]
        want = float(n(arrivals[0])/n(emissions[0]))       # predicted 1+z = n(t_o)/n(t_e)
        t_start = emissions[0] - 6*sigma
        t_end = arrivals[-1] + 5*sigma*want
        # Long enough that nothing emitted after t_start reaches the outflow edge before t_end, so
        # no energy leaves the grid and no boundary reflection can reach the observer.
        length = float(tau(t_end) - tau(t_start)) + 5
        sim = W.simulate(n, t_start, t_end, D, length, list(emissions), nu, sigma, points_per_wavelength=points)
        obs = W.pulse_measurements(sim['t'], sim['signal'], [(a - 4*sigma*want, a + 4*sigma*want) for a in arrivals])
        te_grid = np.linspace(t_start, emissions[-1] + 6*sigma, 40001)
        emi = W.pulse_measurements(te_grid, W.emitter_signal(te_grid, emissions, nu, sigma),
                                   [(e - 4*sigma, e + 4*sigma) for e in emissions])
        carrier = emi[0]['carrier']/obs[0]['carrier']
        carrier2 = emi[1]['carrier']/obs[1]['carrier']
        width = obs[0]['width']/emi[0]['width']
        spacing = (obs[1]['centroid'] - obs[0]['centroid'])/(emi[1]['centroid'] - emi[0]['centroid'])
        tol = .005 if label == 'linear_index' else .002
        passed = all(abs(v/want - 1) < tol for v in (carrier, carrier2, width, spacing))
        E = sim['energy']
        m = E[:, 0] > emissions[-1] + 5*sigma                 # emission finished
        En = E[m, 1]*n(E[m, 0])
        drift = float(np.max(np.abs(En/En[0] - 1)))
        passed &= drift < .005
        ok &= passed
        out[label] = dict(expected_one_plus_z=want, carrier_ratio=carrier, carrier_ratio_second_pulse=carrier2,
                          width_ratio=width, spacing_ratio=spacing, energy_times_n_max_drift=drift,
                          wave_energy_first_last=[float(E[m, 1][0]), float(E[m, 1][-1])],
                          index_first_last=[float(n(E[m, 0][0])), float(n(E[m, 0][-1]))],
                          arrivals=arrivals, observed=obs, emitted=emi, dx=sim['dx'], dt=sim['dt'], length=length,
                          passed=bool(passed))
    out['passed'] = bool(ok)
    return out


def t2():
    """Radiation modes coupled to a dynamical index: free field and a restoring potential."""
    k = np.array([1., 2., 3., 5.])
    A0, P0 = np.array([.3, .2, .1, .05]), np.zeros(4)
    E_tau = float(np.sum(P0*P0 + k*k*A0*A0)/2)
    cases = dict(free=(40., lambda n: 0*n, lambda n: 0*n),
                 restoring=(40., lambda n: .01*(n - 1.5)**2, lambda n: .02*(n - 1.5)))
    out, ok = {}, True
    for name, (M, V, dV) in cases.items():
        r = W.coupled(k, A0, P0, 1., .002, M, V, dV, 300., points=20001)
        H = r['total']
        N = r['mode_energy']*r['n']/k[:, None]
        dEdt = np.gradient(r['wave_energy'], r['t'])
        pred = -(r['ndot']/r['n'])*r['wave_energy']
        inner = slice(5, -5)
        red = solve_ivp(lambda t, y: [y[1], (-dV(y[0]) + E_tau/y[0]**2)/M], (0, 300.), [1., .002], method='DOP853',
                        rtol=1e-12, atol=1e-14, dense_output=True)
        n_red = red.sol(r['t'])[0]
        gain = r['field_energy'] - r['field_energy'][0]
        loss = r['wave_energy'][0] - r['wave_energy']
        res = dict(energy_drift=float(np.max(np.abs(H/H[0] - 1))),
                   photon_number_drift=float(np.max(np.abs(N/N[:, :1] - 1))),
                   exchange_rate_residual=float(np.max(np.abs(dEdt[inner] - pred[inner]))/np.max(np.abs(pred[inner]))),
                   reduced_mechanics_max_dn=float(np.max(np.abs(n_red - r['n']))),
                   field_gain_minus_wave_loss=float(np.max(np.abs(gain - loss))/r['wave_energy'][0]),
                   n_range=[float(r['n'].min()), float(r['n'].max())],
                   wave_energy_range=[float(r['wave_energy'].min()), float(r['wave_energy'].max())])
        passed = (res['energy_drift'] < 1e-9 and res['photon_number_drift'] < 1e-9 and res['reduced_mechanics_max_dn'] < 1e-8
                  and res['exchange_rate_residual'] < 1e-3 and res['field_gain_minus_wave_loss'] < 1e-9)
        res['passed'] = bool(passed)
        ok &= passed
        out[name] = res
    out['passed'] = bool(ok)
    return out


def t3(omega_m=.3):
    """Distance-modulus differences against a flat FLRW comparator with the same low-z slope, H0 = c alpha."""
    H0 = C_KMS*ALPHA
    rows = []
    for z in [.1, .5, 1., 1.5]:
        dl_pf = float(W.luminosity_distance(z, ALPHA))
        dc = C_KMS/H0*quad(lambda s: 1/np.sqrt(omega_m*(1 + s)**3 + 1 - omega_m), 0, z)[0]
        dl_fl = (1 + z)*dc
        rows.append(dict(z=z, DL_PF1_Mpc=dl_pf, DL_FLRW_Mpc=dl_fl, delta_mu_PF1_minus_FLRW=5*np.log10(dl_pf/dl_fl),
                         time_dilation_PF1=1 + z, energy_loss_factor=1 + z))
    return dict(H0_equivalent_kms_Mpc=H0, alpha_per_Mpc=ALPHA, omega_m=omega_m, rows=rows,
                redshift_law='1+z = exp(alpha D) for a linear index; b = 1 exactly', passed=True)


def main():
    args = evidence_io.parse(__doc__)
    t0 = time.time()
    result = dict(T1=t1(), T2=t2(), T3=t3())
    result.update(scope='PF-1 internal-consistency tests of the proposed propagation-field wave law; not an '
                        'observational validation.',
                  protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),
                  source_sha256={p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in ('wave.py', 'pf1.py')},
                  all_passed=bool(all(result[k]['passed'] for k in ('T1', 'T2', 'T3'))), runtime_seconds=time.time() - t0)
    text = json.dumps(result, indent=1) + '\n'
    summary = {k: {kk: vv for kk, vv in v.items() if kk not in ('observed', 'emitted')} if isinstance(v, dict) else v
               for k, v in result['T1'].items()}
    print(json.dumps(dict(T1=summary, T2=result['T2'], T3=result['T3']['rows'], all_passed=result['all_passed']),
                     indent=1, default=float))
    status = evidence_io.finish(args, 'propagation-field-pf1', text, HERE/'pf1-results.json', ignore={'/runtime_seconds'})
    return status if result['all_passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())

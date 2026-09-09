"""Conditional fixed-volume energy ledgers; no microscopic/gravity claim.

All energies and times must use a consistent unit system. Companion and deposit
losses enter an explicit receiving sector. Stellar injection spends finite fuel.
"""
from dataclasses import dataclass, asdict
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Rates:
    photon_loss: float = 1.0
    capture: float = 0.2
    companion_loss: float = 0.0
    deposit_loss: float = 0.0

    def validate(self):
        if any(not np.isfinite(x) or x < 0 for x in asdict(self).values()):
            raise ValueError('Rates must be finite and nonnegative.')


def matrix(rates, injection):
    """State: photons, companions, deposits, receiving sector, fuel, constant."""
    rates.validate()
    h, g, lc, ld = (rates.photon_loss, rates.capture,
                    rates.companion_loss, rates.deposit_loss)
    a = np.zeros((6, 6))
    a[0, 0], a[0, 5] = -h, injection
    a[1, 0], a[1, 1] = h, -(g + lc)
    a[2, 1], a[2, 2] = g, -ld
    a[3, 1], a[3, 2] = lc, ld
    a[4, 5] = -injection
    return a


def solve(times, rates, initial_photons=1.0, fuel=0.0,
          injection=0.0, injection_duration=0.0):
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or not len(times) or not np.isfinite(times).all()
            or np.any(times < 0) or np.any(np.diff(times) < 0)):
        raise ValueError('Times must be finite, nonnegative and sorted.')
    if any(not np.isfinite(x) or x < 0 for x in
           [initial_photons, fuel, injection, injection_duration]):
        raise ValueError('Initial energies and source settings must be finite and nonnegative.')
    end = min(injection_duration, fuel / injection) if injection > 0 else 0.0
    on, off = matrix(rates, injection), matrix(rates, 0.0)
    y0 = np.array([initial_photons, 0., 0., 0., fuel, 1.])
    at_end = expm(on * end) @ y0
    values = np.array([(expm(on*t) @ y0 if t <= end
                        else expm(off*(t-end)) @ at_end)[:5] for t in times])
    return values


def analytic_impulse(times, h, g):
    """Independent no-loss closed form, with zero-rate/equal-rate limits."""
    t = np.asarray(times, dtype=float)
    photon = np.exp(-h*t)
    if h == 0:
        companion = np.zeros_like(t)
    elif g == h:
        companion = h*t*np.exp(-h*t)
    else:
        companion = h/(g-h)*(np.exp(-h*t)-np.exp(-g*t))
    return np.column_stack([photon, companion, 1-photon-companion])


def verify():
    records = []
    t = np.linspace(0, 12, 121)
    for h, g in [(1., .2), (1., 1.), (1., 2.), (1., 0.), (0., 1.)]:
        got = solve(t, Rates(h, g))
        expected = analytic_impulse(t, h, g)
        error = float(np.max(np.abs(got[:, :3]-expected)))
        assert error < 1e-11
        records.append({'check':'analytic_impulse', 'h':h, 'capture':g, 'error':error})
    for rates in [Rates(), Rates(companion_loss=1.), Rates(deposit_loss=.1),
                  Rates(companion_loss=.4, deposit_loss=.3)]:
        got = solve(t, rates, initial_photons=0., fuel=2., injection=.5,
                    injection_duration=100.)
        # Independent adaptive ODE, explicitly split at fuel exhaustion t=4.
        def rhs(inject):
            h,g,lc,ld = rates.photon_loss,rates.capture,rates.companion_loss,rates.deposit_loss
            return lambda _, y: [inject-h*y[0], h*y[0]-(g+lc)*y[1],
                                  g*y[1]-ld*y[2], lc*y[1]+ld*y[2], -inject]
        first = solve_ivp(rhs(.5), [0, 4], [0,0,0,0,2], dense_output=True,
                          rtol=1e-11, atol=1e-13)
        second = solve_ivp(rhs(0.), [4, 12], first.y[:,-1], dense_output=True,
                           rtol=1e-11, atol=1e-13)
        expected = np.array([first.sol(v) if v<=4 else second.sol(v) for v in t])
        err=float(np.max(np.abs(got-expected)))
        residual=float(np.max(np.abs(got.sum(axis=1)-2.)))
        assert err<1e-8 and residual<1e-11 and got.min()>-1e-12
        records.append({'check':'finite_fuel_ode_and_conservation', 'rates':asdict(rates),
                        'ode_error':err,'total_energy_residual':residual})
    shared=solve(t,Rates(capture=0.,companion_loss=1.))
    assert np.max(np.abs(shared[:,1]-t*np.exp(-t)))<1e-11
    records.append({'check':'shared_loss_companion_solution', 'passed':True})
    return records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    verification=verify()
    t=np.array([0.,1.,4.,10.,50.])
    scenarios={}
    for label,rates in [('no_loss_permanent',Rates()),
                        ('shared_secondary_loss',Rates(companion_loss=1.)),
                        ('finite_residence',Rates(deposit_loss=.1))]:
        scenarios[label]={'rates':asdict(rates),'times':t.tolist(),
            'energies':solve(t,rates,initial_photons=0.,fuel=2.,injection=.5,
                             injection_duration=4.).tolist()}
    result={'scope':'dimensionless bookkeeping verification; not a physical interaction or halo fit',
            'columns':['photons','companions','deposits','explicit_receiving_sector','stellar_fuel'],
            'checks':verification,'finite_fuel_scenarios':scenarios}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(f'Passed {len(verification)} checks. Results: {args.output}')


if __name__ == '__main__':
    main()

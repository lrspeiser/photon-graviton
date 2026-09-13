"""Retrospective published-summary test of the stationary conversion branch."""
from pathlib import Path
import json,math

OUT=Path(__file__).resolve().parent
b=1.003
rows=[]
for z in [.1,.2,.5,1.,1.2]:
    stretch=(1+z)**b
    survival=1/(1+z)
    rows.append(dict(z=z,stationary_conversion_duration_ratio=1.,published_fitted_duration_ratio=stretch,
        fraction_of_fitted_duration_missing=1-1/stretch,
        stationary_flux_over_same_distance_unshifted_flux=survival,
        flux_if_measured_timing_also_supplied=survival/stretch,
        required_arrival_delay_derivative=stretch-1,
        required_intrinsic_width_ratio_if_all_stretch_assigned_to_source=stretch))
    # Common delays cannot change event spacing; differential delay is needed.
    dt=12.;T=1e10
    assert (T+dt)-T==dt
    assert math.isclose(dt+dt*(stretch-1),dt*stretch)
result=dict(candidate='S0: static Euclidean paths, stationary per-photon energy conversion, unchanged c and local clocks, photon count conserved, no source duration evolution.',
    provenance='Known exponential-loss solution and arrival-time bookkeeping; companion interpretation hypothetical. Published timing exponent is an external inferred summary, not our raw-data fit.',
    published_source='https://arxiv.org/html/2406.05050v2',published_b=b,published_statistical_uncertainty=.005,published_systematic_estimate=.010,
    prediction_b=0.,exponent_difference=b,
    verdict='Incompatible with the published duration relation conditional on its source-comparison and reduction assumptions. No calibrated significance from these summary numbers. Not an exclusion of all nonexpanding models.',
    formulas=['1+z=exp(integral alpha ds)',
    't_arrival=t_emission+T(D), so dt_arrival/dt_emission=1 for stationary T',
    'F=L/[4*pi*D^2*(1+z)*S] for photon number conservation and achromatic energy shift',
    'A revised delay law must give dT/dt_emission=S-1; setting S=(1+z)^b by hand is a fit prescription, not a causal derivation'],
    scope='Retrospective benchmark. Rows sample a published fitted relation; they are not independent supernova measurements or held-out tests.',rows=rows)
(OUT/'timing-verdict.json').write_bytes((json.dumps(result,indent=2)+'\n').encode())
print(json.dumps(rows,indent=2))

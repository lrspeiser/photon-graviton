"""Exact energy-population recursion for a fixed-dump candidate until boundary failure."""
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
step=.05
n=np.arange(96,105)
ein=2*n*step
a=np.exp(-((ein-10)/.22)**2/4)
weights=a*a/np.sum(a*a)
loss=n*step
dump=float(weights@loss)
increments=n-100
variance=float(weights@(loss-dump)**2)
# Same initial receiver: U=50, mean excess 5, sigma=.5, initial support 0..10.
initial=np.exp(-((step*np.arange(201)-5)/.5)**2/2)
initial/=initial.sum()
p=np.zeros(2001)
p[:len(initial)]=initial
lower=upper=0.
rows=[]
for use in range(1,10001):
    nxt=np.zeros_like(p)
    for shift,w in zip(increments,weights):
        if shift<0:
            lower+=float(w*p[:-shift].sum())
            nxt[:shift]+=w*p[-shift:]
        elif shift>0:
            upper+=float(w*p[-shift:].sum())
            nxt[shift:]+=w*p[:-shift]
        else:
            nxt+=w*p
    p=nxt
    if use in [1,100,1000,10000]:
        rows.append(dict(uses=use,permitted_branch_weight=float(p.sum()),
                         lower_boundary_failure_weight=lower,upper_grid_escape_weight=upper,
                         unbounded_reference_mean_excess=5.,
                         unbounded_reference_excess_variance=float(initial@((step*np.arange(201)-5)**2)+use*variance)))
assert abs(p.sum()+lower+upper-1)<1e-10
result=dict(scope='Conditional fixed export per use; failure tracks leaving its allowed energy domain',
    units='Arbitrary energy, c=hbar=1',receiver_ground_energy=50.,
    initial_mean_excess=5.,initial_sigma=.5,export_energy_per_use=dump,
    incoming_mean_photon_energy=float(weights@ein),photon_loss_variance=variance,
    increment_mean=float(weights@increments)*step,
    source_population_is_repeated_independent_identical_photon_packets=True,
    boundary_completion_is_absorbing_failure_flag_not_a_physical_repair=True,
    rows=rows)
(HERE/'partial-export-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result,indent=2))

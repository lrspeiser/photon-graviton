"""Declared grid follow-up; previous failed evidence remains unchanged."""
import numpy as np
from common import HERE,hashes,read,save
from spatial import simulate
from one_dimensional import summary

def run(out):
    previous=HERE/'evidence/spatial1d-v1/results.json';prior=read(previous)
    old={r['label']:r for r in prior['cases']}
    rows=[]
    for label,kwargs in [('primary',{}),('vacuum',dict(g=0)),
                         ('color-1',dict(wavelength=1.)),('vacuum-color-1',dict(wavelength=1.,g=0)),
                         ('color-4',dict(wavelength=4.)),('vacuum-color-4',dict(wavelength=4.,g=0))]:
        result,_=simulate(n=4096,**kwargs);save(out/(label+'.json'),result)
        rows.append(summary(result,label))
        print('REFINE',label,result['energy_relative_error'],flush=True)
    cases={r['label']:r for r in rows}
    for label,base in [('primary','vacuum'),('color-1','vacuum-color-1'),('color-4','vacuum-color-4')]:
        row=cases[label];metrics=row['pulse_metrics'];baseline=cases[base]['pulse_metrics']
        row['spectral_stretch']=baseline['frequency_mean']/metrics['frequency_mean']
        row['event_stretch']=metrics['clock_separation']/baseline['clock_separation']
        row['timing_spectral_discrepancy']=abs(row['event_stretch']/row['spectral_stretch']-1)
    fine=cases['primary'];coarse=old['primary-2048']
    arrival=max(abs(a['coordinate_centroid']-b['coordinate_centroid']) for a,b in
                zip(fine['pulse_metrics']['pulses'],coarse['pulse_metrics']['pulses']))
    receiving=abs(fine['receiving_gain']-coarse['receiving_gain'])
    colors=[cases[k]['spectral_stretch'] for k in ['color-1','primary','color-4']]
    gates=dict(energy=fine['energy_relative_error']<1e-3,arrival_refinement=arrival<.02,
               receiving_refinement=receiving<.02,momentum=fine['momentum_error_over_initial_em']<.02)
    return dict(stage='refinement',input_sha256=hashes([previous]),cases=rows,gates=gates,
                numerical_pass=all(gates.values()),mechanism_pass=False,
                arrival_change=arrival,receiving_change=receiving,
                color_stretch_fractional_spread=np.ptp(colors)/np.mean(colors),
                note='Fixed original tolerances and windows. Colors refined, not retuned. Other physical failures and original failure remain.')

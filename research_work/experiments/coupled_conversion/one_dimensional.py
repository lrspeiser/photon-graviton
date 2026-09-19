"""Stage2: finite spatial conversion, clocks, colors and reciprocal recoil."""
import numpy as np
from common import save
from spatial import simulate

def summary(result,label):
    return dict(label=label,**{k:v for k,v in result.items()
                if k not in ['snapshots','receiver','time','energy_history','momentum_history','central_receiving_energy']})

def run(out):
    rows=[];cases={}
    def execute(label,**kwargs):
        result,model=simulate(**kwargs)
        save(out/(label+'.json'),result)
        row=summary(result,label);rows.append(row);cases[label]=row
        print('1D',label,'energy',row['energy_relative_error'],'receiving',row['receiving_gain'],flush=True)
        return row
    for n in [512,1024,2048]:
        execute('vacuum-'+str(n),n=n,g=0)
        execute('primary-'+str(n),n=n)
    for wavelength in [1.,4.]:
        execute('vacuum-color-'+str(wavelength),wavelength=wavelength,g=0)
        execute('color-'+str(wavelength),wavelength=wavelength)
    for g in [.05,.2,.5]:
        for b in [0,1,2]:
            if g==.2 and b==2:continue
            execute(f'coupling-{g}-clock-{b}',g=g,b=b)
    execute('time-refined',n=1024,step_factor=.04)
    for mass in [0.,1.]:execute('mass-'+str(mass),m=mass)
    execute('slow-receiving-wave',v=.5)
    for energy in [.25,4.]:execute('intensity-'+str(energy),energy=energy)
    execute('source-off',energy=0)
    execute('source-off-protected',energy=0,b=0)
    execute('zero-clock-energy',clock_energy=0)
    execute('reversed',reverse=True)
    execute('free-wave',energy=0,g=0,clock_energy=0,free=True)
    for row in rows:
        metrics=row['pulse_metrics'];config=row['config']
        if metrics is None:continue
        basekey='vacuum-'+str(config['n'])
        if config['wavelength']!=2.:basekey='vacuum-color-'+str(config['wavelength'])
        base=cases[basekey]['pulse_metrics']
        row['spectral_stretch']=base['frequency_mean']/metrics['frequency_mean']
        row['event_stretch']=metrics['clock_separation']/base['clock_separation']
        row['timing_spectral_discrepancy']=abs(row['event_stretch']/row['spectral_stretch']-1)
        row['width_ratios']=[a['pulse_width']/b['pulse_width'] for a,b in zip(metrics['pulses'],base['pulses'])]
    coarse=cases['primary-1024'];fine=cases['primary-2048'];half=cases['time-refined']
    receiving_change=abs(fine['receiving_gain']-coarse['receiving_gain'])
    arrival_change=max(abs(fine['pulse_metrics']['pulses'][j]['coordinate_centroid']-
                           coarse['pulse_metrics']['pulses'][j]['coordinate_centroid']) for j in [0,1])
    colors=[cases['color-1.0'],coarse,cases['color-4.0']]
    color_spread=np.ptp([r['spectral_stretch'] for r in colors])/np.mean([r['spectral_stretch'] for r in colors])
    gates=dict(energy=fine['energy_relative_error']<1e-3,
               time_refinement=half['energy_relative_error']<coarse['energy_relative_error'],
               spatial_receiving=receiving_change<.02,spatial_arrival=arrival_change<.02,
               momentum=fine['momentum_error_over_initial_em']<.02,
               no_coupling=cases['vacuum-2048']['max_final_phi']==0,
               no_source_protected=cases['source-off-protected']['max_final_phi']==0,
               free_transport=cases['free-wave']['energy_relative_error']<1e-3)
    physics=dict(net_conversion=fine['receiving_gain']>.001,
                 positive_redshift=fine['spectral_stretch']>1.001,
                 timing=abs(fine['event_stretch']/fine['spectral_stretch']-1)<.01,
                 achromatic=color_spread<.01,
                 fixed_ruler_speed=fine['local_clock_speed_max_change']<.001,
                 photon_exclusive=cases['source-off']['max_final_phi']<1e-10,
                 universal_matter_force=False)
    return dict(stage='one_dimensional',gates=gates,numerical_pass=all(gates.values()),
                mechanism_gates=physics,mechanism_pass=all(physics.values()),cases=rows,
                receiving_refinement_change=receiving_change,arrival_refinement_change=arrival_change,
                color_stretch_fractional_spread=color_spread,
                note='Initial finite pulse train, effective clocks, no halo or expanding geometry. Individual JSON files retain raw histories.')

"""Source-width comparison using the previously derived closed Hamiltonian."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import os
import numpy as np
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PARENT=HERE.parent/'generated-wave-access'
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'source-timescale'
spec=importlib.util.spec_from_file_location('generated_wave_solver',PARENT/'check.py')
solver=importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)


def main():
    protocol=json.loads((HERE/'protocol.json').read_text())
    baseline=json.loads((PARENT/'protocol.json').read_text())
    runs=[];convergence=[];signals={};weak_limits=[]
    for a in protocol['conversion_half_widths']:
        cfg=copy.deepcopy(baseline)
        cfg['parameters']['conversion_half_width']=a
        for m in protocol['restoring_frequencies']:
            for cells in protocol['cells']:
                row,signal=solver.simulate(cells,protocol['smoothing_width'],m,cfg)
                row['conversion_half_width']=a
                row['photon_fraction_lost']=row['photon_energy_lost']/row['photon_initial_energy']
                row['transmitted_over_initial_photon_energy']=row['transmitted_energy']/row['photon_initial_energy']
                if m==0: signals[(a,cells)]=signal
                parameters=cfg['parameters']
                dt=parameters['end_time']/(parameters['samples']-1)
                slab_half_width=(parameters['slab_edges'][1]-parameters['slab_edges'][0])/2
                row['source_spectrum']=solver.spectral_prediction(signals[(a,cells)],dt,m,parameters['field_speed'],slab_half_width)
                runs.append(row)
                print(f"a={a} m={m} N={cells}: photon loss={row['photon_fraction_lost']:.6g}; T={row['fractions_of_incident']['transmitted']:.6g}; delivered={row['transmitted_energy']:.6g}",flush=True)
            coarse,fine=runs[-2:]
            delta=max(abs(coarse[k]/fine[k]-1) for k in ['photon_energy_lost','incident_energy','transmitted_energy'])
            fraction_delta=abs(coarse['fractions_of_incident']['transmitted']-fine['fractions_of_incident']['transmitted'])
            convergence.append({'half_width':a,'mass':m,'max_relative_energy_change':delta,'absolute_transmission_fraction_change':fraction_delta})
            assert delta<protocol['maximum_relative_energy_change_on_refinement'],convergence[-1]
            assert fraction_delta<protocol['maximum_absolute_transmission_fraction_change'],convergence[-1]
        # Independent point-source, weak-loading limit: source velocity c=1>v.
        g_squared=quad(lambda t:np.cos(np.pi*t/(2*a))**4,-a,a,epsabs=1e-13)[0]
        assert abs(g_squared-3*a/4)<1e-12
        p=baseline['parameters'];P=p['photon_initial_energy'];K=p['field_inertia'];v=p['field_speed']
        weak_limits.append({'half_width':a,'integral_g_squared':g_squared,
                            'point_source_weak_loading_left_energy':P**2*g_squared/(4*K*(1+v)),
                            'point_source_weak_loading_right_energy':P**2*g_squared/(4*K*(1-v)),
                            'scope':'Zero smoothing and negligible feedback limit, not an exact formula for the finite-smoothing numerical run.'})
    comparisons=[]
    fine=[r for r in runs if r['cells']==max(protocol['cells'])]
    for row in fine:
        ref=next(r for r in fine if r['conversion_half_width']==2.0 and r['restoring_frequency']==row['restoring_frequency'])
        comparisons.append({'half_width':row['conversion_half_width'],'mass':row['restoring_frequency'],
                            'delivered_energy_ratio_to_width_2':row['transmitted_energy']/ref['transmitted_energy'],
                            'photon_loss_ratio_to_width_2':row['photon_energy_lost']/ref['photon_energy_lost']})
    result={'scope':protocol['scope'],'checks_pass':True,'runs':runs,'grid_convergence':convergence,
            'weak_loading_limits':weak_limits,'comparisons':comparisons,
            'full_redshift_timing_clock_law_validated':False,'permanent_deposition_derived':False,
            'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in [HERE/'check.py',HERE/'protocol.json',PARENT/'check.py',PARENT/'protocol.json']}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'source-timescale-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(comparisons,indent=2))


if __name__=='__main__':
    main()

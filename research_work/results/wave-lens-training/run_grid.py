"""Serial grid execution with all failures retained and no held-out data access."""
from pathlib import Path
import subprocess,sys,json
import numpy as np
H=Path(__file__).resolve().parent
directory=H/'parameter-grid';directory.mkdir(exist_ok=True)
states=[]
for i in range(9):
    cell=directory/f'cell-{i}';cell.mkdir(exist_ok=True)
    with (cell/'execution.log').open('w') as log:
        outcome=subprocess.run([sys.executable,str(H/'run.py'),'--grid-index',str(i)],stdout=log,stderr=subprocess.STDOUT)
    logfile=cell/'execution.log'
    logfile.write_text(logfile.read_text(),newline='\n')
    state={'cell':i,'exit_code':outcome.returncode}
    if outcome.returncode==0:
        result=json.loads((cell/'results.json').read_text())
        pred=json.loads((cell/'predictions.json').read_text())
        s=[r for r in pred if r['model']=='stationary_wave']
        assert len(s)==32 and all(r['role']=='training' for r in pred)
        state.update({'field_mass_eV_c2':result['protocol']['field_mass_eV_c2'],
                      'source_fraction':result['protocol']['source_to_stellar_mass_ratio'],
                      'scores':result['scores'],
                      'selection_score':float(.5*np.mean([np.log(r['sigma_pred_km_s']/r['sigma_observed_km_s'])**2+
                                                       np.log(r['theta_pred_arcsec']/r['theta_SIE_arcsec'])**2 for r in s]))})
    states.append(state)
    (directory/'execution-state.json').write_text(json.dumps(states,indent=2)+'\n',newline='\n')
    print(json.dumps(state),flush=True)
successful=[r for r in states if r['exit_code']==0]
summary={'grid_complete':len(successful)==9,'cells':states,
         'best_successful_cell':min(successful,key=lambda r:r['selection_score'])['cell'] if successful else None,
         'selection_status':'Training diagnostic only; no observational uncertainty likelihood or holdout validation'}
(directory/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',newline='\n')

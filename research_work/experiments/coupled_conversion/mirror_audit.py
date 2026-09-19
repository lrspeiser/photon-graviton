"""Declared source-symmetry diagnostic, preserving the original failure."""
import numpy as np
from common import HERE,hashes,read,save
from model import Model
from two_dimensional import rays,bend_metrics

def run(out):
    original=HERE/'evidence/spatial2d-v1/primary-192.json'
    base=read(original);rows=[]
    for label in ['equal-actions','reflected-sources']:
        model=Model(n=192,dim=2,L=32)
        if label=='equal-actions':model.I[:]=.00075
        else:model.X[:,1]*=-1
        model.initialize_light(energy=1.,wavelength=3.)
        initial_I=model.I.copy();initial_X=model.X.copy()
        initial=model.sectors();E0=initial.sum()
        steps=int(np.ceil(18/(.08*model.dx/np.sqrt(2))));dt=18/steps
        history=[];error=0.
        for k in range(steps+1):
            if k%4==0 or k==steps:
                sec=model.sectors();error=max(error,abs(sec.sum()/E0-1))
                if k%max(4,steps//100//4*4)==0 or k==steps:history.append([k*dt,*sec])
            if k<steps:model.step(dt)
        rr=rays(model)
        row=dict(label=label,initial_actions=initial_I,initial_positions=initial_X,
                 initial_sectors=initial,final_sectors=model.sectors(),energy_relative_error=error,
                 history=history,rays=rr,bend=bend_metrics(rr),
                 axis=model.axis[::2],phi=model.phi[::2,::2])
        save(out/(label+'.json'),row)
        rows.append({k:v for k,v in row.items() if k not in ['history','axis','phi']})
        print('MIRROR',label,row['bend']['mirror_relative'],flush=True)
    scale=max(abs(r['angle']) for r in base['rays'])
    covariance=max(abs(a['angle']+b['angle']) for a,b in zip(rows[1]['rays'],base['rays'][::-1]))/scale
    gates=dict(symmetric_sources=rows[0]['bend']['mirror_relative']<.001,
               reflected_fixture=covariance<.001,energy=max(r['energy_relative_error'] for r in rows)<.001)
    return dict(stage='mirror_audit',cases=rows,input_sha256=hashes([original]),
                reflected_covariance_error=covariance,gates=gates,numerical_pass=all(gates.values()),
                original_primary_mirror_pass=False,note='Diagnostic only; original gate unchanged.')

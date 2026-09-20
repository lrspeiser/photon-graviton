"""SE-TB1: unchanged LR3 replay with RK-stage face-energy transport."""
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from local_transfer import rhs
from discrete_flux import site_energy,face_flux,outward_power

ROOT=Path(__file__).resolve().parent


def main():
    out=ROOT/'transport-budget-v1';out.mkdir(exist_ok=False)
    files=('transport_budget.py','discrete_flux.py','local_transfer.py','local_rotor.py')
    manifest=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},radii=[1.2,1.5,2.],
                  trace_columns=['time','energy_r1.2','energy_r1.5','energy_r2','net_out_r1.2','net_out_r1.5','net_out_r2'])
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');rows=[]
    for name in ('combined','time','space'):
        earlier=ROOT/'local-transfer-v1';cfg=json.loads((earlier/(name+'.json')).read_text())
        original=np.load(earlier/(name+'.npz'));initial=original['initial'];state=initial.copy()
        n=cfg['n'];h=8/n;dt=cfg['dt'];x=np.stack(np.meshgrid(*(np.arange(n)*h-4 for _ in range(3)),indexing='ij'),axis=-1)
        masks=[np.sum(x*x,axis=-1)<r*r for r in (1.2,1.5,2.)]
        def powers(s):
            flux=face_flux(s,h,2.,1.)
            return np.array([outward_power(flux,m,h) for m in masks])
        def energies(s):
            values=site_energy(s,h,2.,1.)
            return np.array([np.sum(values[m])*h**3 for m in masks])
        total=float(np.sum(site_energy(state,h,2.,1.))*h**3)
        accumulated=np.zeros(3);trace=[[0.,*energies(state),*accumulated]];stages=[]
        p0=powers(state);print('start',name,flush=True)
        for step in range(round(1/dt)):
            k1=rhs(state,h,2.,1.);s2=state+.5*dt*k1;k2=rhs(s2,h,2.,1.)
            s3=state+.5*dt*k2;k3=rhs(s3,h,2.,1.);s4=state+dt*k3;k4=rhs(s4,h,2.,1.)
            power=np.stack([powers(s) for s in (state,s2,s3,s4)]);stages.append(power)
            accumulated+=dt*(power[0]+2*power[1]+2*power[2]+power[3])/6
            state+=dt*(k1+2*k2+2*k3+k4)/6
            trace.append([(step+1)*dt,*energies(state),*accumulated])
        trace=np.array(trace);error=float(np.max(abs(trace[:,1:4]-trace[0,1:4]+trace[:,4:7]))/total)
        replay=float(np.max(abs(state-original['final'])))
        row=dict(name=name,n=n,dt=dt,initial_total_energy=total,balance_error=error,replay_error=replay,
                 net_outward_energy=accumulated.tolist(),net_outward_fractions=(accumulated/total).tolist(),
                 initial_power=p0.tolist(),final_power=powers(state).tolist(),passed=error<1e-6 and replay<1e-12)
        np.savez_compressed(out/(name+'.npz'),initial=initial,final=state,trace=trace,stage_powers=np.array(stages))
        (out/(name+'.json')).write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(json.dumps(row),flush=True)
    error=abs(rows[0]['net_outward_energy'][1]-rows[1]['net_outward_energy'][1])/max(abs(rows[1]['net_outward_energy'][1]),1e-12)
    summary=dict(runs=rows,time_relative_difference=error,passed=all(r['passed'] for r in rows) and error<.001)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(dict(passed=summary['passed'],time_relative_difference=error)),flush=True)


if __name__=='__main__':main()

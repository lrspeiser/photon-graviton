"""Test persistence at matched endpoints, separating local and reference energy."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent


def field(mode,t,x,contrast):
    if mode=='static_fast':sign,e,de=1,.2,0.
    elif mode=='static_slow':sign,e,de=-1,.2,0.
    elif mode=='relaxing_fast':sign,e,de=1,.2*math.exp(-.2*t),-.04*math.exp(-.2*t)
    elif mode=='growing_slow':sign,e,de=-1,.2-.1*math.exp(-.2*t),.02*math.exp(-.2*t)
    elif mode=='relaxing_slow':sign,e,de=-1,.2*math.exp(-.2*t),-.04*math.exp(-.2*t)
    else:raise ValueError(mode)
    I=math.sin(math.pi*x)**2 if 0<x<1 else 0.
    Ix=math.pi*math.sin(2*math.pi*x) if 0<x<1 else 0.
    B=math.sin(math.pi*(x-3))**2 if 3<x<4 else 0.
    Bx=math.pi*math.sin(2*math.pi*(x-3)) if 3<x<4 else 0.
    return 1+sign*e*I+contrast*B,sign*e*Ix+contrast*Bx,sign*de*I


def trace(mode,launch,contrast,tol,step):
    checkpoints=np.array([-1.,0.,.5,1.,2.,3.,3.5,4.,6.])
    def rhs(x,y):
        q,qx,qt=field(mode,y[0],x,contrast)
        # y = time, log momentum, log reference energy, reference work, log dt_arr/dt_emit
        return [1/q,-qx/q,qt/(q*q),math.exp(y[2])*qt/(q*q),-qt/(q*q)]
    state=[launch,0,0,0,0]
    samples=[(-1.,np.array(state))]
    boundaries=(-1.,0.,1.,3.,4.,6.)
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        sol=solve_ivp(rhs,(lo,hi),state,method='DOP853',rtol=tol,atol=tol*.01,
                      max_step=step,dense_output=True)
        assert sol.success
        samples.extend((x,sol.sol(x)) for x in checkpoints if lo<x<=hi)
        state=sol.y[:,-1]
    rows=[]
    for x,y in samples:
        q,_,_=field(mode,y[0],x,contrast);p=math.exp(y[1]);H=math.exp(y[2])
        assert abs(q*p-H)<2e-9
        assert abs(H-1-y[3])<2e-9
        rows.append({'x':float(x),'reference_time':float(y[0]),'q':q,'local_photon_energy':p,
                     'reference_photon_energy':H,'reference_work_on_photon':float(y[3]),
                     'instantaneous_arrival_stretch':math.exp(y[4])})
    tail=[r['reference_photon_energy'] for r in rows if r['x']>=1]
    assert max(tail)-min(tail)<2e-9
    return rows


if __name__=='__main__':
    old=json.loads((HERE/'results.json').read_text())
    records=[]
    for prior in old['records']:
        mode,launch=prior['mode'],prior['launch_time']
        for contrast in (0.,.2,-.2):
            rows=trace(mode,launch,contrast,2e-12,.01)
            coarse=trace(mode,launch,contrast,2e-10,.04)
            S=1/rows[-1]['local_photon_energy']
            assert abs(S-prior['wavelength_stretch'])<2e-8
            assert abs(S-rows[-1]['instantaneous_arrival_stretch'])<2e-9
            difference=max(abs(a['local_photon_energy']-b['local_photon_energy']) for a,b in zip(rows,coarse))
            assert difference<2e-8
            records.append({'mode':mode,'launch_time':launch,'static_downstream_contrast':contrast,
                'received_redshift':S-1,'refinement_max_energy_difference':difference,'checkpoints':rows})
    out={'status':'prescribed-field arrival retention verified; companion receiver and source dynamics unresolved',
         'formula_provenance':'established Hamiltonian and clock relations applied to the prior hypothetical q profiles; no originality claim',
         'prior_results_sha256':hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest(),
         'protocol_sha256':hashlib.sha256((HERE/'arrival-protocol.md').read_bytes()).hexdigest(),
         'cases':len(records),'records':records}
    (HERE/'arrival-results.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n')
    print(json.dumps({'cases':len(records),'max_refinement_difference':max(r['refinement_max_energy_difference'] for r in records),
        'launch_zero_contrast_point_two':[{k:r[k] for k in ('mode','received_redshift')} for r in records if r['launch_time']==0 and r['static_downstream_contrast']==.2]},indent=2))

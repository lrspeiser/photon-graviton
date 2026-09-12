from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path[:0]=[str(HERE.parent/'bar-field-foundation'),str(HERE.parent/'rotating-bar-orbits')]
from fast_multipole import FastMultipole

grid=sys.argv[1];order=int(sys.argv[2])
data=json.loads((HERE/f'{grid}-L{order}.json').read_text())
assert len(data['records'])==(16 if grid=='coarse' else 36)
assert all(r['final_check']['passes'] for r in data['records'])
bar=FastMultipole.load(ROOT/f'research_work/data-cache/bar-field/bar-L{order}.npz')
times=(.05,.1,.25)
def cross_z(v):return np.array([-v[1],v[0],0.])

def follow(start,y0,tol):
    def rhs(t,y):
        _,a=bar.evaluate(y[:3])
        return np.r_[y[3:]-37.5*cross_z(y[:3]),a[0]-37.5*cross_z(y[3:])]
    def boundary(t,y):return np.linalg.norm(y[:3])-.1
    boundary.terminal=False;boundary.direction=0
    sol=solve_ivp(rhs,(start,.25),y0,method='DOP853',rtol=tol,atol=tol*.01,max_step=.0002,events=boundary,dense_output=True)
    assert sol.success
    crossing=[float(t) for t in sol.t_events[0] if t>start+1e-8]
    residence=[]
    for T in times:
        previous=start;inside=True;duration=0.
        for t in crossing+[.25]:
            if inside:duration+=max(0.,min(t,T)-previous)
            previous=t;inside=not inside
        residence.append(duration)
    ts=np.linspace(start,.25,201);y=sol.sol(ts).T
    E=.5*np.sum(y[:,3:]**2,axis=1)+bar.evaluate(y[:,:3])[0]
    L=y[:,0]*y[:,4]-y[:,1]*y[:,3];J=E-37.5*L
    return dict(crossings=crossing,inside_durations=residence,invariant_error=float(max(abs(J-J[0]))/220**2),evaluations=sol.nfev)

records=[]
for i,r in enumerate(data['records']):
    initial=r['extra_refined'] or r['refined'];start=initial['entry_time']
    record={k:r[k] for k in ('R','mu','phi','raw_weight')}
    if start is None:
        record.update(inside_durations=[0.,0.,0.],status='no entry during interval',passes=True)
    else:
        try:
            a=follow(start,initial['endpoint'],2e-9);b=follow(start,initial['endpoint'],2e-11)
            error=max(abs(x-y) for x,y in zip(a['inside_durations'],b['inside_durations']))
            record.update(coarse=a,refined=b,inside_durations=b['inside_durations'],duration_difference=error,
                          status='continued through center',passes=error<1e-5 and b['invariant_error']<1e-5)
        except ValueError as error:
            record.update(status='field domain unresolved',error=str(error),inside_durations=None,passes=False)
    records.append(record)
    print(grid,order,i+1,record['inside_durations'],record['passes'],flush=True)
    (HERE/f'residence-{grid}-L{order}.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf8',newline='\n')

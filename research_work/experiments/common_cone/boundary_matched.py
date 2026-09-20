"""OB-2 independently advected, equally damped wave characteristics."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import subprocess
import numpy as np
from run_audit import HERE,ROOT,save,digest


def integrate(n,length,width=0,strength=0):
    dx=length/n;x=np.arange(n)*dx-length/2
    z=np.maximum(1-(x/.6)**2,0)
    gradient=-8*x/.6**2*z**3
    gamma=strength*np.maximum((abs(x)-(length/2-width))/width,0)**2 if width else np.zeros(n)
    y=np.concatenate((-gradient,gradient,[0.]))
    initial=dx*np.sum(gradient**2)
    drift=0;opposite=0;steps=int(np.ceil(9/(.2*dx)));dt=9/steps
    def rhs(state):
        p=state[:n];v=state[n:2*n]
        dp=(np.roll(v,-1)-np.roll(v,1))/(2*dx)-gamma*p
        dv=(np.roll(p,-1)-np.roll(p,1))/(2*dx)-gamma*v
        return np.concatenate((dp,dv,[dx*np.sum(gamma*(p*p+v*v))]))
    for _ in range(steps):
        k1=rhs(y);k2=rhs(y+dt*k1/2);k3=rhs(y+dt*k2/2);k4=rhs(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        p=y[:n];v=y[n:2*n]
        drift=max(drift,abs(dx*np.sum(p*p+v*v)/2+y[-1]-initial)/initial)
        opposite=max(opposite,dx*np.sum((p+v)**2)/4/initial)
    return dict(x=x,state=y,initial_energy=initial,ledger_drift=drift,opposite_characteristic_energy=opposite)


def main():
    out=HERE/'boundary-matched-v1';out.mkdir(exist_ok=False)
    files=[HERE/name for name in ('boundary-protocol.md','boundary-amendment.md','boundary_matched.py','run_audit.py')]
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                 hashes={p.relative_to(ROOT).as_posix():digest(p) for p in files}))
    rows=[];profiles={};references=[]
    for n in (256,512):
        reference=integrate(3*n,36)
        profiles[f'reference-n{n}']=np.column_stack((reference['x'],reference['state'][:3*n],reference['state'][3*n:6*n]))
        references.append(dict(n=n,ledger_drift=reference['ledger_drift'],passed=bool(reference['ledger_drift']<=1e-5)))
        pref=reference['state'][n:2*n];vref=reference['state'][4*n:5*n]
        for width in (1,2,3):
            for strength in (2,8,32,128):
                result=integrate(n,12,width,strength);x=result['x'];y=result['state'];dx=12/n
                dp=y[:n]-pref;dv=y[n:2*n]-vref
                error=np.sqrt(dx*np.sum((dp*dp+dv*dv)[abs(x)<3.5])/2/result['initial_energy'])
                rows.append(dict(n=n,width=width,strength=strength,relative_interior_state_error=error,
                                 ledger_drift=result['ledger_drift'],opposite_characteristic_energy=result['opposite_characteristic_energy'],
                                 energy_removed_fraction=y[-1]/result['initial_energy'],initial_energy=result['initial_energy'],
                                 passed=bool(error<=.01 and result['ledger_drift']<=1e-5 and result['opposite_characteristic_energy']<=1e-12)))
                profiles[f'n{n}-w{width}-g{strength}']=np.column_stack((x,y[:n],y[n:2*n],dp,dv))
                save(out/'runs.json',rows)
    choices=[]
    for width in (1,2,3):
        for strength in (2,8,32,128):
            if all(r['passed'] for r in rows if r['width']==width and r['strength']==strength):choices.append(dict(width=width,strength=strength))
    summary=dict(cases=len(rows),passes=sum(r['passed'] for r in rows),references=references,
                 selected=choices[0] if choices else None,selected_rows=[r for r in rows if choices and r['width']==choices[0]['width'] and r['strength']==choices[0]['strength']],
                 maximum_opposite_characteristic_energy=max(r['opposite_characteristic_energy'] for r in rows),
                 maximum_ledger_drift=max(r['ledger_drift'] for r in rows),
                 scope='Auxiliary 1D vacuum boundary prototype; not yet integrated into CC-2.')
    np.savez_compressed(out/'profiles.npz',**profiles);save(out/'summary.json',summary)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(encoding='utf8'))


if __name__=='__main__':main()

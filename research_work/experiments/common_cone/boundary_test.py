"""OB-1 independent 1D outgoing-layer/reflection campaign."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_audit import HERE,ROOT,save,digest


def energy(f,p,dx):return .5*dx*np.sum(p*p+((np.roll(f,-1)-f)/dx)**2)


def integrate(n,length,width=0,strength=0):
    dx=length/n;x=np.arange(n)*dx-length/2
    z=np.maximum(1-(x/.6)**2,0)
    gamma=strength*np.maximum((abs(x)-(length/2-width))/width,0)**2 if width else np.zeros(n)
    y=np.concatenate((z**4,8*x/.6**2*z**3,[0.]))
    initial=energy(y[:n],y[n:2*n],dx);drift=0
    steps=int(np.ceil(9/(.2*dx)));dt=9/steps
    def rhs(state):
        f=state[:n];p=state[n:2*n]
        lap=(np.roll(f,-1)+np.roll(f,1)-2*f)/dx**2
        return np.concatenate((p,lap-gamma*p,[dx*np.sum(gamma*p*p)]))
    for _ in range(steps):
        k1=rhs(y);k2=rhs(y+dt*k1/2);k3=rhs(y+dt*k2/2);k4=rhs(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
        drift=max(drift,abs(energy(y[:n],y[n:2*n],dx)+y[-1]-initial)/initial)
    return dict(x=x,state=y,initial_energy=initial,ledger_drift=drift,dt=dt)


def main():
    out=HERE/'boundary-v1';out.mkdir(exist_ok=False)
    files=[HERE/name for name in ('boundary-protocol.md','boundary_test.py','run_audit.py')]
    save(out/'manifest.json',dict(git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                                hashes={p.relative_to(ROOT).as_posix():digest(p) for p in files}))
    rows=[];profiles={};references=[]
    for n in (256,512):
        reference=integrate(3*n,36)
        profiles[f'reference-n{n}']=np.column_stack((reference['x'],reference['state'][:3*n],reference['state'][3*n:6*n]))
        references.append(dict(n=n,ledger_drift=reference['ledger_drift'],passed=bool(reference['ledger_drift']<=1e-5)))
        fref=reference['state'][n:2*n];pref=reference['state'][4*n:5*n]
        for width in (1,2,3):
            for strength in (2,8,32,128):
                result=integrate(n,12,width,strength);y=result['state'];x=result['x'];dx=12/n
                df=y[:n]-fref;dp=y[n:2*n]-pref
                error_density=.5*(dp*dp+((np.roll(df,-1)-df)/dx)**2)
                error=np.sqrt(dx*np.sum(error_density[abs(x)<3.5])/result['initial_energy'])
                rows.append(dict(n=n,width=width,strength=strength,relative_interior_state_error=error,
                                 ledger_drift=result['ledger_drift'],energy_removed_fraction=y[-1]/result['initial_energy'],
                                 max_remaining_amplitude=np.max(abs(y[:n])),
                                 passed=bool(error<=.01 and result['ledger_drift']<=1e-5)))
                profiles[f'n{n}-w{width}-g{strength}']=np.column_stack((x,y[:n],y[n:2*n],df,dp))
                print(f'N={n}, width={width}, gamma={strength}: interior error={error:.4g}, ledger={result["ledger_drift"]:.3g}',flush=True)
                save(out/'runs.json',rows)
    choices=[]
    for width in (1,2,3):
        for strength in (2,8,32,128):
            cases=[r for r in rows if r['width']==width and r['strength']==strength]
            if all(r['passed'] for r in cases):choices.append(dict(width=width,strength=strength))
    summary=dict(cases=len(rows),passes=sum(r['passed'] for r in rows),reference_runs=references,
                 selected=choices[0] if choices else None,
                 current_cc2=[r for r in rows if r['width']==1 and r['strength']==2],
                 best_error=min(r['relative_interior_state_error'] for r in rows),
                 scope='1D vacuum compact pulse only; not nonlinear 3D boundary certification.')
    np.savez_compressed(out/'profiles.npz',**profiles)
    save(out/'summary.json',summary)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for n,marker in ((256,'o'),(512,'s')):
        for width in (1,2,3):
            cases=[r for r in rows if r['n']==n and r['width']==width]
            axes[0].loglog([r['strength'] for r in cases],[r['relative_interior_state_error'] for r in cases],marker=marker,label=f'N={n}, width={width}')
    axes[0].axhline(.01,ls='--',color='black',label='Declared limit')
    axes[0].set(xlabel='Maximum damping',ylabel='Interior relative state error',title='Reflected / wrapped pulse contamination');axes[0].legend(fontsize=7)
    current=profiles['n512-w1-g2'];axes[1].plot(current[:,0],current[:,1],label='Current CC-2 layer')
    best=min((r for r in rows if r['n']==512),key=lambda r:r['relative_interior_state_error'])
    better=profiles[f'n512-w{best["width"]}-g{best["strength"]}']
    axes[1].plot(better[:,0],better[:,1],label=f'Lowest measured error: width={best["width"]}, gamma={best["strength"]}')
    axes[1].axvspan(-3.5,3.5,color='gray',alpha=.1)
    axes[1].set(xlabel='x at time 9',ylabel='Field amplitude',title='Returning signal inside the finite box');axes[1].legend(fontsize=8)
    fig.suptitle('OB-1: outgoing-layer test against a larger-domain reference')
    fig.savefig(HERE/'boundary.png',dpi=170);plt.close(fig)
    save(out/'hashes.json',{p.name:digest(p) for p in sorted(out.iterdir()) if p.is_file()})
    print((out/'summary.json').read_text(encoding='utf8'))


if __name__=='__main__':main()

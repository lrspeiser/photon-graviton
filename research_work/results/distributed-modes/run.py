from pathlib import Path
import itertools,json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
HERE=Path(__file__).resolve().parent
V=.5;U=.003;N=48
positions=np.array(list(itertools.product([.25,.75],repeat=3)))
X0=np.repeat(positions,6,axis=0)
directions=np.tile(np.vstack([np.eye(3),-np.eye(3)]),(8,1))
grid=np.array(list(itertools.product(np.arange(8)/8,repeat=3)))
def modes(cut):
    ks=[]
    for k in itertools.product(range(-cut,cut+1),repeat=3):
        if not any(k) or np.linalg.norm(k)>cut:continue
        if next(x for x in k if x)!=abs(next(x for x in k if x)):continue
        ks.append(k)
    return 2*np.pi*np.array(ks,dtype=float)
def basis(x,k):
    phase=x@k.T
    B=np.column_stack([np.ones(len(x)),np.sqrt(2)*np.cos(phase),np.sqrt(2)*np.sin(phase)])
    G=np.concatenate([np.zeros((len(x),1,3)),-np.sqrt(2)*np.sin(phase)[:,:,None]*k[None,:,:],np.sqrt(2)*np.cos(phase)[:,:,None]*k[None,:,:]],axis=1)
    return B,G
records=[]
for amp in [0.,.5]:
  for cut in [1,2]:
    k=modes(cut);M=1+2*len(k)
    omega2=np.r_[0.,np.sum(k*k,axis=1)*V*V,np.sum(k*k,axis=1)*V*V]
    w=1+amp*np.sin(2*np.pi*X0[:,0]);w=U*w/w.sum()
    init=np.r_[np.zeros(2*M),X0.ravel(),(directions*w[:,None]).ravel()]
    def unpack(y):return y[:M],y[M:2*M],y[2*M:2*M+3*N].reshape(N,3),y[2*M+3*N:].reshape(N,3)
    def rhs(t,y):
        q,p,x,P=unpack(y);B,G=basis(x,k);n=1+B@q;mag=np.linalg.norm(P,axis=1)
        assert n.min()>0
        force=mag/n**2
        return np.r_[p,-omega2*q+B.T@force,(P/(mag*n)[:,None]).ravel(),(force[:,None]*np.einsum('ima,m->ia',G,q)).ravel()]
    runs=[]
    for tol in ([1e-9,1e-11] if cut==2 else [1e-9]):
        sol=solve_ivp(rhs,[0,40],init,method='DOP853',rtol=tol,atol=tol/1000,max_step=.1,dense_output=True)
        assert sol.success
        runs.append(sol)
    times=np.linspace(0,40,401);states=runs[-1].sol(times)
    energies=[];momenta=[];minn=1.
    Bg,_=basis(grid,k)
    for y in states.T:
        q,p,x,P=unpack(y);B,_=basis(x,k);n=1+B@q
        energies.append(.5*np.sum(p*p+omega2*q*q)+np.sum(np.linalg.norm(P,axis=1)/n))
        L=len(k)
        field_momentum=(p[1+L:]*q[1:1+L]-p[1:1+L]*q[1+L:])@k
        momenta.append(P.sum(axis=0)+field_momentum)
        minn=min(minn,float(n.min()),float((1+Bg@q).min()))
    q,p,x,P=unpack(states[:,-1]);B,_=basis(x,k)
    error=float(max(abs(np.array(energies)-U))/U)
    difference=float(np.max(abs(runs[0].sol(times)-states)))
    momentum_error=float(np.max(np.linalg.norm(np.array(momenta)-momenta[0],axis=1))/U)
    assert error<1e-6 and difference<1e-5 and minn>0 and momentum_error<1e-7
    records.append(dict(amplitude=amp,cutoff=cut,modes=M,relative_energy_error=error,tolerance_repeated=cut==2,tolerance_state_difference=difference if cut==2 else None,momentum_error_over_initial_energy=momentum_error,min_n_sampled=minn,zero_mode_rate=float(p[0]),rate_spatial_rms=float(np.std(Bg@p)),nonzero_mode_energy=float(.5*np.sum(p[1:]**2+omega2[1:]*q[1:]**2)),final_photon_energy=float(np.sum(np.linalg.norm(P,axis=1)/(1+B@q)))))
(HERE/'results.json').write_text(json.dumps(dict(protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),cases=records,scope='Finite-mode periodic 3D diagnostic; no continuum or observational validation'),indent=2)+'\n',newline='\n')
print(json.dumps(records,indent=2))

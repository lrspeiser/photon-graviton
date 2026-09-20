"""SV-1 conservative mechanical proxy, not a massless graviton theory."""
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
M=100.;INERTIA=100.;EPS=.2;N=12
PACKETS=("radial","rotating","reverse","counter","hot","stream")
def unpack(y):
    n=len(y)//6
    return y[:3*n].reshape(n,3),y[3*n:].reshape(n,3)
def pack(q,p):return np.r_[q.ravel(),p.ravel()]
def pair(q,ell):
    r=q[:,None,:]-q[None,:,:]
    w=np.exp(-np.sum(r*r,axis=2)/(2*ell**2))
    np.fill_diagonal(w,0)
    return r,w
def components(y,spec):
    q,p=unpack(y);n=len(q)-1
    r,w=pair(q[:-1],spec["ell"]);k=spec["eta"]/(n-1)
    kinetic=.5*np.sum(p[:-1]**2)
    interaction=-.5*k*np.sum(w*(p[:-1]@p[:-1].T))
    recoil=np.sum(p[-1]**2)/(2*M)
    rel=q[:-1]-q[-1]
    potential=-spec["mu"]*np.sum(1/np.sqrt(np.sum(rel*rel,axis=1)+EPS**2))
    return np.array([kinetic,interaction,recoil,potential])
def energy(y,spec):return float(components(y,spec).sum())
def rhs(t,y,spec):
    q,p=unpack(y);n=len(q)-1
    r,w=pair(q[:-1],spec["ell"]);k=spec["eta"]/(n-1)
    rel=q[:-1]-q[-1]
    source=-spec["mu"]*rel/(np.sum(rel*rel,axis=1)+EPS**2)[:,None]**1.5
    qd=np.vstack((p[:-1]-k*w@p[:-1],p[-1]/M))
    pd=-k/spec["ell"]**2*np.sum((w*(p[:-1]@p[:-1].T))[:,:,None]*r,axis=1)+source
    return pack(qd,np.vstack((pd,-source.sum(axis=0))))
def initial(packet):
    a=np.arange(N)*2*np.pi/N
    er=np.column_stack((np.cos(a),np.sin(a),np.zeros(N)))
    et=np.column_stack((-np.sin(a),np.cos(a),np.zeros(N)))
    q=er.copy();p=.8*er
    if packet=="rotating":p+=.6*et
    elif packet=="reverse":p-=.6*et
    elif packet=="counter":p+=.6*(-1.)**np.arange(N)[:,None]*et
    elif packet=="hot":
        z=1-2*(np.arange(N)+.5)/N
        a=np.arange(N)*np.pi*(3-np.sqrt(5));rr=np.sqrt(1-z*z)
        q=np.column_stack((rr*np.cos(a),rr*np.sin(a),z))
        rng=np.random.default_rng(20260919);p=rng.normal(size=(N,3));p*=.9/np.linalg.norm(p,axis=1)[:,None]
    elif packet=="stream":
        q=np.column_stack((-3+.12*np.arange(N),np.full(N,.7),np.zeros(N)))
        p=np.tile([.9,0.,0.],(N,1))
    elif packet!="radial":raise ValueError(packet)
    spin=-np.cross(q,p).sum(axis=0)
    return pack(np.vstack((q,np.zeros(3))),np.vstack((p,-p.sum(axis=0)))),spin
def metric(y,spec):
    q,p=unpack(y);v,_=unpack(rhs(0,y,spec));rel=q[:-1]-q[-1]
    vel=v[:-1]-v[-1];rad=np.linalg.norm(rel,axis=1)
    cent=q[:-1]-q[:-1].mean(axis=0);speeds=np.linalg.norm(vel,axis=1)
    _,w=pair(q[:-1],spec["ell"]);den=speeds[:,None]*speeds[None,:]
    cos=np.divide(vel@vel.T,den,out=np.zeros_like(den),where=den>1e-30)
    total=float(w.sum())
    return dict(rms_radius=float(np.sqrt(np.mean(rad**2))),
        packet_width=float(np.sqrt(np.mean(np.sum(cent**2,axis=1)))),
        retained=float(np.mean(rad<3)),
        circulation=float(np.linalg.norm(np.cross(rel,vel).sum(axis=0))/max(np.sum(rad*speeds),1e-30)),
        heading=float(np.sum(w*cos)/total) if total>1e-12 else 0.,
        neighbor_weight=total,
        minimum_kinetic_eigenvalue=float(np.linalg.eigvalsh(np.eye(len(rad))-spec["eta"]/(len(rad)-1)*w)[0]))
def invariants(y,spin):
    q,p=unpack(y)
    return p.sum(axis=0),np.cross(q,p).sum(axis=0)+spin
def controls():
    rows=[];pairs=[]
    def add(name,val,lim):rows.append(dict(name=name,value=float(val),limit=lim,passed=bool(val<=lim)))
    rng=np.random.default_rng(19)
    q=rng.normal(size=(5,3));p=rng.normal(size=(5,3));y=pack(q,p)
    for eta in (-.8,.8):
        for mu in (0.,1.):
            s=dict(eta=eta,ell=1.,mu=mu);eps=1e-6
            grad=np.array([(energy(y+eps*np.eye(len(y))[i],s)-energy(y-eps*np.eye(len(y))[i],s))/(2*eps) for i in range(len(y))])
            half=len(y)//2
            canonical=np.r_[grad[half:],-grad[:half]]
            add(f"Hamilton gradient eta={eta} mu={mu}",np.max(abs(rhs(0,y,s)-canonical))/max(1,np.max(abs(canonical))),1e-6)
            dy=rhs(0,y,s);qd,pd=unpack(dy)
            add(f"reciprocal force eta={eta} mu={mu}",np.linalg.norm(pd.sum(axis=0)),1e-12)
            add(f"angular derivative eta={eta} mu={mu}",np.linalg.norm((np.cross(qd,p)+np.cross(q,pd)).sum(axis=0)),1e-12)
            add(f"translation invariance eta={eta} mu={mu}",abs(energy(pack(q+np.array([2,3,4]),p),s)-energy(y,s)),1e-11)
            rot=np.linalg.qr(rng.normal(size=(3,3)))[0]
            add(f"rotation invariance eta={eta} mu={mu}",abs(energy(pack(q@rot,p@rot),s)-energy(y,s)),1e-11)
    s=dict(eta=0.,ell=1.,mu=0.);qd,pd=unpack(rhs(0,y,s))
    add("zero coupling free velocities",np.max(abs(qd[:-1]-p[:-1])),0)
    add("zero source and coupling free momenta",np.max(abs(pd)),0)
    w=np.ones((N,N))-np.eye(N)
    add("global positive kinetic bound",abs(np.linalg.eigvalsh(np.eye(N)-.8/(N-1)*w)[0]-.2),1e-12)
    for eta in (-.8,0.,.8):
        for heading,p2 in (("parallel",[1,0,0]),("opposite",[-1,0,0]),("orthogonal",[0,1,0])):
            q=np.array([[-.5,0,0],[.5,0,0],[0,0,0]],float)
            p=np.array([[1,0,0],p2,[0,0,0]],float);p[-1]=-p[:-1].sum(axis=0)
            z=pack(q,p);s=dict(eta=eta,ell=1.,mu=0.)
            f=rhs(0,z,s);v,dp=unpack(f)
            eps=1e-6
            ac,_=unpack((rhs(0,z+eps*f,s)-rhs(0,z-eps*f,s))/(2*eps))
            rr=q[1]-q[0];vv=v[1]-v[0];aa=ac[1]-ac[0];d=np.linalg.norm(rr)
            sepdd=(vv@vv+rr@aa)/d-(rr@vv)**2/d**3
            expected=eta*np.exp(-.5)*np.dot(p[0],p[1])
            add(f"canonical pair sign {heading} eta={eta}",abs(dp[0,0]-expected),1e-12)
            pairs.append(dict(eta=eta,heading=heading,canonical_force_left_x=float(dp[0,0]),physical_separation_acceleration=float(sepdd)))
    return rows,pairs

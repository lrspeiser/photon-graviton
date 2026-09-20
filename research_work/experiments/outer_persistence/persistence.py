"""OP-1 radial retention closure; hypothesis, not a field-energy density."""
import sys,json,hashlib
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/"coherent_memory_fit"))
import laws as L
D=L.D
BASE=D.read(L.HERE/"summary.json")["selected"]["joint"]
BASE_FORCE=L.candidate_force(BASE)
BASELINE=dict(id="CMF-baseline",A=0.,delta=0.,L0=100.,eta=0.)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare(obj,n=512):
    r=np.asarray(obj["r"]);gb=np.asarray(obj["gb"])
    if len(r)<2 or np.any(np.diff(r)<=0):raise ValueError("Need increasing positive radii")
    grid=np.unique(np.r_[r,np.geomspace(r[0],r[-1],n)])
    bg=np.interp(np.log(grid),np.log(r),gb)
    g0=BASE_FORCE(dict(r=grid,gb=bg,mass=obj["mass"])) if np.any(bg) else np.zeros_like(bg)
    return dict(r=grid,q=np.maximum(grid*(g0-bg),0),g0=g0,
                release=1/(1+(bg*D.CODE_TO_SI/1e-7)**2),
                indices=np.searchsorted(grid,r),mass=obj["mass"])
def envelope(r,q,delta,length):
    h=np.empty_like(q);h[0]=q[0]
    decay=np.exp(-delta*np.diff(np.log(r))-np.diff(r)/length)
    for i in range(1,len(q)):h[i]=max(q[i],h[i-1]*decay[i-1])
    return h
def direct_envelope(r,q,delta,length):
    return np.array([max(q[:i+1]*(r[:i+1]/r[i])**delta*np.exp(-(r[i]-r[:i+1])/length)) for i in range(len(r))])
def correction(p,spec):
    if p["mass"]<=0:return np.zeros(len(p["indices"]))
    length=spec["L0"]*(p["mass"]/1e11)**spec["eta"]
    h=envelope(p["r"],p["q"],spec["delta"],length)
    return (p["release"]*(h-p["q"])/p["r"])[p["indices"]]
def force(spec,n=512):
    def evaluate(obj):
        if spec["A"]==0:return BASE_FORCE(obj) if np.any(obj["gb"]) else np.zeros_like(obj["gb"])
        p=prepare(obj,n)
        return p["g0"][p["indices"]]+spec["A"]*correction(p,spec)
    return evaluate
def specs():
    for d in (-.5,0.,.25,.5,1.):
        for ell in (30.,100.,300.):
            for eta in (0.,1/3,.5):
                for a in (-.5,.25,.5,1.,2.):
                    yield dict(id=f"d{d:g}-L{ell:g}-eta{eta:.6g}-A{a:g}",delta=d,L0=ell,eta=eta,A=a)
def controls(gals,cs):
    checks=[]
    def add(name,value,limit):checks.append(dict(name=name,value=float(value),limit=limit,passed=bool(value<=limit)))
    rng=np.random.default_rng(20260919)
    r=np.geomspace(.1,10000,113);q=rng.uniform(0,10,len(r))
    for d in (-.5,0,1):
        a=envelope(r,q,d,300);b=direct_envelope(r,q,d,300)
        add("independent envelope delta "+str(d),np.max(abs(a-b))/max(abs(b)),1e-11)
    q=np.zeros(len(r));q[0]=7
    expected=7*(r[0]/r)**.25*np.exp(-(r-r[0])/100)
    add("analytic isolated peak decay",np.max(abs(envelope(r,q,.25,100)-expected))/7,1e-12)
    add("constant support no correction",np.max(abs(envelope(r,np.ones_like(r),.25,100)-1)),0)
    add("duplicate evaluation no accumulation",np.max(abs(envelope(r,q,0,100)-envelope(r,q,0,100))),0)
    add("radial length rescaling",np.max(abs(envelope(7*r,q,.25,700)-envelope(r,q,.25,100))),1e-11)
    add("same galaxy count",abs(len(gals)-149),0)
    add("same galaxy row count",abs(sum(len(g["r"]) for g in gals)-3152),0)
    add("same pressure row count",abs(sum(len(c["y"]) for c in cs)-246),0)
    spec=dict(A=2,delta=-.5,L0=300,eta=.5)
    z=dict(r=np.geomspace(1,100,30),gb=np.zeros(30),mass=0.)
    add("zero source",np.max(abs(force(spec)(z))),0)
    for g in gals:
        add("A zero "+g["name"],np.max(abs(force(BASELINE)(g)-BASE_FORCE(g))),1e-12)
    hi=dict(r=np.geomspace(1,100,30),gb=np.full(30,1/D.CODE_TO_SI),mass=1e11)
    add("large acceleration ordinary limit",np.max(abs(force(spec)(hi)/hi["gb"]-1)),1e-12)
    for c in cs:
        add("independent pressure "+c["name"],np.max(abs(c["Pmat"]@c["gb"]-D.pressure(c,c["gb"])))/np.max(D.pressure(c,c["gb"])),1e-10)
    return checks

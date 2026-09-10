"""Resolve a fixed smooth pulse in the finite-feedback Hamiltonian."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss

HERE=Path(__file__).resolve().parent

def solve(energy,nodes,modes=32,split=1):
    q,w=leggauss(nodes);x0=.15+.2*q
    measure=w*np.cos(np.pi*q/2)**2;measure/=measure.sum()
    x0=np.repeat(x0,split);measure=np.repeat(measure/split,split)
    N=len(x0)+1;weights=np.r_[.01,energy*measure]
    L=8.;mass0=8.;mass=4.;v=.5
    k=2*np.pi*np.arange(1,modes+1)/L;win=np.exp(-.5*(.2*k)**2)
    off=2+4*modes;duration=4.5
    def optical(y,x):
        co=np.cos(np.atleast_1d(x)[:,None]*k);si=np.sin(np.atleast_1d(x)[:,None]*k)
        n=y[0]+(co*win)@y[2:2+modes]+(si*win)@y[2+modes:2+2*modes]
        nx=(-si*win*k)@y[2:2+modes]+(co*win*k)@y[2+modes:2+2*modes]
        return n,nx,co,si
    def rhs(t,y):
        n,nx,co,si=optical(y,y[off:off+N]);a=weights*y[off+N:]/n**2
        if min(n)<=0:raise ValueError('Nonpositive optical factor')
        spring=mass*v*v*k*k
        return np.r_[y[1]/mass0,a.sum(),y[2+2*modes:off]/mass,
          -spring*y[2:2+modes]+win*(a@co),-spring*y[2+modes:2+2*modes]+win*(a@si),
          1/n,y[off+N:]*nx/n**2]
    initial=np.zeros(off+2*N);initial[0]=1;initial[1]=mass0*.05
    initial[off:off+N]=np.r_[.4,x0];initial[off+N:]=1.
    sol=solve_ivp(rhs,[0,duration],initial,dense_output=True,method='DOP853',rtol=2e-11,atol=1e-13,max_step=.02)
    assert sol.success
    te=[];to=[];fe=[];fo=[]
    for i in range(1,N):
        a=brentq(lambda t:sol.sol(t)[off+i]-.4,0,duration,xtol=1e-13)
        b=brentq(lambda t:sol.sol(t)[off+i]-2.4,a,duration,xtol=1e-13)
        ya,yb=sol.sol(a),sol.sol(b)
        te.append(a);to.append(b);fe.append(ya[off+N+i]/optical(ya,[.4])[0][0]);fo.append(yb[off+N+i]/optical(yb,[2.4])[0][0])
    te,to,fe,fo=map(np.array,(te,to,fe,fo));S=fe/fo
    def avg(x):return float(measure@x)
    et=avg(te);ot=avg(to)
    slope=avg((te-et)*(to-ot))/avg((te-et)**2)
    residual=to-(ot+slope*(te-et))
    width_e=np.sqrt(avg((te-et)**2));width_o=np.sqrt(avg((to-ot)**2))
    def energy_width(t,f):
        ww=measure*f;ww/=ww.sum();mean=ww@t
        return float(np.sqrt(ww@((t-mean)**2)))
    energy_ratio=energy_width(to,fo)/energy_width(te,fe)
    totals=[];momenta=[];field=[];bound=[]
    for t in np.linspace(0,duration,181):
        y=sol.sol(t);qc=y[2:2+modes];qs=y[2+modes:2+2*modes];pc=y[2+2*modes:2+3*modes];ps=y[2+3*modes:off]
        f=y[1]**2/(2*mass0)+np.sum((pc**2+ps**2)/(2*mass)+mass*v*v*k*k*(qc**2+qs**2)/2)
        photon=float(weights@(y[off+N:]/optical(y,y[off:off+N])[0]))
        totals.append(f+photon);field.append(f)
        momenta.append(float(weights@y[off+N:]+np.sum(k*(qc*ps-qs*pc))))
        bound.append(float(y[0]-np.sum(win*np.hypot(qc,qs))))
    ee=float(np.max(abs(np.array(totals)-totals[0]))/totals[0]);pe=float(np.max(abs(np.array(momenta)-momenta[0]))/momenta[0])
    assert ee<1e-8 and pe<1e-8 and min(bound)>0
    summary={'signal_energy':energy,'quadrature_nodes':nodes,'coincident_split':split,'modes':modes,
       'mean_carrier_stretch':avg(S),'carrier_stretch_std':float(np.sqrt(avg((S-avg(S))**2))),
       'count_weighted_duration_ratio':float(width_o/width_e),'energy_weighted_duration_ratio':energy_ratio,
       'best_affine_slope':float(slope),'affine_residual_rms_over_received_width':float(np.sqrt(avg(residual**2))/width_o),
       'source_count_time_width':float(width_e),'detector_count_time_width':float(width_o),
       'field_energy_gain':float(field[-1]-field[0]),'max_relative_energy_error':ee,'max_relative_momentum_error':pe,
       'sampled_global_optical_lower_bound':min(bound)}
    curves={'initial_x':x0.tolist(),'normalized_count_weights':measure.tolist(),'source_times':te.tolist(),'detector_times':to.tolist(),'source_frequency':fe.tolist(),'detector_frequency':fo.tolist()}
    return summary,curves

def main():
    runs=[];curves={};checks=[]
    for energy in [0.,.001,.01,.1]:
        group=[]
        for nodes in [8,16,32]:
            a,b=solve(energy,nodes);runs.append(a);group.append(a)
            if nodes==32:curves[str(energy)]=b
            print(energy,nodes,a['energy_weighted_duration_ratio'],a['affine_residual_rms_over_received_width'],flush=True)
        keys=['mean_carrier_stretch','count_weighted_duration_ratio','energy_weighted_duration_ratio','best_affine_slope','affine_residual_rms_over_received_width']
        delta={k:abs(group[-1][k]-group[-2][k]) for k in keys}
        assert max(delta.values())<1e-6
        checks.append({'energy':energy,'16_to_32_absolute_changes':delta})
    base=next(x for x in runs if x['signal_energy']==.01 and x['quadrature_nodes']==32)
    split,_=solve(.01,32,split=4);mode,_=solve(.01,32,modes=64)
    for control in [split,mode]:
        assert max(abs(control[k]-base[k]) for k in ['mean_carrier_stretch','energy_weighted_duration_ratio','field_energy_gain'])<1e-8
    out={'runs':runs,'resolved_curves':curves,'quadrature_checks':checks,'coincident_split_control':split,'field_mode_control':mode,
       'formula_provenance':'Known Hamiltonian ray dynamics, Gaussian quadrature and weighted moments. Smooth-pulse application is a project diagnostic; no novelty or observational validation claim.',
       'scope':'Fixed smooth pulse shape and energy; reference-clock convention only. No matter-clock closure, physical normalization, astronomical fit or holdout scoring.',
       'run_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'results.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')

if __name__=='__main__':main()

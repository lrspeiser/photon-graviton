"""Finite emission band from a positive fourth-spatial-derivative scalar action."""
from pathlib import Path
import json,os
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'dispersive-companion'

def omega(q,v):return q*np.sqrt(v*v+q*q) # M=1

def qlimit(e,v):
    qc=np.sqrt(1-v*v)
    return qc if e>=qc else brentq(lambda q:q+omega(q,v)-2*e,0,qc,xtol=1e-15)

def kernel(q,e,v):
    w=omega(q,v);ep=e-w
    delta=(q*q-w*w)/(2*e*ep)
    mu=1-delta
    p=(1+mu*mu)/2
    rate=ep*ep/(16*np.pi*np.sqrt(v*v+q*q))*p
    return rate,w/e,delta,2*np.arcsin(np.sqrt(np.clip(delta/2,0,1)))

def moments(e,v):
    end=qlimit(e,v)
    def integ(index):
        def f(y):
            rate,loss,delta,theta=kernel(y*end,e,v)
            return end*rate*[1,loss,loss*loss,delta,theta*theta][index]
        return quad(f,0,1,epsabs=1e-25,epsrel=1e-10)[0]
    n,loss,loss2,turn,angle2=[integ(i) for i in range(5)]
    return {'v':v,'E_over_M':e,'E_over_band_edge':e/np.sqrt(1-v*v),
            'maximum_Q_over_M':end,'Gamma_over_g_squared_M_cubed':n,
            'alpha_over_g_squared_M_cubed':loss,'mean_fractional_loss':loss/n,
            'local_loss_second_to_first_moment':loss2/loss,
            'mean_one_minus_cos_theta':turn/n,'rms_deflection_degrees':np.sqrt(angle2/n)*180/np.pi}

def main():
    previous=json.loads((ROOT/'research_work/results/slow-companion/slow-companion-results.json').read_text(encoding='utf-8'))
    rows=[];crosschecks=[];limits=[]
    for v in [.1,.5,.9]:
        qc=np.sqrt(1-v*v)
        for ratio in [.0001,.01,.1,1,10,100,10000]:
            e=ratio*qc;row=moments(e,v);double=moments(2*e,v)
            row['alpha_ratio_at_double_E']=double['alpha_over_g_squared_M_cubed']/row['alpha_over_g_squared_M_cubed']
            if ratio>=1:assert row['alpha_ratio_at_double_E']>=2-1e-10
            if ratio==10000:assert abs(row['alpha_ratio_at_double_E']-2)<.001
            rows.append(row)
        small=rows[-7];big=rows[-1]
        ref=next(r for r in previous['scenarios'] if r['v_over_photon_speed']==v)
        low=small['alpha_over_g_squared_M_cubed']/small['E_over_M']**3/ref['fractional_loss_rate_over_g_squared_E_cubed']
        high=big['alpha_over_g_squared_M_cubed']/(big['E_over_M']*qc*qc/(32*np.pi))
        assert abs(low-1)<1e-5 and abs(high-1)<.001
        limits.append({'v':v,'low_energy_ratio_to_linear_dispersion':low,'high_energy_ratio_to_E_Qc_squared_over_32pi':high})
        # Independent integration in emitted frequency rather than momentum.
        for ratio in [.1,1,10]:
            e=ratio*qc;qmax=qlimit(e,v);wmax=omega(qmax,v)
            def in_frequency(y):
                w=y*wmax
                q=np.sqrt(2*w*w/(v*v+np.sqrt(v**4+4*w*w)))
                dqdw=np.sqrt(v*v+q*q)/(v*v+2*q*q)
                return kernel(q,e,v)[0]*dqdw*wmax
            n=quad(in_frequency,0,1,epsabs=1e-18,epsrel=1e-10)[0]
            direct=moments(e,v)['Gamma_over_g_squared_M_cubed']
            crosschecks.append(abs(n/direct-1))
    assert max(crosschecks)<1e-8
    rng=np.random.default_rng(29110);errors=[]
    for _ in range(200):
        v=rng.uniform(.05,.95);e=10**rng.uniform(-2,2);q=qlimit(e,v)*rng.uniform(.001,.999)
        w=omega(q,v);ep=e-w
        cos_q=(q*q+2*e*w-w*w)/(2*e*q)
        assert -1-1e-10<=cos_q<=1+1e-10
        recovered=np.sqrt(e*e+q*q-2*e*q*cos_q)
        errors.append(abs(recovered-ep)/e)
    assert max(errors)<1e-10
    result={'scope':'Tree-level preferred-frame scalar/electric action with positive (Laplacian phi)^2 energy, M-scaled dispersion and initially empty companion modes. No astronomical fit or material/UV completion.',
            'dispersion':'omega^2=v^2 Q^2+Q^4/M^2',
            'band_edge':'Qc=M sqrt(1-v^2); emission also requires Q+omega(Q)<=2E',
            'rate':'dGamma/dQ=g^2 (E-omega)^2 Q (1+mu^2)/(32 pi omega)',
            'scenarios':rows,'limit_checks':limits,
            'checks':{'momentum_vs_frequency_integration_max_relative_error':max(crosschecks),'random_energy_momentum_max_relative_error':max(errors)},
            'high_energy_fractional_loss':'alpha=g^2 E Qc^2/(32 pi) to leading order; alpha(2E)/alpha(E)>=2 for E>=Qc within this full-band model',
            'verdict':'The dispersion scale permits both small packets and small deflections for E much greater than Qc, but it retains at least linear color dependence in that regime. It is not an achromatic conversion law.',
            'limits':['Need Lambda_UV above photon energies and Qc to trust full-band result','Fourth-spatial-derivative theory is preferred-frame, not a demonstrated causal UV completion','No nonlinear illuminated background, matter clocks, capture or gravitational stress coupling derived']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'dispersive-companion-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'limits':limits,'selected':[r for r in rows if r['v']==.5]},indent=2))

if __name__=='__main__':main()

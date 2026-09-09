"""Forward/inverse finite-band photon transitions in an isotropic scalar bath."""
from pathlib import Path
import importlib.util,json,os
import numpy as np
from scipy.integrate import quad

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(os.environ.get('PHOTON_GRAVITON_RESULTS',ROOT/'research_work/generated'))/'companion-bath'
spec=importlib.util.spec_from_file_location('balance',ROOT/'research_work/results/electromagnetic-balance/check_electromagnetic_balance.py')
balance=importlib.util.module_from_spec(spec);spec.loader.exec_module(balance)

def occupation(w,t):
    if t==0 or w/t>700:return 0.
    return 1/np.expm1(w/t)

def bare_rates(q,e,v,A,B):
    w=q*np.sqrt(v*v+q*q);d=q*q-w*w
    def channel(sign):
        ep=e+sign*w
        # ep^2 times polarization factor, avoiding cancellation near mu=1.
        return q/(64*np.pi*w)*(A*(2*ep-d/(2*e))**2+B*d*d/(4*e*e))
    return channel(-1),channel(1),w

def compute(e,v,t,A,B):
    qc,plus,minus=balance.coefficients(v)
    def integrals(kind):
        def f(y):
            q=qc*y;down,up,w=bare_rates(q,e,v,A,B);n=occupation(w,t)
            terms={'down_power':w*down*(1+n),'up_power':w*up*n,
                   'noise':w*w*(down*(1+n)+up*n),'noise0':w*w*down}
            return qc*terms[kind]
        return quad(f,0,1,epsabs=1e-20,epsrel=1e-10)[0]
    loss,gain,noise,noise0=[integrals(k) for k in ['down_power','up_power','noise','noise0']]
    J0=quad(lambda y:qc*(qc*y)*occupation((qc*y)*np.sqrt(v*v+(qc*y)**2),t)*(qc*y)*np.sqrt(v*v+(qc*y)**2),0,1,epsabs=1e-18)[0]
    def j2(y):
        q=qc*y;w=q*np.sqrt(v*v+q*q)
        return qc*q*occupation(w,t)*(q*q-w*w)*w
    J2=quad(j2,0,1,epsabs=1e-18)[0]
    expected=A*(plus(e)-J0/(4*np.pi)+J2/(16*np.pi*e*e))+B*minus(e)
    actual=(loss-gain)/e
    assert abs(actual-expected)<1e-9*max(abs(expected),(loss+gain)/e,1e-20)
    assert noise>=noise0*(1-1e-12)
    if A==0:assert abs(actual/(B*minus(e))-1)<1e-9
    return {'v':v,'E_over_Qc':e/qc,'temperature_over_Qc':t/qc,'A':A,'B':B,
            'gross_photon_loss_power':loss,'gross_photon_gain_power':gain,
            'net_fractional_loss_rate':actual,'empty_bath_fractional_loss_rate':A*plus(e)+B*minus(e),
            'local_energy_second_jump_moment':noise,'noise_over_empty_bath':noise/noise0,
            'net_companion_energy_gain_per_incident_photon_per_time':loss-gain,
            'analytic_net_rate':expected}

def main():
    rows=[]
    for A,B,label in [(1,1,'electric only'),(0,1,'minus channel only')]:
        for eratio in [1,10,100]:
            v=.5;qc=np.sqrt(1-v*v)
            for tratio in [0,.1,1,10,100]:
                row=compute(eratio*qc,v,tratio*qc,A,B);row['case']=label;rows.append(row)
    # Photon-state density E^2 is needed for detailed balance of per-photon rates.
    detailed=[]
    for v in [.1,.5,.9]:
        qc=np.sqrt(1-v*v);e=3*qc;t=.7*qc
        for q in qc*np.array([.1,.4,.8]):
            down,_,w=bare_rates(q,e,v,1,1)
            _,reverse,_=bare_rates(q,e-w,v,1,1)
            n=occupation(w,t)
            lhs=e*e*np.exp(-e/t)*down*(1+n)
            rhs=(e-w)**2*np.exp(-(e-w)/t)*reverse*n
            detailed.append(abs(lhs/rhs-1))
    assert max(detailed)<1e-12
    assert any(r['net_fractional_loss_rate']<0 for r in rows if r['case']=='electric only')
    # Thermal soft occupation makes the total event count IR sensitive even
    # though the energy-transfer and energy-noise moments remain finite.
    v=.5;qc=np.sqrt(1-v*v);e=10*qc;t=qc;ir=[]
    for cut in [1e-5,1e-7,1e-9]:
        def integrand(logq):
            q=np.exp(logq);down,_,w=bare_rates(q,e,v,1,1)
            return q*down*(1+occupation(w,t))
        count=quad(integrand,np.log(cut*qc),np.log(qc),epsabs=1e-10)[0]
        ir.append({'minimum_Q_over_Qc':cut,'downward_event_rate':count})
    predicted=e*e*t/(16*np.pi*v*v)*np.log(100)
    assert abs((ir[-1]['downward_event_rate']-ir[-2]['downward_event_rate'])/predicted-1)<1e-5
    result={'scope':'Dilute photons with isotropic incoherent companion occupation fixed independently of photon energy; tree-level constant couplings and full finite emission band E>=Qc. Temperatures are illustrative, not inferred cosmic inputs.',
            'net_rate':'alpha=A[F_plus(E)-J0/(4 pi)+J2/(16 pi E^2)]+B F_minus(E); J0=int Q n omega dQ; J2=int Q n (Q^2-omega^2) omega dQ',
            'scenarios':rows,'infrared_event_count_example':ir,
            'checks':{'dilute_photon_state_density_detailed_balance_max_relative_error':max(detailed)},
            'units':'hbar=c=M=1; A and B are coupling-squared weights in the displayed model. Absolute physical rates and temperature are not adopted.',
            'conclusions':['A populated bath adds both stimulated downward and inverse upward transitions','Minus-channel net drift is exactly unchanged by isotropic occupation in the full-band dilute-photon limit, while energy-jump noise increases','The forward-channel bath correction lowers net photon energy loss and can reverse its sign','Fixed finite bath moments change only constant and inverse-square terms, so they cannot cancel the positive E term and produce nonzero exact achromatic drift over a full-band interval'],
            'limits':['No finite-distance line-width or image prediction','Thermal soft occupation gives logarithmic forward-channel event-count sensitivity; finite-volume cutoff or inclusive treatment needed despite finite energy moments','Bath source/depletion, anisotropy, coherence and finite photon occupation need coupled transport','No permanent store, gravitational response or background history derived']}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'companion-bath-results.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'checks':result['checks'],'examples':[r for r in rows if r['E_over_Qc']==10 and r['temperature_over_Qc'] in [0,1,10]]},indent=2))

if __name__=='__main__':main()

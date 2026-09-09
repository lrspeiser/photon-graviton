from pathlib import Path
import json,hashlib
import mpmath as mp
mp.mp.dps=50
HERE=Path(__file__).resolve().parent
a=mp.mpf(1)/137
def level(N,j,n):
    b=mp.mpf(j)+mp.mpf('.5');x=a/n
    return 1/mp.sqrt(1+(x/(N-b+mp.sqrt(b*b-x*x)))**2)
def lines(n):
    return {'lyman_alpha':level(2,1.5,n)-level(1,.5,n),
            'lyman_beta':level(3,1.5,n)-level(1,.5,n),
            'fine_2p':level(2,1.5,n)-level(2,.5,n)}
initial=lines(mp.mpf(1));rows=[]
for value in ['1','1.001','1.01','1.1','2']:
    n=mp.mpf(value);current=lines(n)
    for name,energy in current.items():
        q=energy/initial[name];stretch=n*q
        power=4 if name=='fine_2p' else 2
        approximation=n**(-power)
        assert abs(q/approximation-1)<mp.mpf('.001')
        if n==1:assert stretch==1
        rows.append(dict(n=float(n),transition=name,atomic_frequency_ratio=float(q),propagation_stretch=float(n),measured_one_plus_z=float(stretch),measured_z=float(stretch-1),leading_scaling_ratio=float(approximation)))
    if n>1:assert abs(n*current['lyman_alpha']/initial['lyman_alpha']-n*current['fine_2p']/initial['fine_2p'])>mp.mpf('0.0001')
(HERE/'results.json').write_text(json.dumps(dict(protocol_sha256=hashlib.sha256((HERE/'protocol.md').read_bytes()).hexdigest(),rows=rows,checks_passed=True,scope='Illustrative unscreened Dirac-Coulomb completion; no observational data or novelty claim'),indent=2)+'\n',newline='\n')
print(json.dumps([r for r in rows if r['n']==1.1],indent=2))

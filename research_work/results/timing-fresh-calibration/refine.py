"""Execute the inherited eight-case refinement rule without replacing base results."""
from pathlib import Path
import ast,json,hashlib
import numpy as np

H=Path(__file__).resolve().parent
source=H/'run.py';body=[]
for node in ast.parse(source.read_text()).body:
    if isinstance(node,ast.For) and isinstance(node.target,ast.Tuple) and any(isinstance(t,ast.Name) and t.id=='seed' for t in node.target.elts):break
    body.append(node)
ns={'__file__':str(source)}
exec(compile(ast.Module(body=body,type_ignores=[]),str(source),'exec'),ns)
base=json.loads((H/'results.json').read_text())
assert len(base['cases'])==160 and base['hashes']==ns['hashes']
selected=[]
for family in ns['c']['families']:
    for snr in ns['c']['nominal_snr']:
        for truth in ns['c']['truth_b']:
            eligible=[r for r in base['cases'] if (r['family'],r['snr'],r['truth_b'])==(family,snr,truth)]
            selected.append(max(eligible,key=lambda r:(abs(r['fit']['best']['parameters'][1]-truth),r['label'])))
out=ns['R']/'research_work/generated/timing-fresh-refinement';out.mkdir(exist_ok=True)
hashes={**ns['hashes'],str(source.with_name('refine.py').relative_to(ns['R'])):ns['sha'](source.with_name('refine.py')),
        str((H/'results.json').relative_to(ns['R'])):ns['sha'](H/'results.json')}
cp=out/'checkpoint.json'
state=json.loads(cp.read_text()) if cp.exists() else dict(hashes=hashes,selection=[r['label'] for r in selected],cases=[])
assert state['hashes']==hashes and state['selection']==[r['label'] for r in selected]
ns['save'](cp,state)
quadrature=dict(ns['c']['quadrature'],sobol_power=13,seed=1402)
for row in selected:
    label=row['label']
    if any(r['label']==label for r in state['cases']):continue
    oldpath=ns['out']/(label+'.npz');assert ns['sha'](oldpath)==row['array_sha256']
    old=np.load(oldpath);x=old['log_width'];z=old['redshift'];oldlogs=old['event_log_likelihoods']
    rng=np.random.default_rng(row['seed']);logs=[]
    for i,s in enumerate(ns['slots']):
        scale=np.exp(rng.normal(0,ns['c']['intrinsic_log_scatter']))*(1+s['z'])**row['truth_b']
        offset=rng.uniform(-5,5);e=s['e']/row['snr']
        flux=ns['injection'].injected(s['t'],scale,offset,row['family'])+rng.normal(0,e)
        if i==0:
            check=ns['event_log_likelihood'](s['t'],flux,e,x,ns['c']['quadrature'],ns['p'])
            assert np.max(abs(check-oldlogs[0]))<1e-10
        logs.append(ns['event_log_likelihood'](s['t'],flux,e,x,quadrature,ns['p']))
        if (i+1)%10==0:print(json.dumps(dict(label=label,events=i+1,total_events=len(ns['slots']))),flush=True)
    logs=np.array(logs)
    a,b,sigma=row['fit']['best']['parameters']
    fitted=ns['fit'](x,logs,z,dict(a=a,b=b,sigma=sigma))
    aa,bb,ss=fitted['best']['parameters'];mu=aa+bb*np.log1p(z)
    centered=logs-logs.max(axis=1)[:,None];oldcentered=oldlogs-oldlogs.max(axis=1)[:,None]
    L=np.exp(centered);difference=abs(centered-oldcentered)
    per_event=ns['averaged_likelihood'](x,L*difference,mu,ss)/ns['averaged_likelihood'](x,L,mu,ss)
    delta=abs(bb-b);metric=float(np.mean(per_event))
    path=out/(label+'.npz');np.savez_compressed(path,log_width=x,event_log_likelihoods=logs,redshift=z)
    state['cases'].append(dict(label=label,base_b=b,refined_b=bb,truth_b=row['truth_b'],absolute_b_change=delta,
        posterior_weighted_mean_event_curve_change=metric,fit=fitted,array_sha256=ns['sha'](path),
        numerical_gate_pass=bool(delta<=.05 and metric<=.1)))
    ns['save'](cp,state)
    print(json.dumps(dict(completed=len(state['cases']),total=8,label=label,delta_b=delta,curve_change=metric)),flush=True)
state['scope']='Selected worst-case integration sensitivity; not a replacement coverage sample or physical validation. Width grid unchanged.'
ns['save'](H/'refinement-results.json',state)

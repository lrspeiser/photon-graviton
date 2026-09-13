"""Re-bin archived predictions; never import fitting modules or refit parameters."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
hashes={}
def read(rel):
    p=ROOT/rel
    hashes[rel]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p.read_text(encoding='utf-8')
pred=json.loads(read('research_work/results/isotropic-galaxy-transfer/model-comparison-predictions.json'))
old=json.loads(read('research_work/results/isotropic-galaxy-transfer/model-comparison-results.json'))
cat={}
for line in read('temporal_candidate_audit/data/SPARC_Lelli2016c.mrt').splitlines():
    f=line.split()
    if len(f)!=19:continue
    try:
        L,rd,sb,hi=map(float,(f[7],f[11],f[12],f[13]))
    except ValueError:continue
    cat[f[0]]=(rd,sb,1.33*hi/(1.33*hi+.5*L))
models=['baryons','companion_third','MOND_simple_fitted','NFW_shared_scaling']
pred=[p for p in pred if p['model'] in models]
rng=np.random.default_rng(13092026)
def stats(values):
    n=len(values)
    if not n:return None
    bias=np.array([x['bias'] for x in values])
    boot=np.mean(bias[rng.integers(0,n,(1000,n))],axis=1)
    return dict(galaxies=n,points=sum(x['n'] for x in values),bias_kms=float(bias.mean()),
        bias_bootstrap95_kms=np.quantile(boot,[.025,.975]).tolist(),
        mean_fractional_residual=float(np.mean([x['frac'] for x in values])),
        RMSE_kms=float(np.sqrt(np.mean([x['mse'] for x in values]))),
        log_RMS=float(np.sqrt(np.mean([x['log_mse'] for x in values]))))
groups={}; galaxy=[]; paired=[]
for p in pred:
    r=np.array(p['R_kpc']);y=np.array(p['observed_kms']);v=np.array(p['predicted_kms']);err=v-y
    rd,sb,gas=cat[p['galaxy']]
    radial=np.digitize(r/rd,[1,3]);bs=int(np.digitize(sb,[100,500]));bg=int(np.digitize(gas,[.2,.5]))
    masks=[('all','all',np.ones(len(r),bool))]
    masks += [('radius',str(k),radial==k) for k in range(3)]
    masks += [('surface_brightness',str(bs),np.ones(len(r),bool)),('gas_proxy',str(bg),np.ones(len(r),bool))]
    for axis,b,mask in masks:
        if not mask.any():continue
        value=dict(n=int(mask.sum()),bias=float(err[mask].mean()),frac=float((err[mask]/y[mask]).mean()),mse=float(np.mean(err[mask]**2)),log_mse=float(np.mean(np.log10(v[mask]/y[mask])**2)))
        for split in [p['split'],'all_exposed']:
            groups.setdefault((p['model'],split,axis,b),[]).append(value)
        if axis=='all':galaxy.append(dict(model=p['model'],galaxy=p['galaxy'],split=p['split'],**value))
    if p['model']=='companion_third' and (radial==0).any() and (radial==2).any():
        inner=float(np.mean(err[radial==0]/y[radial==0]));outer=float(np.mean(err[radial==2]/y[radial==2]))
        paired.append(dict(galaxy=p['galaxy'],split=p['split'],inner=inner,outer=outer,outer_minus_inner=outer-inner))
summary=[dict(model=k[0],split=k[1],axis=k[2],bin=k[3],**stats(v)) for k,v in groups.items()]
for prev in old['summary']:
    if prev['model'] not in models:continue
    new=next(x for x in summary if (x['model'],x['split'],x['axis'])==(prev['model'],prev['split'],'all'))
    assert np.isclose(new['RMSE_kms'],prev['RMSE_kms'],rtol=1e-12)
    assert np.isclose(new['log_RMS'],prev['log_RMS'],rtol=1e-12)
assert sum(x['n'] for x in galaxy if x['model']=='companion_third')==3150
assert sum(x['model']=='companion_third' for x in galaxy)==149
lens=[]
comp=json.loads(read('research_work/results/companion-extensions/orbit-transition-results.json'))['rows']
star=json.loads(read('research_work/results/companion-extensions/stellar-only-lens-control-results.json'))['rows']
for c in comp:
    s=next(x for x in star if (x['Name'],x['population'])==(c['Name'],c['population']))
    assert c['observed_stellar_vrms']==s['observed_stellar_vrms']
    assert np.isclose(s['matched_companion_all_motion_chi2'],c['all_motion_chi2'])
    lens.append(dict(name=c['Name'],population=c['population'],companion_chi2=c['all_motion_chi2'],stellar_chi2=s['all_motion_chi2'],difference=c['all_motion_chi2']-s['all_motion_chi2'],companion_bias_kms=float(np.mean(np.array(c['predicted_stellar_vrms'])-c['observed_stellar_vrms'])),stellar_bias_kms=float(np.mean(np.array(s['predicted_stellar_vrms'])-s['observed_stellar_vrms']))))
contributions={}
for model in models:
    gg=[x for x in galaxy if x['model']==model]
    contributions[model]={}
    for loss in ['mse','log_mse']:
        ordered=sorted(gg,key=lambda x:x[loss],reverse=True);total=sum(x[loss] for x in ordered)
        contributions[model][loss]=dict(top5_fraction=sum(x[loss] for x in ordered[:5])/total,top10_fraction=sum(x[loss] for x in ordered[:10])/total,top10=[dict(galaxy=x['galaxy'],split=x['split'],share=x[loss]/total,bias_kms=x['bias']) for x in ordered[:10]])
read('research_work/results/companion-extensions/shared-residual-map-protocol.md')
out=dict(scope='All data exposed; zero fitting; fixed binning; bootstrap is across-galaxy descriptive uncertainty only.',source_sha256=hashes,summary=summary,per_galaxy=galaxy,paired_radial=paired,loss_contributions=contributions,lens=lens)
(HERE/'shared-residual-map-results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')
for x in summary:
    if x['model']=='companion_third' and x['split']=='all_exposed':print(x)
print('PAIRED',len(paired),np.mean([x['outer_minus_inner'] for x in paired]),sum(x['outer']>x['inner'] for x in paired))
print('LOSS',contributions['companion_third'])
print('LENS',lens)

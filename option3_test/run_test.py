"""Reproduce with Python + numpy/scipy/matplotlib: run_test.py freeze; run_test.py score.
Original data, protocol, frozen outputs and checksum records accompany this script.
"""
from pathlib import Path
import csv, hashlib, io, json, sys, zipfile
import numpy as np
from scipy.optimize import minimize_scalar, minimize

ROOT = Path(__file__).resolve().parent
C = 299792458.0
LY = 9.4607304725808e15
KPC = 3.085677581491367e19
K = 0.000077315 / (1e6 * LY)
AK = C*C*K
CANDIDATES = ['constant', 'radial', 'square_root', 'power']
MODELS = ['ordinary'] + CANDIDATES + ['RAR']

def dump(name, data):
    (ROOT/name).write_text(json.dumps(data, indent=2, allow_nan=False))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def writecsv(name, rows):
    with (ROOT/name).open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def load_data():
    meta={}
    for line in (ROOT/'raw/SPARC_Lelli2016c.mrt').read_text().splitlines():
        try:
            # Downloaded table uses whitespace separators beyond advertised byte widths.
            fields=line.split()
            if len(fields)!=19: continue
            name=fields[0]
            meta[name]=dict(distance=float(fields[2]), method=int(fields[4]),
                            inc=float(fields[5]), rd=float(fields[11]), q=int(fields[17]))
        except (ValueError, IndexError):
            continue
    data={}; exclusions=[]; counts=[]
    with zipfile.ZipFile(ROOT/'raw/Rotmod_LTG.zip') as z:
        for fn in sorted(z.namelist()):
            name=fn.replace('_rotmod.dat',''); m=meta[name]
            a=np.atleast_2d(np.loadtxt(io.BytesIO(z.read(fn))))
            if m['q']>2 or m['inc']<30 or m['rd']<=0:
                exclusions.append(dict(galaxy=name,reason='quality/inclination/disk-scale cut')); continue
            vb2=a[:,3]*abs(a[:,3])+0.5*a[:,4]*abs(a[:,4])+0.7*a[:,5]*abs(a[:,5])
            good=np.isfinite(a).all(axis=1)&(a[:,0]>0)&(a[:,1]>0)&(a[:,2]>0)&(vb2>0)
            if good.sum()<5:
                exclusions.append(dict(galaxy=name,reason='fewer than five valid rows')); continue
            a=a[good]; vb2=vb2[good]; r=a[:,0]*KPC
            data[name]=dict(**m, name=name,r=r,v=a[:,1],err=a[:,2],gb=vb2*1e6/r)
            counts.append(dict(galaxy=name,input_rows=len(good),retained_rows=int(good.sum())))
    return data,exclusions,counts

def predict(model, params, d):
    x=d['gb']/AK; A=10**params[0] if len(params) else 0
    if model=='ordinary': f=np.zeros_like(x)
    elif model=='constant': f=np.full_like(x,A)
    elif model=='radial': f=A/(1+d['r']/(d['rd']*KPC))
    elif model=='square_root': f=A*np.sqrt(x)
    elif model=='power': f=A*x**params[1]
    elif model=='RAR':
        return np.sqrt(d['r']*d['gb']/(-np.expm1(-np.sqrt(x/A))))/1000
    return np.sqrt(d['r']*(d['gb']+AK*f))/1000

def loss(model,p,ds):
    return np.mean([np.mean(np.log10(predict(model,p,d)/d['v'])**2) for d in ds])

def metrics(ds,preds,outer=False):
    lm=[]; km=[]; frac=[]; n=0
    for d in ds:
        mask=d['r']>=d['r'].min()+2/3*(d['r'].max()-d['r'].min()) if outer else np.ones(len(d['r']),bool)
        y=d['v'][mask]; v=preds[d['name']][mask]
        lm.append(np.mean(np.log10(v/y)**2)); km.append(np.mean((v-y)**2)); frac.append(np.mean(abs(v/y-1))); n+=len(y)
    return dict(galaxies=len(ds),points=n,log_speed_rmse_dex=float(np.sqrt(np.mean(lm))),
                rmse_kms=float(np.sqrt(np.mean(km))),mean_abs_fraction=float(np.mean(frac)))

def freeze():
    if (ROOT/'frozen.json').exists(): raise RuntimeError('Refusing to overwrite existing freeze')
    data,exc,counts=load_data()
    names=sorted(data,key=lambda n:hashlib.sha256(('option3-v1:'+n).encode()).hexdigest())
    nt=int(.6*len(names)); nv=int(.2*len(names))
    split=dict(train=names[:nt],validation=names[nt:nt+nv],test=names[nt+nv:])
    train=[data[n] for n in split['train']]; val=[data[n] for n in split['validation']]
    fits={'ordinary':[]}; diagnostics={}
    for model in MODELS[1:]:
        if model=='power':
            opts=[minimize(lambda p:loss(model,p,train),[a,p],bounds=[(-5,1),(0,1)],method='L-BFGS-B',
                          options={'ftol':1e-14,'gtol':1e-9}) for a in [-2,-.5] for p in [.2,.5,.8]]
            opt=min(opts,key=lambda o:o.fun); fits[model]=opt.x.tolist()
        else:
            opt=minimize_scalar(lambda a:loss(model,[a],train),bounds=(-5,1),method='bounded',options={'xatol':1e-10})
            fits[model]=[float(opt.x)]
        diagnostics[model]=dict(success=bool(opt.success),message=str(opt.message))
    vals={m:float(np.sqrt(loss(m,fits[m],val))) for m in MODELS}
    selected=min(CANDIDATES,key=vals.get)
    frozen=dict(k_per_Mly=.000077315,aK_SI=AK,split=split,params=fits,selected=selected,
                validation_scores=vals,optimizer=diagnostics,
                hashes={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'protocol.md',ROOT/'run_test.py',*sorted((ROOT/'raw').iterdir())]})
    # Target-free file: generate test predictions before scoring or displaying test residuals.
    rows=[]
    for n in split['test']:
        d=data[n]; preds={m:predict(m,fits[m],d) for m in MODELS}
        for i,r in enumerate(d['r']):
            rows.append(dict(galaxy=n,index=i,radius_kpc=r/KPC,gbar_SI=d['gb'][i],**{m:preds[m][i] for m in MODELS}))
    writecsv('frozen_test_predictions.csv',rows)
    frozen['prediction_sha256']=sha(ROOT/'frozen_test_predictions.csv')
    dump('frozen.json',frozen); dump('exclusions.json',exc); writecsv('row_audit.csv',counts)
    print(json.dumps(dict(eligible=len(names),split_counts={s:len(v) for s,v in split.items()},selected=selected,
                         parameters=fits,validation_scores=vals,freeze_sha256=sha(ROOT/'frozen.json')),indent=2))

def score():
    frozen=json.loads((ROOT/'frozen.json').read_text()); data,_,_=load_data()
    for name,h in frozen['hashes'].items():
        assert sha(ROOT/name)==h, 'Changed frozen input or code: '+name
    assert sha(ROOT/'frozen_test_predictions.csv')==frozen['prediction_sha256']
    selected=frozen['selected']; fits=frozen['params']; splits=frozen['split']
    test=[data[n] for n in splits['test']]
    # Read and score frozen predictions, not fresh fitted values.
    preds={m:{d['name']:np.empty(len(d['r'])) for d in test} for m in MODELS}
    for row in csv.DictReader((ROOT/'frozen_test_predictions.csv').open()):
        for m in MODELS: preds[m][row['galaxy']][int(row['index'])]=float(row[m])
    summary={}; pergal=[]
    independent=[d for d in test if d['method'] in [2,3,5]]
    for m in MODELS:
        summary[m]={s:metrics([data[n] for n in names],{n:predict(m,fits[m],data[n]) for n in names})
                    for s,names in splits.items() if s!='test'}
        summary[m]['test']=metrics(test,preds[m]); summary[m]['test_outer']=metrics(test,preds[m],True)
        summary[m]['test_independent_distance']=metrics(independent,preds[m]) if independent else None
        for d in test: pergal.append(dict(galaxy=d['name'],model=m,distance_method=d['method'],**metrics([d],preds[m])))
    rng=np.random.default_rng(20260908); idx=rng.integers(0,len(test),(2000,len(test)))
    loss_by={m:np.array([np.mean(np.log10(preds[m][d['name']]/d['v'])**2) for d in test]) for m in MODELS}
    bootstrap={}
    for m in ['ordinary','RAR']:
        delta=np.sqrt(loss_by[m][idx].mean(axis=1))-np.sqrt(loss_by[selected][idx].mean(axis=1))
        bootstrap[m+'_minus_selected']=dict(delta_dex=float(np.sqrt(loss_by[m].mean())-np.sqrt(loss_by[selected].mean())),
                                              ci95=np.quantile(delta,[.025,.975]).tolist())
    clock=[]
    for d in test:
        name=d['name']; extra=(preds[selected][name]*1000)**2/d['r']-d['gb']
        order=np.argsort(d['r']); integral=np.trapezoid(extra[order]/C**2,d['r'][order])
        clock.append(dict(galaxy=name,rmin_ly=d['r'].min()/LY,rmax_ly=d['r'].max()/LY,
                          delta_ln_T=float(integral),clock_change_ppm=float(np.expm1(integral)*1e6)))
    result=dict(selected=selected,summary=summary,bootstrap=bootstrap,
                test_distance_method_counts={str(m):sum(d['method']==m for d in test) for m in range(1,6)},
                clock_ppm_median=float(np.median([r['clock_change_ppm'] for r in clock])),
                freeze_sha256=sha(ROOT/'frozen.json'))
    dump('results.json',result);writecsv('test_galaxy_scores.csv',pergal);writecsv('test_clock_gradients.csv',clock)
    scored=[]
    for d in test:
        for i,r in enumerate(d['r']):
            scored.append(dict(galaxy=d['name'],radius_kpc=r/KPC,radius_ly=r/LY,observed_kms=d['v'][i],error_kms=d['err'][i],
                               **{m:preds[m][d['name']][i] for m in MODELS}))
    writecsv('test_predictions_with_observations.csv',scored)
    make_plots(test,preds,selected,result)
    print(json.dumps(result,indent=2))

def make_plots(test,preds,selected,result):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    rows=int(np.ceil(len(test)/4)); fig,axs=plt.subplots(rows,4,figsize=(15,rows*2.7),squeeze=False)
    for ax,d in zip(axs.flat,sorted(test,key=lambda d:d['name'])):
        name=d['name'];x=d['r']/LY/1000
        ax.errorbar(x,d['v'],yerr=d['err'],fmt='.',color='#263343',ms=3,alpha=.8,label='Observed')
        for m,color,label in [('ordinary','#999999','Ordinary matter'),(selected,'#007e87','Time candidate'),('RAR','#bb6b20','RAR benchmark')]:
            ax.plot(x,preds[m][name],color=color,lw=1.6,label=label)
        ax.set_title(name,loc='left',fontsize=10);ax.set_ylim(bottom=0);ax.set_xlabel('Radius (thousand light-years)');ax.set_ylabel('km/s')
    for ax in list(axs.flat)[len(test):]: ax.set_visible(False)
    handles,labels=axs.flat[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',ncol=4,bbox_to_anchor=(.5,.975))
    fig.suptitle('Every held-out galaxy · frozen predictions',y=.996,fontsize=17)
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(ROOT/'all_holdout_curves.png',dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4.5))
    labels=['Ordinary matter','Constant effect','Radius effect','Square-root effect','Power-law effect','Established RAR']
    values=[result['summary'][m]['test']['mean_abs_fraction']*100 for m in MODELS]
    ax.barh(labels,values,color=['#9da5af']+['#007e87' if m==selected else '#88bcc0' for m in MODELS[1:-1]]+['#bb6b20'])
    ax.invert_yaxis();ax.set_xlabel('Mean absolute speed error (%) · each galaxy weighted equally')
    ax.set_title(f"Prediction errors on {len(test)} held-out galaxies",loc='left',pad=16)
    for i,v in enumerate(values):ax.text(v+.5,i,f'{v:.1f}%',va='center')
    ax.set_xlim(0,max(values)*1.2);fig.tight_layout();fig.savefig(ROOT/'holdout_comparison.png',dpi=160);plt.close(fig)

if __name__=='__main__':
    {'freeze':freeze,'score':score}[sys.argv[1]]()

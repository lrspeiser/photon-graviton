from pathlib import Path
import sys,json,hashlib,csv
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.integrate import quad
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'option3_test'))
import run_test as old

def dump(name,value): (ROOT/name).write_text(json.dumps(value,indent=2,allow_nan=False))
def csvout(name,rows):
    with (ROOT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def cliff(d,width):
    x=np.log(d['r']);y=np.log(d['gb']);n=len(x);out=[]
    assert np.all(np.diff(x)>0)
    for i in range(n):
        lo=max(0,min(i-width//2,n-width));xs=x[lo:lo+width]-x[i]
        slope=np.polynomial.polynomial.polyfit(xs,y[lo:lo+width],2)[1]
        s=max(0.,-slope);out.append(s/(1+s))
    return np.array(out)

def pred(model,p,d):
    if model in ['ordinary','RAR']:return old.predict(model,p,d)
    A=10**p[0];exponent=p[1]
    effect=A*old.AK*(d['gb']/old.AK)**exponent
    if model.startswith('cliff'):effect*=np.exp(p[2]*d[model])
    return np.sqrt(d['r']*(d['gb']+effect))/1000

def objective(model,p,ds):
    return np.mean([np.mean(np.log10(pred(model,p,d)/d['v'])**2) for d in ds])

def fit(model,ds):
    if model=='ordinary':return [],True
    if model=='RAR':
        opt=minimize_scalar(lambda p:objective(model,[p],ds),bounds=(-5,1),method='bounded',options={'xatol':1e-10})
        return [float(opt.x)],bool(opt.success)
    bounds=[(-5,1),(0,1)]+([(-np.log(2),np.log(2))] if model.startswith('cliff') else [])
    starts=[[-.6,p]+([b] if model.startswith('cliff') else []) for p,b in [(.3,-.2),(.5,0),(.7,.2)]]
    opts=[minimize(lambda p:objective(model,p,ds),start,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-8}) for start in starts]
    opt=min(opts,key=lambda o:o.fun)
    return opt.x.tolist(),bool(opt.success)

def main():
    data,_,_=old.load_data();data={n:d for n,d in data.items() if d['method'] in [2,3,5]}
    names=sorted(data,key=lambda n:hashlib.sha256(('option3-cliff-v1:'+n).encode()).hexdigest())
    folds={n:i%5 for i,n in enumerate(names)}
    for d in data.values():
        d['cliff5']=cliff(d,5);d['cliff3']=cliff(d,3)
    models=['ordinary','power','cliff5','cliff3','RAR'];fits={};predictions={m:{} for m in models}
    for fold in range(5):
        train=[data[n] for n in names if folds[n]!=fold];test=[data[n] for n in names if folds[n]==fold]
        fits[str(fold)]={}
        for m in models:
            p,success=fit(m,train);fits[str(fold)][m]=dict(params=p,success=success)
            for d in test:predictions[m][d['name']]=pred(m,p,d)
    # All predictions precede aggregate evaluation.
    dump('fold_models.json',dict(folds=folds,fits=fits,protocol_sha=old.sha(ROOT/'protocol.md'),code_sha=old.sha(ROOT/'run_cliff.py')))
    rows=[]
    for n in names:
        d=data[n]
        for i,r in enumerate(d['r']):rows.append(dict(galaxy=n,fold=folds[n],r_kpc=r/old.KPC,gbar_SI=d['gb'][i],cliff_feature5=d['cliff5'][i],cliff_feature3=d['cliff3'][i],**{m:predictions[m][n][i] for m in models}))
    csvout('oof_predictions.csv',rows)
    ds=[data[n] for n in names]
    summary={m:dict(all=old.metrics(ds,predictions[m]),outer=old.metrics(ds,predictions[m],True)) for m in models}
    rng=np.random.default_rng(20260909);idx=rng.integers(0,len(ds),(5000,len(ds)))
    losses={m:np.array([np.mean(np.log10(predictions[m][d['name']]/d['v'])**2) for d in ds]) for m in models}
    comparisons={}
    for a,b in [('power','cliff5'),('power','cliff3'),('RAR','cliff5')]:
        deltas=np.sqrt(losses[a][idx].mean(axis=1))-np.sqrt(losses[b][idx].mean(axis=1))
        comparisons[a+'_minus_'+b]=dict(delta=float(np.sqrt(losses[a].mean())-np.sqrt(losses[b].mean())),ci95=np.quantile(deltas,[.025,.975]).tolist())
    pergal=[]
    for d in ds:
        for m in models:pergal.append(dict(galaxy=d['name'],fold=folds[d['name']],model=m,**old.metrics([d],predictions[m])))
    csvout('galaxy_scores.csv',pergal)
    beta={m:[fits[str(i)][m]['params'][2] for i in range(5)] for m in ['cliff5','cliff3']}
    result=dict(galaxies=len(ds),points=sum(len(d['r']) for d in ds),fold_sizes=[sum(v==i for v in folds.values()) for i in range(5)],summary=summary,comparisons=comparisons,beta=beta,
                success=all(v['success'] for fold in fits.values() for v in fold.values()))
    dump('results.json',result)
    controls=consistency();dump('consistency.json',controls)
    plot(result)
    print(json.dumps(dict(cv=result,consistency=controls),indent=2))

def consistency():
    # Analytic flattened baryonic potential Phi=-GM/sqrt(R^2+(a+sqrt(z^2+b^2))^2).
    G=6.67430e-11;M=5e10*1.98847e30;a=3*old.KPC;b=.3*old.KPC
    prior=json.loads((old.ROOT/'frozen.json').read_text())['params']['power'];A=10**prior[0];p=prior[1]
    def gradient(R,z,extra):
        B=np.sqrt(z*z+b*b);den=(R*R+(a+B)**2)**1.5
        v=G*M/den*np.array([R,(a+B)*z/B]);mag=np.linalg.norm(v)
        return v*(A*(mag/old.AK)**(p-1)) if extra else v
    r0,r1=3*old.KPC,15*old.KPC;z0,z1=.3*old.KPC,3*old.KPC
    def integral(extra,tol):
        def radial(z):return quad(lambda u:gradient(u*old.KPC,z,extra)[0]*old.KPC/old.C**2,3,15,epsabs=tol,epsrel=tol)[0]
        def vertical(r):return quad(lambda u:gradient(r,u*old.KPC,extra)[1]*old.KPC/old.C**2,.3,3,epsabs=tol,epsrel=tol)[0]
        path1=radial(z0)+vertical(r1);path2=vertical(r0)+radial(z1)
        return dict(radial_then_vertical=path1,vertical_then_radial=path2,closed_loop=path1-path2,
                    fractional_disagreement=(path1-path2)/((path1+path2)/2))
    result=dict(toy_extra=integral(True,1e-12),toy_extra_tighter=integral(True,1e-15),newtonian_control=integral(False,1e-15))
    # Stationary endpoint factors telescope through arbitrarily many cliff intervals.
    T=np.exp(np.array([0.,-3e-7,2e-7,-1e-7,0.]))
    result['static_clock_chain']=dict(product_1_plus_z=float(np.prod(T[1:]/T[:-1])),endpoint_ratio=float(T[-1]/T[0]),
                                      target_100_Mly=float(np.exp(.000077315*100)))
    return result

def plot(result):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(9,4.5))
    models=['ordinary','power','cliff5','cliff3','RAR'];labels=['Ordinary matter','Prior form, refitted','Cliff extension (primary)','Cliff extension (sensitivity)','Established RAR']
    vals=[result['summary'][m]['all']['mean_abs_fraction']*100 for m in models]
    ax.barh(labels,vals,color=['#a5abb4','#528da3','#007e87','#89b6b7','#bd742d']);ax.invert_yaxis()
    for i,v in enumerate(vals):ax.text(v+.5,i,f'{v:.1f}%',va='center')
    ax.set_xlim(0,max(vals)*1.18);ax.set_xlabel('Mean absolute speed error (%) · equal galaxy weight')
    ax.set_title('42 galaxies with stellar or supernova distances\nExploratory five-fold predictions',loc='left',pad=14)
    ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(ROOT/'cliff_comparison.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()

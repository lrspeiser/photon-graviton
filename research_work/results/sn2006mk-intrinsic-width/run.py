from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import differential_evolution
H=Path(__file__).resolve().parent;R=H.parents[2]
loader=H.parent/'sn2006mk-phase-profile/run.py';s=loader.read_text(encoding='utf-8');cut='\nresults=[]\n';assert s.count(cut)==1
ns={'__file__':str(loader),'__name__':'likelihood_loader'};exec(compile(s.split(cut)[0],str(loader),'exec'),ns);evaluate=ns['evaluate'];cache=ns['cache'];epoch=ns['epoch'];t=ns['t']
bounds=[(-10.,5.),(1.,2.)]
def loss(x):return evaluate(x[0],x[1])['score']
runs=[]
for seed in [6101,6102]:
 opt=differential_evolution(loss,bounds,seed=seed,tol=1e-10,atol=1e-9,popsize=15,polish=True);assert opt.success
 runs.append(dict(seed=seed,phase=float(opt.x[0]),effective_stretch=float(opt.x[1]),score=float(opt.fun)))
best=min(runs,key=lambda x:x['score']);assert max(x['score'] for x in runs)-best['score']<1e-5
q=best['phase'];B=best['effective_stretch'];fit=evaluate(q,B)
# Both propagation A and intrinsic width w affect the same template phase index.
# Preserve exact predictions and anchor covariances by holding A*w fixed.
decompositions=[]
for A in [1.,1.4754,1.5957062624377376]:
 width=B/A;other=evaluate(q,A*width);maxdiff=0.
 for c in cache:
  pp=q+(c['dates']-epoch)/B;aa=q+(c['dates']-epoch)/(A*width);maxdiff=max(maxdiff,float(np.max(abs(pp-aa))))
 assert abs(other['score']-fit['score'])<1e-10 and maxdiff<1e-10
 decompositions.append(dict(propagation_stretch=A,intrinsic_source_width=width,effective_stretch=A*width,score=other['score'],max_phase_difference=maxdiff))
out=dict(scope='Conditional effective width fit and exact source/propagation degeneracy on exposed photometry; not source-population validation',bounds=dict(template_start_phase=bounds[0],effective_stretch=bounds[1]),optimizer_runs=runs,best_fit=fit,effective_stretch=B,decompositions=decompositions,interior=bool(-10<q<5 and 1<B<2),hashes={str(p.relative_to(R)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [loader,Path(__file__)]})
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps(out))

"""Component-specific attached loading: nonnegative common coefficients."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.optimize import linprog
H=Path(__file__).resolve().parent;R=H.parents[2]
p=H.parent/'baryon-attached-deposits/results.json';d=json.loads(p.read_text(encoding='utf-8'))
for name,value in d['source_hashes'].items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==value
source=H.parent/'full-bar-source-audit/results.json';targets=json.loads(source.read_text(encoding='utf-8'))['rows'];rows=d['rows'];assert len(rows)==len(targets)==240
for a,b in zip(rows,targets):assert all(a[k]==b[k] for k in ['R_kpc','z_kpc','phi_rad'])
names=list(rows[0]['components_Msun_kpc3']);X=np.array([[r['components_Msun_kpc3'][k] for k in names] for r in rows]);assert np.all(X>=0)
subsets={'full':np.ones(len(rows),bool),'inner':np.array([r['R_kpc']<=8 and abs(r['z_kpc'])<=1 for r in rows]),'midplane':np.array([r['z_kpc']==0 for r in rows])};results=[]
for resolution in ['fine','finer']:
 y=np.array([r['rho_'+resolution+'_Msun_kpc3'] for r in targets]);assert np.all(y>0)
 for name,mask in subsets.items():
  A=X[mask]/y[mask,None];scale=A.max(axis=0);B=A/scale[None,:];n,m=B.shape
  lhs=np.vstack([np.c_[B,-np.ones(n)],np.c_[-B,-np.ones(n)]]);rhs=np.r_[np.ones(n),-np.ones(n)]
  fit=linprog(np.r_[np.zeros(m),1.],A_ub=lhs,b_ub=rhs,bounds=[(0,None)]*m+[(0,1)],method='highs')
  assert fit.success
  coeff=fit.x[:m]/scale;pred=X@coeff;rel=pred/y-1
  assert abs(np.max(abs(rel[mask]))-fit.fun)<1e-7
  dual=float(rhs@fit.ineqlin.marginals);assert abs(dual-fit.fun)<1e-7
  results.append(dict(resolution=resolution,fit_subset=name,n_fit=int(mask.sum()),coefficient_by_component=dict(zip(names,coeff.tolist())),minimum_worst_relative_error=float(fit.fun),linear_program_dual_objective=dual,inner_worst_error=float(np.max(abs(rel[subsets['inner']]))),full_worst_error=float(np.max(abs(rel))),excluded_worst_error=float(np.max(abs(rel[~mask]))) if np.any(~mask) else None,fit_active_rank=int(np.linalg.matrix_rank(B[:,fit.x[:m]>1e-10])),predictions=[dict(index=i,target=float(y[i]),predicted=float(pred[i]),relative_error=float(rel[i]),used_for_fit=bool(mask[i])) for i in range(len(y))]))
  print(json.dumps({k:results[-1][k] for k in ['resolution','fit_subset','minimum_worst_relative_error','coefficient_by_component','excluded_worst_error']}),flush=True)
out=dict(scope='Exposed source-shape calibration, not an observed gravitational fit or microscopic capture law',postulate='rho_extra=sum_j eta_j rho_b,j with nonnegative constant eta_j per component',hashes={str(v.relative_to(R)):hashlib.sha256(v.read_bytes()).hexdigest() for v in [Path(__file__),p,source]},cases=results)
(H/'results.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8',newline='\n')

"""Finite-time inward transport maps; conditional shared and leave-one-out fits."""
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parent
source=P/'stream-capture-results.json';data=json.loads(source.read_text())
B=float(np.log(100))
specs={'uniform':([0.],[B]),'partial':([0.,0.],[B,1.]),'supported':([0.,0.],[B,.5]),'retention_conditioned':([0.],[B])}
grid=np.linspace(0,1,401)

def mapped(row,branch,model,p,x):
    raw=np.asarray(row['case_cdfs'][branch]);xx=np.asarray(row['r_over_R'])
    base=np.interp(x,xx,raw)
    b=p[0]*(1-row['eta']) if model=='retention_conditioned' else p[0]
    s=np.exp(-b)
    inv=x/s
    if model=='supported':
        h=p[1];inv=np.where(x<=h,x,h+(x-h)/s)
    moved=np.interp(inv,xx,raw)
    return (1-p[1])*base+p[1]*moved if model=='partial' else moved

def target(row,x):return np.interp(x,row['r_over_R'],row['target_cdf'])

def fit(rows,branch,model):
    lo,hi=map(np.array,specs[model]); sols=[]
    def residual(p):return np.concatenate([mapped(r,branch,model,p,grid)-target(r,grid) for r in rows])
    for f in [.05,.35,.8]:
        result=least_squares(residual,lo+f*(hi-lo),bounds=(lo,hi),ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=400)
        sols.append((float(np.mean(residual(result.x)**2)),result.x))
    # Explicitly retain the zero-migration baseline as an allowed endpoint.
    sols.append((float(np.mean(residual(lo)**2)),lo))
    score,p=min(sols,key=lambda t:t[0])
    return p,score, bool(np.any((p-lo)<1e-4)|np.any((hi-p)<1e-4))

results=[];max_grid_drift=0
for geometry in ['companion_regular','standard_flat_FLRW']:
 for population in ['Chabrier','Salpeter']:
  rows=[r for r in data['rows'] if r['geometry']==geometry and r['population']==population]
  assert len(rows)==6
  for branch,label in [(-1,'distant'),(0,'near_1.1R')]:
   base_res=np.array([np.interp(grid,r['r_over_R'],r['case_cdfs'][branch])-target(r,grid) for r in rows])
   for model in specs:
    p,score,boundary=fit(rows,branch,model);per=[];loo_all=[]
    for i,r in enumerate(rows):
     lp,_,lb=fit(rows[:i]+rows[i+1:],branch,model)
     y=mapped(r,branch,model,p,grid);ly=mapped(r,branch,model,lp,grid);truth=target(r,grid)
     baseline=np.interp(grid,r['r_over_R'],r['case_cdfs'][branch])
     assert abs(y[0])<1e-10 and abs(y[-1]-1)<1e-10
     assert np.min(np.diff(y))>=-1e-10 and np.min(y-baseline)>=-1e-10
     assert np.max(np.abs(mapped(r,branch,model,[0.]*len(p),grid)-baseline))<1e-10
     fine=np.linspace(0,1,801);fr=mapped(r,branch,model,p,fine)-target(r,fine)
     drift=abs(float(np.sqrt(np.mean(fr*fr)))-float(np.sqrt(np.mean((y-truth)**2))))
     max_grid_drift=max(max_grid_drift,drift);assert drift<.002
     loo_all.extend((ly-truth).tolist())
     per.append(dict(Name=r['Name'],baseline_max_error=float(np.max(np.abs(baseline-truth))),shared_max_error=float(np.max(np.abs(y-truth))),loo_max_error=float(np.max(np.abs(ly-truth))),loo_parameters=lp.tolist(),loo_boundary=lb,shared_cdf=y.tolist(),loo_cdf=ly.tolist(),target_cdf=truth.tolist(),baseline_cdf=baseline.tolist()))
    results.append(dict(geometry=geometry,population=population,source=label,model=model,shared_parameters=p.tolist(),parameter_boundary=boundary,baseline_rms=float(np.sqrt(np.mean(base_res**2))),shared_rms=float(np.sqrt(score)),loo_rms=float(np.sqrt(np.mean(np.array(loo_all)**2))),rows=per))
    print(geometry,population,label,model,'RMS pp baseline/shared/LOO',*[round(100*v,3) for v in [results[-1]['baseline_rms'],results[-1]['shared_rms'],results[-1]['loo_rms']]],flush=True)
out=dict(scope='Normalized inventory-conserving inward migration; conditional halo-shape diagnostic, not full energy dynamics',input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),models_bounds=specs,grid=grid.tolist(),max_grid_rms_drift=max_grid_drift,results=results)
(P/'migration-variants-results.json').write_text(json.dumps(out,separators=(',',':'),allow_nan=False)+'\n',encoding='utf-8',newline='\n')
lines=['# Several inward-migration variations','',
'The exact one-third retention and existing capture constants remain unchanged. Four post-capture inward transport maps were executed using both distant illumination and the nearest tested source shell, for both geometries and both population proxies. These maps conserve deposited inventory; they do not yet derive binding, support or the energy released during settling.','',
'## What was tested','',
'Let x=r/R, where R=5 Re is the previous diagnostic region; b is a dimensionless migration rate times duration, and s=exp(-b). No universe age is assumed.','',
'| Variation | Where a captured element moves | Interpretation |','|---|---|---|',
'| Uniform | y=s*x | Every element relaxes inward |',
'| Partial | Fraction f moves to s*x; remainder stays | Only some captured energy is mobile |',
'| Supported | y=x below h; y=h+s*(x-h) above h | Settling slows toward a support radius |',
'| Retention-conditioned | y=exp[-b*(1-eta)]*x | Migration rate depends on the existing retention factor |','',
'These are known exponential-relaxation and distribution-transport constructions. Their application to companion deposits, the mobile fraction and stopping prescription are proposed phenomenology, not unique first-principles physics. The supported map solves dy/db=-(y-h) above h; the uniform map solves dy/db=-y. R-based scaling is a diagnostic assumption, not a gravitational force law. b=0 reproduces no migration.','',
'## Shared rule and transfer to an omitted galaxy','',
'Parameters are shared across six targets within each geometry/population/source branch. A separate leave-one-out exercise fits five and predicts the sixth without changing that fit. These are already-inspected halo targets, not fresh astronomical data. Errors below are RMS cumulative-fraction differences, in percentage points, averaged equally over radius and galaxies; they are not lensing or rotation residuals.','',
'| Geometry | Source | Variation | No migration | Shared fit | Omitted-galaxy transfer |','|---|---|---|---:|---:|---:|']
for r in results:
 if r['population']=='Chabrier':lines.append(f"| {r['geometry']} | {r['source']} | {r['model']} | {100*r['baseline_rms']:.2f} | {100*r['shared_rms']:.2f} | {100*r['loo_rms']:.2f} |")
lines+=['','All Salpeter cases, fitted parameters, boundary flags, individual targets and curves are retained in the JSON. The nearest-source branch assumes an isotropic source shell at 1.1 R; its existence has not been measured. It is not selected as the actual source distribution.','',
'## What this can establish','',
'Every strictly inward map increases enclosed deposit fraction at fixed radius. It can help a target requiring more central material, but cannot reduce an already excessive central fraction. Shared fitting therefore tests whether one movement rule can balance these incompatible demands. Closely fitting a target-derived halo shape still does not establish the source budget, its true shape, or correct projected lensing. Standard geometry is a comparison only.','',
'The finite-time maps track locations and retain total positive inventory. They do not conserve the full gravitational, kinetic and reservoir energy automatically: settling releases energy that requires an explicit receiving channel, and stopping needs stress or a binding mechanism. The support-radius candidate imposes where motion slows; it does not explain that support. Self-gravity and source histories have not been evolved.','',
'## Verification','',f"Zero-migration identity, endpoint inventory, monotonicity and inward cumulative ordering pass. Doubling evaluation-grid resolution changes fitted RMS by at most {max_grid_drift:.3g}. Three fitting starts and the explicit no-migration endpoint are compared. Bounds were declared before execution and were not expanded after fits. All six broader project goals remain open."]
summary=['## Findings in plain language','',
'Partial migration is the strongest of the four tested variations. For distant input under our retained geometry, a shared fit moves about 35.8% of the deposited inventory to 15.0% of its former radius, leaving the rest in place. Its cumulative-profile RMS falls from 26.74 to 16.80 percentage points in the shared fit and to 20.01 in omitted-galaxy transfer. These are substantial placement improvements, but not a solution: the improvement comes with worse central overconcentration in three other targets. This new fraction is fitted to normalized halo shapes and must not be substituted for the earlier 0.37% rotation-trained redistribution result.','',
'Uniform settling and retention-conditioned settling provide smaller improvements under the retained geometry. The supported model chooses h=0, so its added stopping-radius parameter collapses to ordinary uniform contraction rather than establishing nonzero support. No bounds were expanded.','',
'Under standard comparison geometry and distant input, none of the variants improves aggregate omitted-galaxy RMS: the no-migration value is 23.55 points and partial migration gives 24.32. The near-source partial case improves that branch from 26.23 to 24.32, still worse than the distant-input no-migration benchmark. Thus this is not a geometry-independent or source-independent solution.','',
'For the retained geometry and distant input, the following maximum cumulative errors show the tradeoff (percentage points):','',
'| Galaxy | No migration | Partial migration, omitted-galaxy prediction |','|---|---:|---:|']
chosen=next(r for r in results if r['geometry']=='companion_regular' and r['population']=='Chabrier' and r['source']=='distant' and r['model']=='partial')
for h in chosen['rows']:summary.append(f"| {h['Name']} | {100*h['baseline_max_error']:.2f} | {100*h['loo_max_error']:.2f} |")
summary+=['',
'The next mechanism question is what independently measured property controls the mobile fraction and where movement stops. Fitting a separate fraction to each desired halo would describe the answer rather than predict it. A gravity-dependent mobility or support law must be specified from source/matter inputs, then transferred without using the held-out halo shape. These results do not yet establish such a law.','']
pos=lines.index('## Verification');lines[pos:pos]=summary
(P/'migration-variants-report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
selected=[r for r in results if r['geometry']=='companion_regular' and r['population']=='Chabrier' and r['source']=='distant']
for i,ax in enumerate(axes.flat):
 h=selected[0]['rows'][i]
 ax.plot(grid,h['target_cdf'],'k',label='Fitted halo target');ax.plot(grid,h['baseline_cdf'],color='gray',label='No migration')
 for result in selected:ax.plot(grid,result['rows'][i]['loo_cdf'],label=result['model'])
 ax.set(title=h['Name'],xlabel='Radius / diagnostic region',ylabel='Enclosed fraction',xlim=(0,1),ylim=(0,1));ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=7)
fig.suptitle('Inward migration: omitted-galaxy predictions, distant input, retained geometry')
fig.savefig(P/'migration-variants.png',dpi=150)

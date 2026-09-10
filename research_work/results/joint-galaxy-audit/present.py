"""Create a figure and a self-contained, filterable observed/predicted table."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

H=Path(__file__).resolve().parent
load=lambda f:json.loads((H/f).read_text())
red=load('redshift-predictions.json'); sp=load('galaxy-rotation-predictions.json'); mw=load('milky-way-predictions.json'); repair=load('repair-predictions.json')
plt.rcParams.update({'font.size':10,'svg.hashsalt':'joint-galaxy-audit-2026-09-09'})
fig,ax=plt.subplots(2,3,figsize=(16,9),layout='constrained')
for split,marker in [('train','.'),('validation','x'),('test','o')]:
    r=[x for x in red if x['model']=='fixed_alpha' and x['split']==split]
    ax[0,0].scatter([x['D_Mpc'] for x in r],[x['residual_kms'] for x in r],marker=marker,label=split,s=20,alpha=.75)
ax[0,0].axhline(0,color='black',lw=.8); ax[0,0].set(title='164 groups: fixed redshift law',xlabel='Adopted distance (Mpc)',ylabel='c × (predicted z − observed z), km/s'); ax[0,0].legend(fontsize=8)
for model,label in [('baryons','Ordinary matter'),('power','Empirical extra gravity')]:
    r=[x for x in sp if x['model']==model and x['split']=='test']
    ax[0,1].scatter([x['observed_kms'] for x in r],[x['predicted_kms'] for x in r],s=6,alpha=.35,label=label)
ax[0,1].plot([0,350],[0,350],'k--',lw=.8); ax[0,1].set(title='31 reserved galaxies: 642 points',xlabel='Observed rotation (km/s)',ylabel='Predicted rotation (km/s)'); ax[0,1].legend(fontsize=8)
obs=[x for x in mw if x['baryons']=='I' and x['model']=='baryons' and x['observable']=='vc_kms']
source=json.loads((H.parent/'milky-way-capture/inputs.json').read_text())
ax[0,2].errorbar([x['R_kpc'] for x in obs],[x['observed'] for x in obs],yerr=[[x['err_minus_kms'] for x in source['eilers']['rows']],[x['err_plus_kms'] for x in source['eilers']['rows']]],fmt='.',c='black',label='Published speeds (statistical errors)')
for variant in ['I','II']:
    r=[x for x in repair if x['baryons']==variant and x['observable']=='vc_kms']
    ax[0,2].plot([x['R_kpc'] for x in r],[x['predicted'] for x in r],label=f'Capture + adjusted stars {variant}')
ax[0,2].axvline(15,c='gray',ls=':'); ax[0,2].set(title='Milky Way: outer radii withheld',xlabel='Radius (kpc)',ylabel='Circular speed (km/s)'); ax[0,2].legend(fontsize=8)
obs=[x for x in mw if x['baryons']=='I' and x['model']=='baryons' and x['observable']!='vc_kms']
ax[1,0].errorbar([x['R_kpc'] for x in obs],[x['observed'] for x in obs],yerr=[x['Kz_error'] for x in source['bovy']['rows']],fmt='.',c='black',alpha=.65,label='Published total-force inference')
for variant in ['I','II']:
    r=sorted([x for x in repair if x['baryons']==variant and x['observable']=='Kz_over_2piG'],key=lambda x:x['R_kpc'])
    ax[1,0].plot([x['R_kpc'] for x in r],[x['predicted'] for x in r],label=f'Joint diagnostic {variant}')
ax[1,0].set(title='Vertical pull at 1.1 kpc above disk',xlabel='Radius (kpc)',ylabel='|Kz|/(2πG), solar masses/pc²'); ax[1,0].legend(fontsize=8)
for p in load('capture-profiles.json'):
    if p['kind']!='well': continue
    R=np.array(p['R_kpc']); density=np.array(p['rho_Msun_kpc3']); use=(R>=4)&(R<=30)
    ax[1,1].plot(R[use],density[use]/np.interp(5,R,density),label=f'Initial capture model {p["baryons"]}')
ax[1,1].set(title='Fitted capture: density decreases outward',xlabel='Radius (kpc)',ylabel='Deposited density / density at 5 kpc',yscale='log'); ax[1,1].legend(fontsize=8)
summ=load('results.json'); re=load('repair-results.json')
labels=['Fixed stars I','Adjusted stars I','Fixed stars II','Adjusted stars II']
values=[summ['milky_way']['I']['scores']['capture_well']['vertical_diagonal_score'],re['I']['vertical_all_diagonal_score'],summ['milky_way']['II']['scores']['capture_well']['vertical_diagonal_score'],re['II']['vertical_all_diagonal_score']]
ax[1,2].bar(labels,values,color=['#8097b2','#2b7774','#8097b2','#2b7774']); ax[1,2].tick_params(axis='x',labelrotation=20); ax[1,2].set(title='Stellar mass uncertainty changes the result',ylabel='Vertical squared standardized residual sum'); ax[1,2].text(.03,.93,'Exploratory; no full covariance or blind test',transform=ax[1,2].transAxes,fontsize=9)
fig.savefig(H/'summary.png',dpi=150); fig.savefig(H/'summary.svg',metadata={'Date':None})
path=H/'summary.svg'; path.write_text('\n'.join(x.rstrip() for x in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8',newline='\n')

records=[]
for x in red: records.append(dict(sample='Redshift',object=x['group_pgc'],model=x['model'],split=x['split'],position=x['D_Mpc'],unit='redshift z; position Mpc',observed=x['observed_z'],predicted=x['predicted_z']))
for x in sp: records.append(dict(sample='Other galaxies',object=x['galaxy'],model=x['model'],split=x['split'],position=x['R_kpc'],unit='km/s; radius kpc',observed=x['observed_kms'],predicted=x['predicted_kms']))
for x in mw: records.append(dict(sample='Milky Way rotation' if x['observable']=='vc_kms' else 'Milky Way vertical pull',object='Stellar model '+x['baryons'],model=x['model'],split=x['split'],position=x['R_kpc'],unit=x['observable']+'; radius kpc',observed=x['observed'],predicted=x['predicted']))
for x in repair: records.append(dict(sample='Milky Way rotation' if x['observable']=='vc_kms' else 'Milky Way vertical pull',object='Stellar model '+x['baryons'],model='capture_adjusted_stars_POSTHOC',split=x['split'],position=x['R_kpc'],unit=x['observable']+'; radius kpc',observed=x['observed'],predicted=x['predicted']))
for x in load('bulge-inputs.json')['fields']:
    for key,sample in [('mean_heliocentric_RV_kms','Bulge mean motion'),('dispersion_kms','Bulge velocity spread')]:
        records.append(dict(sample=sample,object=x['survey']+' '+x['field'],model='not_derived',split='published_summary',position=x['b_label_deg'],unit='km/s; position is signed Galactic latitude in degrees',observed=x[key],predicted=None))
for r in records:
    if r['sample']=='Milky Way rotation': r['unit']='km/s; radius kpc'
    if r['sample']=='Milky Way vertical pull': r['unit']='|Kz|/(2πG), solar masses/pc²; radius kpc'
html='''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Observed and predicted galaxy data</title>
<style>body{font:16px system-ui;margin:2rem auto;max-width:1200px;padding:0 1rem;color:#20313b;background:#f5f7f7}h1{font-size:1.7rem}select,input{font:inherit;padding:.5rem;margin:.4rem}label{display:inline-block}table{width:100%;border-collapse:collapse;background:white}td,th{padding:.55rem;border-bottom:1px solid #ddd;text-align:right}td:first-child,th:first-child{text-align:left}thead{position:sticky;top:0;background:#e1eeee}p{line-height:1.5}.note{background:#e1eeee;padding:1rem}button{font:inherit;padding:.5rem}#count{font-weight:600}</style>
<h1>Observed and predicted galaxy data</h1><p>All rows from the joint audit. “Train”, “test” and “reserved” preserve analysis splits; these data have already been inspected and are not a new blind test. The formulas do not match every observation.</p>
<p class="note">Redshift is a dimensionless wavelength ratio. Rotation is speed around the galaxy. Vertical pull is gravity toward the disk from above or below, expressed as |Kz|/(2πG). Milky Way stellar profiles and published vertical-force estimates retain model assumptions. The adjusted-stars model is a post-result repair experiment. Capture mass is fitted, not derived from the available photon energy.</p>
<label>Sample <select id="sample"></select></label><label>Model <select id="model"></select></label><label>Object <select id="object"></select></label><label>Subset <select id="split"></select></label>
<p id="count" aria-live="polite"></p><div style="overflow:auto"><table><thead><tr><th>Object</th><th>Model</th><th>Subset</th><th>Distance/radius</th><th>Observed</th><th>Predicted</th><th>Predicted − observed</th></tr></thead><tbody id="body"></tbody></table></div>
<p>Source and formula provenance: see report.md and protocol.json in this delivery. No external network requests are made by this file.</p><script>
const records=DATA;
const labels={fixed_alpha:'Original constant redshift rate',refit_alpha:'Refitted constant redshift rate',quadratic_history:'Modified travel-history rate',baryons:'Ordinary stars and gas',power:'Empirical extra-gravity rule',capture_uniform:'Uniform capture',capture_well:'Centre-weighted capture',sparc_transfer:'Rule calibrated on other galaxies',capture_adjusted_stars_POSTHOC:'Capture with adjusted stellar mass (exploratory)',not_derived:'No dynamical prediction yet',vertical_provisional:'Provisional vertical comparison',published_summary:'Published field summary'};
const ids=['sample','model','object','split'];const el=Object.fromEntries(ids.map(k=>[k,document.getElementById(k)]));
function options(k,values,all){const old=el[k].value;el[k].replaceChildren();for(const v of (all?['All',...values]:values)){const o=document.createElement('option');o.value=v;o.textContent=labels[v]||v;el[k].append(o)}if([...el[k].options].some(o=>o.value===old))el[k].value=old;}
function refreshOptions(){let rows=records.filter(r=>r.sample===el.sample.value);for(const k of ['model','object','split'])options(k,[...new Set(rows.map(r=>r[k]))].sort(),true);draw()}
function draw(){const rows=records.filter(r=>r.sample===el.sample.value&&['model','object','split'].every(k=>el[k].value==='All'||r[k]===el[k].value));document.getElementById('count').textContent=rows.length+' rows (models can repeat observations). Units: '+(rows[0]?.unit||'—');const body=document.getElementById('body');body.replaceChildren();for(const r of rows){const tr=document.createElement('tr');for(const v of [r.object,labels[r.model]||r.model,labels[r.split]||r.split,r.position,r.observed,r.predicted,r.predicted===null?null:r.predicted-r.observed]){const td=document.createElement('td');td.textContent=v===null?'Not derived':typeof v==='number'?v.toPrecision(7):v;tr.append(td)}body.append(tr)}}
options('sample',[...new Set(records.map(r=>r.sample))],false);el.sample.addEventListener('change',refreshOptions);for(const k of ['model','object','split'])el[k].addEventListener('change',draw);refreshOptions();</script></html>'''
html=html.replace('const records=DATA;','const records='+json.dumps(records,ensure_ascii=True).replace('</','<\\/')+';').replace('Distance/radius</th>','Position (see units)</th>')
(H/'comparison.html').write_text(html,encoding='utf-8',newline='\n')
print('Wrote figure and explorer with',len(records),'rows, including bulge observations whose predictions remain undetermined.')

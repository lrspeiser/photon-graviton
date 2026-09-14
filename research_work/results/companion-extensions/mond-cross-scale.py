"""Frozen SPARC-guided rule transferred to all archived Milky Way sensitivities.

python mond-cross-scale.py [--output-dir DIR] [--canonical]
Regenerates into a fresh directory and compares with the archived
mond-cross-scale-results.json; only --canonical overwrites the archive. The
archived file was produced by the pre-entry-point version at 1967ff8.
"""
from pathlib import Path
import json,hashlib,sys
import numpy as np
P=Path(__file__).resolve().parent;OLD=P.parent/'isotropic-galaxy-transfer'
sys.path.insert(0,str(P))
import evidence_io  # noqa: E402
paths=[P/'mond-inventory-results.json',OLD/'milky-way-current-results.json']

def compute():
 rule=json.loads(paths[0].read_text());source=json.loads(paths[1].read_text());f=rule['f'];a0=rule['a0_SI']*3.085677581491367e19/1e6;G=4.30091727003628e-6
 out=dict(scope='Frozen f and a0 transferred without fitting any Milky Way speed; all 18 archived mass/size/luminosity scenarios retained. Cluster and redshift integration specified separately, no new cluster observational fit.',f=f,a0_SI=rule['a0_SI'],runs=[])
 for run in source['runs']:
  rows=[r for r in run['rows'] if r['model']=='exact_third'];R=np.array([r['R_kpc'] for r in rows]);v=np.array([r['predicted_kms'] for r in rows]);vb=np.array([r['ordinary_kms'] for r in rows]);y=np.array([r['observed_kms'] for r in rows])
  gb=vb**2/R;gc=2*a0*gb/(np.sqrt(gb*gb+4*a0*gb)+gb);target=R*R*gc/G
  total=run['total_deposit_Msun'];end=np.minimum(np.maximum.accumulate(target),total);ref=R*(v*v-vb*vb)/G
  newmass=(1-f)*ref+f*end;pred=np.sqrt(vb*vb+G*newmass/R)
  assert np.all(np.diff(end)>=0) and max(end)<=total
  rawmond=np.sqrt(vb*vb+R*gc)
  scores={}
  for name,values in [('reference',v),('MOND_raw',rawmond),('redistributed',pred)]:
   scores[name]=dict(RMSE_kms=float(np.sqrt(np.mean((values-y)**2))),bias_kms=float(np.mean(values-y)))
  out['runs'].append(dict(baryons=run['baryons'],Rd_kpc=run['Rd_kpc'],luminosity_proxy_factor=run['luminosity_proxy_factor'],total_inventory_Msun=total,negative_target_intervals=int(sum(np.diff(target)<0)),capped_radii=int(sum(target>total)),scores=scores,rows=[dict(R_kpc=float(r),observed_kms=float(o),reference_kms=float(b),redistributed_kms=float(n),MOND_raw_kms=float(m)) for r,o,b,n,m in zip(R,y,v,pred,rawmond)]))
 out['fiducial']=[r for r in out['runs'] if r['Rd_kpc']==2.6 and r['luminosity_proxy_factor']==1]
 out['improved_sensitivity_cases']=sum(r['scores']['redistributed']['RMSE_kms']<r['scores']['reference']['RMSE_kms'] for r in out['runs'])
 out['source_sha256']={str(p.relative_to(P.parents[2])).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths+[Path(__file__)]}
 print([(r['baryons'],r['scores'],r['capped_radii']) for r in out['fiducial']]);print('improved',out['improved_sensitivity_cases'],'of',len(out['runs']))
 return out

def main():
 args=evidence_io.parse(__doc__)
 out=compute()
 return evidence_io.finish(args,'mond-cross-scale',json.dumps(out,indent=2)+'\n',P/'mond-cross-scale-results.json',
                           ignore=['/source_sha256/research_work/results/companion-extensions/mond-cross-scale.py'])

if __name__=='__main__':
 raise SystemExit(main())

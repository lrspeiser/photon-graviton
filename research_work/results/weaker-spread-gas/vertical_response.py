#!/usr/bin/env python3
"""Post-fit conditional gas thickness bound; no optimization or new selection.
For a spherical companion, vertical harmonic curvature is gchi(R)/R.
At fixed gas dispersion and column, an external harmonic potential gives
h=cg/nu. With unchanged nonnegative baryonic vertical curvature, the actual
thickness ratio lies between 1 and sqrt(gchi_old/gchi_new) for pure dilation.
Gas self-gravity, multiphase pressure, feedback and radial rearrangement are
not calculated. This bound is not a measured gas thickness or SFR.
"""
import argparse, json, sys
from pathlib import Path
import numpy as np

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--results',type=Path,required=True);ap.add_argument('--jr1',type=Path,required=True);a=ap.parse_args();out=a.results
 sys.path.insert(0,str(a.jr1/'repository/research_work/results/joint-response-iteration'))
 import population_followup as P
 E=P.Experiment();p=json.loads((a.jr1/'evidence/R10_spheroid_stellar_population.json').read_text())['parameters'];u=10**p['logusph'];fs=E.Msph/E.Ms;ufs=u*fs/(1+(u-1)*fs)
 A,rc,rt=P.source(E.Ms,E.Mg,E.Re,p,False,fs);q=p['q']+(p['q_sph']-p['q'])*ufs
 saved=json.loads((out/'summary.json').read_text());sel=json.loads((out/'all-model-predictions.json').read_text())[saved['selected']]
 inv=json.loads((out/'individual-diagnostics.json').read_text())['rows'];rows=[]
 for i,g in enumerate(E.galaxies):
  R=g['Re'];before=float(P.J.companion_force(R,A[i],rc[i],rt[i],q[i]));scenarios={}
  for label,s,lam in [('shared',sel['galaxies'][i]['strength'],sel['galaxies'][i]['spread']),('individual_spread',1.,inv[i]['models']['spreading']['spread'])]:
   after=float(P.J.companion_force(R,s*A[i]/lam,rc[i]*lam,rt[i]*lam,q[i]));ratio=after/before
   scenarios[label]=dict(companion_vertical_curvature_ratio=ratio,maximum_external_field_thickness_ratio=ratio**-.5,
    minimum_midplane_density_ratio_at_fixed_column=ratio**.5,maximum_freefall_time_ratio_at_fixed_G_and_column=ratio**-.25)
  rows.append(dict(name=g['name'],original_over=inv[i]['over'],gas_fraction=inv[i]['fgas'],scenarios=scenarios))
 groups={}
 for label in ('shared','individual_spread'):
  group=[r['scenarios'][label] for r in rows if r['original_over']]
  groups[label]=dict(n=len(group),median_maximum_thickness_ratio=float(np.median([r['maximum_external_field_thickness_ratio'] for r in group])),median_maximum_freefall_ratio=float(np.median([r['maximum_freefall_time_ratio_at_fixed_G_and_column'] for r in group])))
 data=dict(status='Post-primary analytic diagnostic; no changes to fits',scope=__doc__,rows=rows,overpredicted_outlier_groups=groups)
 (out/'vertical-response-bounds.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');print(json.dumps(groups,indent=2))
if __name__=='__main__':main()

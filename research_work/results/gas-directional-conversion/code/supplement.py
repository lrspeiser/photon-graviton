"""Declared post-selection attribution only; does not revise selected model."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from core import AxisModel
from data_screen import prepare,score,predict_one


def main():
    p=argparse.ArgumentParser();p.add_argument('--jr1',type=Path,required=True);p.add_argument('--results',type=Path,required=True);a=p.parse_args()
    target=a.results/'supplement.json'
    if target.exists():raise FileExistsError(target)
    gs,profiles=prepare(a.jr1,profiles=json.loads((a.results/'gas-profiles.json').read_text()))
    mix=json.loads((a.results/'forward_linear_reverse_quadratic.json').read_text())
    full=score(gs,mix['alpha'],mix['beta'],True)
    (a.results/'mixed-full-postselection.json').write_text(json.dumps(full,indent=2)+'\n')
    mono=score(gs,mix['alpha'],mix['beta'],True,monopole=True)
    # Preserved source-match subgroup, using ONLY gas-force mismatch, not Vobs.
    group=[]
    for name in ['selected-full.json','mixed-full-postselection.json','null-full.json']:
      rec=json.loads((a.results/name).read_text());lookup={p['name']:p for p in profiles}
      for good in [True,False]:
        rows=[r for r in rec['rows'] if (lookup[r['name']]['gas_force_relative_RMS']<=.2)==good]
        group.append(dict(model=name,gas_force_reconstruction_within_20pct=good,galaxies=len(rows),
           mean_RMS_kms=float(np.mean([r['RMS_kms'] for r in rows])),
           median_fractional_RMS=float(np.median([r['fractional_RMS'] for r in rows]))))
    # Reverse the order of the same two layers, now with purely photon input.
    layers={}
    for name,parts in {'diffuse_then_dense':[(2.5,2e6),(.1,5e7)],'dense_then_diffuse':[(.1,5e7),(2.5,2e6)]}.items():
      state=np.array([1.,0.]);frac=0.;error=0.
      for width,rho in parts:
        ap=.3*rho/1e6;am=.3*rho/1e6*rho/1e7
        state=expm(np.array([[-ap,am],[ap,-am]])*width)@state
        eq=ap/(ap+am);frac=eq+(frac-eq)*np.exp(-(ap+am)*width);error=max(error,abs(frac-state[1]))
      layers[name]=dict(photon_final=float(state[0]),companion_final=float(state[1]),energy_error=float(abs(state.sum()-1)),closed_form_error=float(error))
    # Store reproducible high-resolution plotting arrays (not observed images).
    plotrows=[]
    for kind,alpha,beta in [('disk',.3,.3),('ring',3.,0.)]:
      m=AxisModel(1,1,20,2,nr=1024,nmu=128,lmax=48).solve(1e9,2,alpha,beta,kind=kind)
      radii=np.geomspace(.3,12,160);d=m.diagnostic(radii)
      base=AxisModel(1,1,20,2,nr=1024,nmu=128,lmax=48).solve(0.,2,0,0)
      d.update(kind=kind,alpha=alpha,beta=beta,lens_ratio={s:m.bend(2,s)/base.bend(2,s) for s in ['face','edge','edge_vertical']})
      plotrows.append(d)
      np.savez_compressed(a.results/f'toy-{kind}-map.npz',r=m.r,mu=m.mu,companion_ratio=2*m.f)
    result=dict(scope='post-selection diagnostic comparison; selection unchanged; no new optimization',
       mixed_groups=full['groups'],mixed_monopole_groups=mono['groups'],gas_reconstruction_subgroups=group,
       pure_photon_layer_order=layers,refined_toy_plot_rows=plotrows)
    target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='refined_toy_plot_rows'}),flush=True)
if __name__=='__main__':main()
